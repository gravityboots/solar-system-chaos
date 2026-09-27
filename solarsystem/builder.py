"""Build a :class:`rebound.Simulation` (plus REBOUNDx physics) from a Config.

The returned :class:`SimContext` holds *strong references* to the REBOUNDx
``Extras`` object and any force objects.  This is essential: if they are garbage
collected, the extra forces silently detach from the simulation.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Any

import rebound
import reboundx

from . import constants
from .config import Config

DEG = math.pi / 180.0


@dataclass
class SimContext:
    sim: rebound.Simulation
    config: Config
    body_names: list[str]                 # ["Sun", planet, ...] in particle order
    planet_names: list[str]               # planets only
    n_bodies: int                         # number of *real* particles (incl. Sun)
    rebx: Any = None                      # reboundx.Extras (keep alive!)
    gr_force: Any = None
    j2_force: Any = None
    gr_mode: str = "none"
    masses: dict[str, float] = field(default_factory=dict)
    frozen_states: dict[int, tuple] = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# Body-table resolution
# --------------------------------------------------------------------------- #
def resolve_body_table(cfg: Config) -> dict[str, dict]:
    """Return {name: {mass, a, e, inc, Omega, pomega, l, ...}} after overrides."""
    src = cfg.bodies.source
    if src == "builtin":
        table = {k: dict(v) for k, v in constants.DEFAULT_BODIES.items()}
    elif src == "file":
        from .horizons import load_body_table
        if not cfg.bodies.file:
            raise ValueError("bodies.source=file requires bodies.file to be set")
        table = load_body_table(cfg.bodies.file)
    elif src == "horizons":
        from .horizons import fetch_body_table
        table = fetch_body_table(cfg.bodies.include, epoch=cfg.bodies.epoch)
    else:
        raise ValueError(f"Unknown bodies.source: {src!r}")

    # Apply per-body overrides (mass scaling/absolute, element edits, include).
    for name, ov in cfg.bodies.overrides.items():
        if name not in table:
            continue
        entry = table[name]
        if ov.mass is not None:
            entry["mass"] = ov.mass
        elif ov.mass_scale != 1.0:
            entry["mass"] = entry["mass"] * ov.mass_scale
        for el in ("a", "e", "inc", "Omega", "pomega", "l"):
            val = getattr(ov, el)
            if val is not None:
                entry[el] = val
    return table


def effective_planets(cfg: Config, table: dict[str, dict]) -> list[str]:
    """Planets to include, in config order, honouring per-body include flags."""
    out = []
    for name in cfg.bodies.include:
        if name not in table:
            continue
        ov = cfg.bodies.overrides.get(name)
        if ov is not None and not ov.include:
            continue
        out.append(name)
    return out


# --------------------------------------------------------------------------- #
# Simulation construction
# --------------------------------------------------------------------------- #
def build_simulation(cfg: Config) -> SimContext:
    table = resolve_body_table(cfg)
    planets = effective_planets(cfg, table)

    sim = rebound.Simulation()
    sim.units = constants.UNITS
    sim.add(m=constants.SUN_MASS, hash="Sun")

    masses = {"Sun": constants.SUN_MASS}
    for name in planets:
        b = table[name]
        sim.add(
            m=b["mass"],
            a=b["a"], e=b["e"],
            inc=b["inc"] * DEG,
            Omega=b["Omega"] * DEG,
            pomega=b["pomega"] * DEG,
            l=b["l"] * DEG,
            primary=sim.particles[0],  # heliocentric J2000 elements
            hash=name,
        )
        masses[name] = b["mass"]

    sim.move_to_com()

    body_names = ["Sun"] + planets
    ctx = SimContext(
        sim=sim, config=cfg, body_names=body_names, planet_names=planets,
        n_bodies=len(body_names), masses=masses,
    )

    _configure_integrator(ctx)
    _setup_physics(ctx)          # before init_megno so forces apply to variationals
    _setup_megno(ctx)
    _setup_freeze(ctx)
    return ctx


def _configure_integrator(ctx: SimContext) -> None:
    sim, cfg = ctx.sim, ctx.config
    name = cfg.integrator.name.lower()
    sim.dt = constants.days_to_code_time(cfg.integrator.dt_days)
    sim.exact_finish_time = cfg.integrator.exact_finish_time

    if name in ("whfast", "whckl"):
        sim.integrator = "whfast"
        w = cfg.integrator.whfast
        sim.ri_whfast.safe_mode = int(w.safe_mode)

        corrector = int(w.corrector)
        kernel = w.kernel
        # MEGNO uses variational particles, which are incompatible with the lazy
        # (WHCKL) kernel and ~30x slower with a symplectic corrector for no gain
        # in the MEGNO value.  So MEGNO runs use the standard kernel and no
        # corrector; clean runs keep the high-precision WHCKL settings.
        if cfg.chaos.megno.enabled:
            if kernel not in ("default", "standard") or corrector != 0:
                warnings.warn(
                    "MEGNO variational particles require the standard kernel and "
                    "run ~30x faster without a symplectic corrector (the MEGNO "
                    "value is unchanged); forcing kernel='default', corrector=0 "
                    "for this run.", stacklevel=2)
            kernel, corrector = "default", 0

        sim.ri_whfast.corrector = corrector
        try:
            sim.ri_whfast.corrector2 = int(w.corrector2)
        except Exception:
            pass
        try:
            sim.ri_whfast.kernel = kernel
        except Exception:
            pass
    else:
        sim.integrator = name  # ias15 / mercurius / trace / saba / ...


def _setup_physics(ctx: SimContext) -> None:
    cfg = ctx.config
    gr_mode = cfg.physics.gr.mode.lower()
    need_rebx = gr_mode != "none" or cfg.physics.j2.enabled
    if not need_rebx:
        ctx.gr_mode = "none"
        return

    rebx = reboundx.Extras(ctx.sim)
    ctx.rebx = rebx
    ctx.gr_mode = gr_mode

    if gr_mode != "none":
        force = rebx.load_force(gr_mode)
        c = cfg.physics.gr.c if cfg.physics.gr.c is not None else constants.SPEED_OF_LIGHT
        force.params["c"] = c
        rebx.add_force(force)
        ctx.gr_force = force
        if gr_mode in ("gr", "gr_potential"):
            for src in cfg.physics.gr.sources:
                try:
                    ctx.sim.particles[src].params["gr_source"] = 1
                except Exception:
                    pass

    if cfg.physics.j2.enabled:
        gh = rebx.load_force("gravitational_harmonics")
        rebx.add_force(gh)
        sun = ctx.sim.particles["Sun"]
        sun.params["J2"] = cfg.physics.j2.J2
        sun.params["R_eq"] = cfg.physics.j2.R_sun_AU
        ctx.j2_force = gh


def _setup_megno(ctx: SimContext) -> None:
    if ctx.config.chaos.megno.enabled:
        # init_megno must be called before any integration has happened.
        ctx.sim.init_megno(seed=ctx.config.chaos.megno.seed)


def _setup_freeze(ctx: SimContext) -> None:
    frozen = list(ctx.config.freeze.bodies)
    if not frozen:
        return
    sim = ctx.sim
    idx = {name: i for i, name in enumerate(ctx.body_names)}
    for name in frozen:
        if name not in idx:
            continue
        p = sim.particles[idx[name]]
        ctx.frozen_states[idx[name]] = (p.x, p.y, p.z, p.vx, p.vy, p.vz)

    states = ctx.frozen_states

    def heartbeat(reb_sim):
        s = reb_sim.contents
        for i, (x, y, z, vx, vy, vz) in states.items():
            pi = s.particles[i]
            pi.x, pi.y, pi.z = x, y, z
            pi.vx, pi.vy, pi.vz = vx, vy, vz
        s.ri_whfast.recalculate_coordinates_this_timestep = 1

    sim.heartbeat = heartbeat
