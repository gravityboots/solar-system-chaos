"""Integration driver: advance the simulation and record diagnostics.

The integration is done in chunks so we can sample energy, angular momentum,
MEGNO and orbital elements on a regular grid (uniform in time by default, which
the secular FFT needs).  With WHFast safe_mode=0 we ``synchronize()`` once per
chunk before reading state, then continue — cheap, and it keeps the trajectory
symplectic.  A bit-reproducible SimulationArchive is written alongside.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from tqdm import tqdm

from . import constants, diagnostics
from .builder import SimContext


@dataclass
class RunResult:
    body_names: list[str]
    planet_names: list[str]
    times_years: np.ndarray
    energy: np.ndarray
    dE_rel: np.ndarray
    angmom: np.ndarray
    dL_rel: np.ndarray
    amd: np.ndarray
    megno: np.ndarray | None
    lyap_time_years: np.ndarray | None
    elements: dict[str, dict[str, np.ndarray]]      # name -> field -> series
    positions: dict[str, np.ndarray]                # name -> (n, 3) AU
    meta: dict[str, Any] = field(default_factory=dict)

    def save_npz(self, path: str | Path) -> None:
        flat: dict[str, Any] = dict(
            times_years=self.times_years, energy=self.energy, dE_rel=self.dE_rel,
            angmom=self.angmom, dL_rel=self.dL_rel, amd=self.amd,
        )
        if self.megno is not None:
            flat["megno"] = self.megno
        if self.lyap_time_years is not None:
            flat["lyap_time_years"] = self.lyap_time_years
        for name, series in self.elements.items():
            for fld, arr in series.items():
                flat[f"elem__{name}__{fld}"] = arr
        for name, arr in self.positions.items():
            flat[f"pos__{name}"] = arr
        np.savez_compressed(path, **flat)


def _step_targets(total_steps: int, n_outputs: int, spacing: str) -> np.ndarray:
    if spacing == "log":
        first = max(1.0, total_steps / (n_outputs * 50))
        raw = np.logspace(math.log10(first), math.log10(total_steps), n_outputs)
    else:
        raw = np.linspace(total_steps / n_outputs, total_steps, n_outputs)
    targets = np.unique(np.round(raw).astype(np.int64))
    return targets[targets > 0]


def run(ctx: SimContext, outdir: str | Path, progress: bool = True) -> RunResult:
    sim, cfg = ctx.sim, ctx.config
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dt = sim.dt
    t_end_code = constants.years_to_code_time(cfg.time.t_end_years)
    total_steps = max(int(round(t_end_code / dt)), cfg.time.n_outputs)
    targets = _step_targets(total_steps, cfg.time.n_outputs, cfg.time.spacing)

    # Set up the automatic SimulationArchive (snapshots during integrate()).
    archive_path = outdir / "archive.bin"
    if cfg.archive.enabled:
        iv_years = cfg.archive.interval_years or (cfg.time.t_end_years / cfg.time.n_outputs)
        sim.save_to_file(str(archive_path), interval=constants.years_to_code_time(iv_years),
                         delete_file=True)

    megno_on = cfg.chaos.megno.enabled
    amd_on = cfg.diagnostics.compute_amd

    # Buffers (start with the t=0 sample).
    t_list = [constants.code_time_to_years(sim.t)]
    e0 = diagnostics.total_energy(ctx)
    l0 = diagnostics.angular_momentum_magnitude(ctx)
    energy_list = [e0]
    angmom_list = [l0]
    amd_list = [diagnostics.angular_momentum_deficit(ctx) if amd_on else math.nan]
    megno_list = [sim.megno() if megno_on else math.nan]
    lyap_list = [diagnostics.lyapunov_time_years(ctx) if megno_on else math.nan]
    elem_buf: dict[str, dict[str, list]] = {
        n: {f: [] for f in diagnostics.ELEMENT_FIELDS} for n in ctx.planet_names}
    pos_buf: dict[str, list] = {n: [] for n in ctx.planet_names}
    _record_elements_positions(ctx, elem_buf, pos_buf)

    status = "completed"
    iterator = tqdm(targets, desc=f"Integrating {cfg.run.name}", unit="chk",
                    disable=not progress, mininterval=5.0)
    for target_step in iterator:
        try:
            sim.integrate(float(target_step) * dt)
            sim.synchronize()
        except Exception as exc:  # Escape / Collision / numerical blow-up
            status = f"stopped: {type(exc).__name__}: {exc}"
            break

        energy = diagnostics.total_energy(ctx)
        if not math.isfinite(energy):
            status = "stopped: non-finite energy (close encounter?)"
            break

        t_list.append(constants.code_time_to_years(sim.t))
        energy_list.append(energy)
        angmom_list.append(diagnostics.angular_momentum_magnitude(ctx))
        amd_list.append(diagnostics.angular_momentum_deficit(ctx) if amd_on else math.nan)
        if megno_on:
            megno_list.append(sim.megno())
            lyap_list.append(diagnostics.lyapunov_time_years(ctx))
        _record_elements_positions(ctx, elem_buf, pos_buf)

    # Assemble arrays.
    times = np.asarray(t_list)
    energy = np.asarray(energy_list)
    angmom = np.asarray(angmom_list)
    dE_rel = np.abs(energy - e0) / abs(e0) if e0 != 0 else np.abs(energy)
    dL_rel = np.abs(angmom - l0) / l0 if l0 != 0 else np.abs(angmom)
    amd = np.asarray(amd_list)
    megno = np.asarray(megno_list) if megno_on else None
    lyap = np.asarray(lyap_list) if megno_on else None
    elements = {n: {f: np.asarray(v) for f, v in series.items()}
                for n, series in elem_buf.items()}
    positions = {n: np.asarray(v) for n, v in pos_buf.items()}

    meta = dict(
        run_name=cfg.run.name, status=status,
        n_samples=int(times.size), n_bodies=ctx.n_bodies,
        dt_days=cfg.integrator.dt_days, integrator=cfg.integrator.name,
        t_end_years=cfg.time.t_end_years, t_reached_years=float(times[-1]),
        total_steps=int(total_steps), gr_mode=ctx.gr_mode,
        j2=cfg.physics.j2.enabled, megno=megno_on,
        archive_path=str(archive_path) if cfg.archive.enabled else None,
        E0=float(e0), L0=float(l0),
        megno_final=float(megno[-1]) if megno is not None else None,
        lyap_time_years_final=float(lyap[-1]) if lyap is not None else None,
        masses={n: ctx.masses[n] for n in ctx.body_names},
    )

    return RunResult(
        body_names=ctx.body_names, planet_names=ctx.planet_names,
        times_years=times, energy=energy, dE_rel=dE_rel, angmom=angmom,
        dL_rel=dL_rel, amd=amd, megno=megno, lyap_time_years=lyap,
        elements=elements, positions=positions, meta=meta,
    )


def _record_elements_positions(ctx, elem_buf, pos_buf) -> None:
    els = diagnostics.orbital_elements(ctx)
    pos = diagnostics.positions(ctx)
    for name in ctx.planet_names:
        for f in diagnostics.ELEMENT_FIELDS:
            elem_buf[name][f].append(els[name][f])
        pos_buf[name].append(pos[name])
