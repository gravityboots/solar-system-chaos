"""Physical diagnostics read from a (synchronized) simulation state.

These are the "standard statistical outputs": conserved quantities (energy,
angular momentum), the angular-momentum deficit (a secular-stability measure),
and per-planet osculating orbital elements.
"""

from __future__ import annotations

import math

from .builder import SimContext

# Orbital-element fields extracted at every output sample.
ELEMENT_FIELDS = ("a", "e", "inc", "Omega", "omega", "pomega", "M", "l")


def total_energy(ctx: SimContext) -> float:
    """Total energy including the active GR contribution where available.

    For ``gr_potential`` the modified potential is added to the Newtonian
    energy; for ``gr``/``gr_full`` REBOUNDx provides the full Hamiltonian.
    Falls back to the Newtonian energy if a helper is unavailable.
    """
    sim = ctx.sim
    e_newton = sim.energy()
    mode = ctx.gr_mode
    if mode == "none" or ctx.gr_force is None or ctx.rebx is None:
        return e_newton
    rebx = ctx.rebx
    try:
        if mode == "gr_potential":
            return e_newton + rebx.gr_potential_potential(ctx.gr_force)
        if mode == "gr":
            return rebx.gr_hamiltonian(ctx.gr_force)
        if mode == "gr_full":
            return rebx.gr_full_hamiltonian(ctx.gr_force)
    except Exception:
        pass
    return e_newton


def angular_momentum_magnitude(ctx: SimContext) -> float:
    lx, ly, lz = list(ctx.sim.angular_momentum())
    return math.sqrt(lx * lx + ly * ly + lz * lz)


def angular_momentum_deficit(ctx: SimContext) -> float:
    """Total AMD = sum_i m_i sqrt(a_i) [1 - sqrt(1-e_i^2) cos i_i]  (G=M_sun=1).

    AMD is conserved by the secular (averaged) dynamics; its growth is a clean
    flag that a planet's eccentricity/inclination is being pumped.
    """
    sim, primary = ctx.sim, ctx.sim.particles[0]
    amd = 0.0
    for i in range(1, ctx.n_bodies):
        p = sim.particles[i]
        o = p.orbit(primary=primary)
        if o.a <= 0 or o.e >= 1.0:
            continue
        circ = math.sqrt(o.a)
        amd += ctx.masses[ctx.body_names[i]] * circ * (
            1.0 - math.sqrt(1.0 - o.e * o.e) * math.cos(o.inc))
    return amd


def orbital_elements(ctx: SimContext) -> dict[str, dict[str, float]]:
    """Osculating elements of each planet relative to the Sun."""
    sim, primary = ctx.sim, ctx.sim.particles[0]
    out: dict[str, dict[str, float]] = {}
    for i in range(1, ctx.n_bodies):
        name = ctx.body_names[i]
        o = sim.particles[i].orbit(primary=primary)
        out[name] = {f: getattr(o, f) for f in ELEMENT_FIELDS}
    return out


def positions(ctx: SimContext) -> dict[str, tuple[float, float, float]]:
    """Heliocentric Cartesian positions (AU) of each planet."""
    sim, primary = ctx.sim, ctx.sim.particles[0]
    out = {}
    for i in range(1, ctx.n_bodies):
        p = sim.particles[i]
        out[ctx.body_names[i]] = (p.x - primary.x, p.y - primary.y, p.z - primary.z)
    return out


def lyapunov_time_years(ctx: SimContext) -> float:
    """Estimate the Lyapunov time (years) from the current Lyapunov number.

    Returns ``inf`` for regular motion / unconverged estimates.
    """
    if not ctx.config.chaos.megno.enabled:
        return float("inf")
    try:
        lcn = ctx.sim.lyapunov()
    except Exception:
        return float("inf")
    if lcn is None or lcn <= 0 or not math.isfinite(lcn):
        return float("inf")
    # lcn is in 1/(yr/2pi); T_lyap[yr] = (1/lcn)/(2*pi).
    return (1.0 / lcn) / (2.0 * math.pi)
