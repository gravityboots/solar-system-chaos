"""Optional initial conditions from NASA JPL Horizons.

The builtin J2000 table (``constants.DEFAULT_BODIES``) keeps the framework fully
offline and reproducible.  This module refreshes the *elements* for a chosen
epoch from Horizons and caches them to a YAML body-table that ``bodies.source:
file`` can consume.  Masses are taken from the builtin table (they are physical
constants, not epoch-dependent), so only the time-varying elements are fetched.
"""

from __future__ import annotations

import math
from pathlib import Path

import yaml

from . import constants

DEG = 180.0 / math.pi


def fetch_body_table(names: list[str], epoch: str = "2000-01-01 12:00",
                     earth_moon: str = "barycenter") -> dict[str, dict]:
    """Query Horizons for each body's state at ``epoch`` and return a table.

    Requires network access. Raises on failure (callers should fall back to the
    builtin table if they want offline robustness).
    """
    import rebound

    sim = rebound.Simulation()
    sim.units = constants.UNITS
    sim.add("Sun", date=epoch, hash="Sun")

    table: dict[str, dict] = {}
    for name in names:
        query = name
        if name == "Earth" and earth_moon == "barycenter":
            query = "3"  # Earth-Moon barycentre (Horizons body 3)
        try:
            sim.add(query, date=epoch, hash=name)
        except Exception:
            if name == "Earth":
                sim.add("Geocenter", date=epoch, hash=name)
            else:
                raise

    sun = sim.particles[0]
    for name in names:
        o = sim.particles[name].orbit(primary=sun)
        base = constants.DEFAULT_BODIES.get(name, {})
        table[name] = dict(
            mass=base.get("mass", sim.particles[name].m),
            a=o.a, e=o.e, inc=o.inc * DEG,
            Omega=o.Omega * DEG, pomega=o.pomega * DEG, l=o.l * DEG,
            group=base.get("group", "outer"), color=base.get("color", "#cccccc"),
            epoch=epoch,
        )
    return table


def save_body_table(table: dict[str, dict], path: str | Path) -> None:
    Path(path).write_text(yaml.safe_dump(table, sort_keys=False))


def load_body_table(path: str | Path) -> dict[str, dict]:
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"Body table {path} is not a mapping")
    return data
