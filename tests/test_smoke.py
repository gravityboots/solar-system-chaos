"""Fast sanity checks for the framework. Run with pytest or directly:

    .venv/bin/python tests/test_smoke.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from solarsystem import build_simulation, load_config, run_pipeline
from solarsystem.config import Config
from solarsystem import diagnostics


def _base_config(**over) -> Config:
    ov = [f"{k}={v}" for k, v in over.items()]
    return load_config(None, ov)


def test_build_default_has_nine_bodies():
    cfg = _base_config(**{"chaos.megno.enabled": "false"})
    ctx = build_simulation(cfg)
    assert ctx.n_bodies == 9            # Sun + 8 planets
    assert ctx.planet_names[0] == "Mercury"
    assert ctx.sim.particles[0].m == 1.0
    print("OK build: 9 bodies")


def test_two_body_energy_conservation():
    # Sun + Jupiter only, no GR, no MEGNO: WHFast must conserve energy tightly.
    cfg = _base_config(**{
        "bodies.include": "[Jupiter]",
        "physics.gr.mode": "none",
        "chaos.megno.enabled": "false",
        "time.t_end_years": "1e4",
        "time.n_outputs": "200",
    })
    ctx = build_simulation(cfg)
    e0 = diagnostics.total_energy(ctx)
    ctx.sim.integrate(ctx.sim.t + 1000.0)
    ctx.sim.synchronize()
    e1 = diagnostics.total_energy(ctx)
    rel = abs(e1 - e0) / abs(e0)
    assert rel < 1e-9, f"energy drift too large: {rel:.2e}"
    print(f"OK 2-body energy conservation: |dE/E| = {rel:.2e}")


def test_quick_pipeline_runs(tmp_path: Path | None = None):
    root = tmp_path or Path(tempfile.mkdtemp())
    cfg = _base_config(**{
        "run.name": "pytest_quick",
        "run.output_root": str(root),
        "time.t_end_years": "2e4",
        "time.n_outputs": "300",
        "plots.enabled": "false",
    })
    out = run_pipeline(cfg, progress=False)
    cons = out.summary["conservation"]
    assert (out.outdir / "timeseries.npz").exists()
    assert (out.outdir / "summary.json").exists()
    assert (out.outdir / "archive.bin").exists()
    assert cons["max_dE_rel"] < 1e-6, cons["max_dE_rel"]
    assert cons["max_dL_rel"] < 1e-8, cons["max_dL_rel"]
    # MEGNO of a 20-kyr inner-system run should sit near the regular value 2.
    assert out.result.megno is not None
    assert 1.0 < out.result.megno[-1] < 4.0
    # Eccentricities stay physical.
    for name in out.result.planet_names:
        e = out.result.elements[name]["e"]
        assert np.all((e >= 0) & (e < 1))
    print(f"OK quick pipeline: max|dE/E|={cons['max_dE_rel']:.2e}, "
          f"MEGNO={out.result.megno[-1]:.2f}")


if __name__ == "__main__":
    test_build_default_has_nine_bodies()
    test_two_body_energy_conservation()
    test_quick_pipeline_runs()
    print("\nAll smoke tests passed.")
