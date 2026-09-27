"""Measure the inner-Solar-System maximum Lyapunov time (Laskar's ~5 Myr result).

Method: a parallel **MEGNO ensemble**.

We tried the two-trajectory (Benettin) shadow method first and found it useless
here: the position-space separation of two nearby orbits is dominated by regular
**Keplerian phase-shearing** (nearby orbits with slightly different periods drift
apart on a ~10^3 yr orbital timescale and saturate at ~orbit scale), which
completely masks the slow ~5 Myr secular chaos.  MEGNO is built precisely to
filter that: its time-averaged ⟨Y⟩ → 2 for regular/quasi-periodic motion (even
when |δ| grows linearly) and grows linearly, slope λ_max/2, only for genuine
exponential divergence.  So MEGNO is the correct estimator for this problem.

We integrate an ensemble of ``n_traj`` full 8-planet+GR systems, each with a
tiny independent IC kick and its own variational seed, in parallel threads
(REBOUND releases the GIL).  For each we read REBOUND's running Lyapunov number
and also fit the late-time slope of ⟨Y⟩(t).  The Lyapunov time is T = 1/λ_max.

Note (documented framework limitation): REBOUNDx GR is not propagated through the
variational equations, so MEGNO measures the Newtonian tangent-flow growth around
the GR-corrected reference orbit.  The inner-system Lyapunov *time* is known to be
robust to GR (~5 Myr with or without it; GR mainly changes the long-term
instability probability), so this does not bias the timescale being verified.

Usage:
    .venv/bin/python experiments/lyapunov_inner.py --t-max-myr 25 --name lyap
"""

from __future__ import annotations

import argparse
import csv
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from solarsystem import constants, load_config
from solarsystem.builder import build_simulation

TWO_PI = 2.0 * math.pi


def lyap_time_years(sim) -> float:
    """Running Lyapunov time (yr) from REBOUND's Lyapunov number (yr2pi units)."""
    try:
        lcn = sim.lyapunov()
    except Exception:
        return math.inf
    if not lcn or lcn <= 0 or not math.isfinite(lcn):
        return math.inf
    return (1.0 / lcn) / TWO_PI


def main() -> int:
    ap = argparse.ArgumentParser(description="Inner-SS Lyapunov time via MEGNO ensemble")
    ap.add_argument("--t-max-myr", type=float, default=25.0)
    ap.add_argument("--chunk-myr", type=float, default=0.25)
    ap.add_argument("--n-traj", type=int, default=8)
    ap.add_argument("--kick", type=float, default=1e-9, help="IC perturbation [AU]")
    ap.add_argument("--dt-days", type=float, default=4.0)
    ap.add_argument("--gr", default="gr_potential")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--name", default="lyapunov_inner")
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    contexts = []
    for i in range(a.n_traj):
        cfg = load_config("config/default.yaml", overrides=[
            "chaos.megno.enabled=true",
            f"chaos.megno.seed={1000 + i}",
            f"integrator.dt_days={a.dt_days}",
            f"physics.gr.mode={a.gr}",
            "archive.enabled=false",
        ])
        ctx = build_simulation(cfg)
        # tiny independent IC kick on Mercury so trajectories are independent
        p = ctx.sim.particles[1]
        d = rng.standard_normal(3); d /= np.linalg.norm(d)
        p.x += a.kick * d[0]; p.y += a.kick * d[1]; p.z += a.kick * d[2]
        ctx.sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        contexts.append(ctx)
    sims = [c.sim for c in contexts]

    chunk = constants.years_to_code_time(a.chunk_myr * 1e6)
    n_chunks = int(round(a.t_max_myr / a.chunk_myr))

    outdir = Path("outputs") / a.name
    outdir.mkdir(parents=True, exist_ok=True)
    fcsv = open(outdir / "lyapunov.csv", "w", newline="")
    writer = csv.writer(fcsv)
    writer.writerow(["t_Myr", "megno_mean", "megno_std",
                     "Tlyap_mean_Myr", "Tlyap_std_Myr", "Tlyap_slope_Myr"]
                    + [f"megno_{i}" for i in range(a.n_traj)])
    fcsv.flush()

    print(f"[lyapunov-MEGNO] {contexts[0].n_bodies - 1} planets + GR({a.gr}), "
          f"dt={a.dt_days} d, {a.n_traj} trajectories, target {a.t_max_myr:g} Myr "
          f"(chunk {a.chunk_myr:g} Myr)", flush=True)

    t_hist: list[float] = []
    megno_hist: list[list[float]] = []
    pool = ThreadPoolExecutor(max_workers=len(sims))

    def advance(sm):
        sm.integrate(sm.t + chunk, exact_finish_time=0)
        sm.synchronize()

    t0 = time.time()
    for k in range(1, n_chunks + 1):
        list(pool.map(advance, sims))
        t_yr = constants.code_time_to_years(sims[0].t)
        megnos = [s.megno() for s in sims]
        Tls = [lyap_time_years(s) for s in sims]
        finite = [T for T in Tls if math.isfinite(T)]

        t_hist.append(t_yr)
        megno_hist.append(megnos)

        # slope-based estimate: <Y> ~ (lambda/2) t  => lambda = 2 d<Y>/dt
        T_slope = math.inf
        if len(t_hist) >= 8:
            tarr = np.array(t_hist); yarr = np.array(megno_hist).mean(axis=1)
            half = len(tarr) // 2
            slope = np.polyfit(tarr[half:], yarr[half:], 1)[0]
            if slope > 0:
                T_slope = 1.0 / (2.0 * slope)

        megno_mean = float(np.mean(megnos)); megno_std = float(np.std(megnos))
        T_mean = float(np.mean(finite)) if finite else math.inf
        T_std = float(np.std(finite)) if len(finite) > 1 else 0.0
        writer.writerow([t_yr / 1e6, megno_mean, megno_std,
                         T_mean / 1e6 if math.isfinite(T_mean) else "",
                         T_std / 1e6,
                         T_slope / 1e6 if math.isfinite(T_slope) else ""]
                        + megnos)
        fcsv.flush()

        rate = (t_yr / 1e6) / max((time.time() - t0) / 3600.0, 1e-9)
        print(f"  t={t_yr/1e6:6.2f} Myr | <Y>={megno_mean:5.2f}±{megno_std:4.2f} | "
              f"T_lyap(run)={T_mean/1e6 if math.isfinite(T_mean) else float('nan'):6.2f} Myr | "
              f"T_lyap(slope)={T_slope/1e6 if math.isfinite(T_slope) else float('nan'):6.2f} Myr | "
              f"{rate:4.1f} Myr/hr", flush=True)

    fcsv.close()
    _final(outdir, np.array(t_hist), np.array(megno_hist), T_slope)
    return 0


def _final(outdir, t_hist, megno_hist, T_slope) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = t_hist / 1e6
    ymean = megno_hist.mean(axis=1)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for j in range(megno_hist.shape[1]):
        ax[0].plot(t, megno_hist[:, j], lw=0.6, alpha=0.5)
    ax[0].plot(t, ymean, color="#ffb454", lw=2.0, label="ensemble mean")
    ax[0].axhline(2.0, color="#5ad1c2", ls="--", lw=1, label="regular ⟨Y⟩→2")
    ax[0].set_xlabel("time [Myr]"); ax[0].set_ylabel("MEGNO ⟨Y⟩")
    ax[0].set_title("MEGNO growth (chaos ⇒ linear rise)"); ax[0].legend(); ax[0].grid(alpha=0.3)

    # running Lyapunov-time convergence (slope-based, expanding window)
    Truns = []
    for n in range(4, len(t) + 1):
        half = n // 2
        sl = np.polyfit(t_hist[half:n], ymean[half:n], 1)[0]
        Truns.append(1.0 / (2.0 * sl) / 1e6 if sl > 0 else np.nan)
    ax[1].plot(t[3:], Truns, color="#5ad1c2", lw=1.8)
    ax[1].axhline(5.0, color="#ff6b6b", ls="--", lw=1.2, label="Laskar ~5 Myr")
    ax[1].set_ylim(0, 15)
    ax[1].set_xlabel("time [Myr]"); ax[1].set_ylabel("Lyapunov time [Myr]")
    ax[1].set_title("Lyapunov-time estimate (MEGNO slope)"); ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(outdir / "lyapunov_convergence.png", dpi=140)

    T = T_slope / 1e6 if math.isfinite(T_slope) else float("nan")
    verdict = "MATCH" if 3.0 <= T <= 8.0 else "outside 3-8 Myr"
    msg = (f"Final inner-SS Lyapunov time (MEGNO slope): {T:.2f} Myr at {t[-1]:.1f} Myr "
           f"integration  ->  {verdict} vs Laskar ~5 Myr\n"
           f"Final ensemble ⟨Y⟩ = {ymean[-1]:.2f}\n")
    (outdir / "verdict.txt").write_text(msg)
    print("\n" + msg)


if __name__ == "__main__":
    raise SystemExit(main())
