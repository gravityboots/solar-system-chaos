# CEAM Solar-System Chaos & Stability Framework

A configuration-driven **REBOUND + REBOUNDx** framework for symplectic N-body
integration of the Solar System, built for the CEAM research project
*"Numerical Investigation of Chaotic Dynamics and Stability in the Solar
System."* It produces the standard statistical outputs (conserved quantities,
orbital-element statistics, MEGNO/Lyapunov chaos indicators, secular spectra)
and a full set of publication-style figures and a Markdown report for every run.

The design philosophy follows the project brief: *this is a science-with-existing-tools
project, so the intellectual content lives in the experimental design.* Every
scientifically meaningful knob — which bodies, their masses and elements, GR,
the timestep, the integrator — is **config, not code**, so the planned chaos
experiments (vary Jupiter's mass, remove a planet, toggle GR, …) are one-liners.

> Status: the simulation + diagnostics + visualization framework is complete and
> verified. The chaos *experiments* themselves (parameter sweeps, ensembles) are
> the next phase and slot directly onto this framework — see
> [Experiment recipes](#experiment-recipes).

---

## Installation

Requires Python ≥ 3.10 and a C compiler (REBOUND/REBOUNDx build from source).

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip setuptools wheel
.venv/bin/python -m pip install -e .          # installs deps + the `ceam-sim` CLI
```

Verify:

```bash
.venv/bin/python tests/test_smoke.py
```

You should see machine-precision energy conservation (|ΔE/E| ~ 1e-15 for a clean
2-body run) and a passing quick pipeline.

---

## Quick start

```bash
# Default 1-Myr, 8-planet run with GR + MEGNO -> outputs/solar_default/
.venv/bin/ceam-sim

# Fast pipeline check (~20 kyr, low-res)
.venv/bin/ceam-sim --quick --name smoke

# A longer reference run
.venv/bin/ceam-sim --t-end 5e6 --name inner_5Myr
```

Or from Python:

```python
from solarsystem import load_config, run_pipeline

cfg = load_config("config/default.yaml",
                  overrides=["time.t_end_years=2e6", "physics.gr.mode=none"])
out = run_pipeline(cfg)
print(out.summary["conservation"]["max_dE_rel"])
```

Every run writes to `outputs/<run.name>/` (see [Outputs](#outputs)).

---

## Configuration

All parameters live in [`config/default.yaml`](config/default.yaml) and are fully
commented. Override any of them three ways (low → high precedence): the YAML file,
then `--set a.b.c=value` flags, then the convenience flags below. Units in the
config are human-friendly (days, Julian years, degrees); the framework converts
to REBOUND's `('AU', 'yr2pi', 'Msun')` system (G = 1) internally.

| Section | Key | Meaning |
|---|---|---|
| `integrator` | `name` | `whfast` (production), `ias15` (validation), `trace`, `saba` |
| | `dt_days` | timestep in days (4 d production, 1 d for convergence checks) |
| | `whfast.corrector` | symplectic corrector order (11 = WHCKL-class) |
| | `whfast.kernel` | `lazy` (WHCKL) — auto-falls back to `default` when MEGNO is on |
| `time` | `t_end_years` | integration span |
| | `n_outputs` | number of recorded snapshots (time-series resolution) |
| | `spacing` | `linear` (required for the secular FFT) or `log` |
| `bodies` | `source` | `builtin` (offline J2000 table), `horizons`, or `file` |
| | `include` | which planets, in order |
| | `overrides.<Planet>` | per-body `mass_scale`, `mass`, `a`, `e`, `inc`, `Omega`, `pomega`, `l`, `include` |
| `physics.gr` | `mode` | `none`, `gr_potential` (recommended), `gr`, `gr_full` |
| `physics.j2` | `enabled` | Sun oblateness term (tiny; for sensitivity checks) |
| `chaos.megno` | `enabled` | MEGNO + Lyapunov via variational particles |
| `freeze` | `bodies` | hold listed bodies on their initial orbit (experimental) |
| `diagnostics` | `energy_tol`, `angmom_tol` | pass/fail thresholds for conservation |
| `plots` | `style`, `dpi`, `animate` | figure appearance and an optional orbit GIF |

Convenience flags: `--t-end`, `--dt`, `--gr`, `--name`, `--no-megno`,
`--no-plots`, `--animate`, `--quick`.

### Refreshing initial conditions from JPL Horizons

The builtin J2000 table keeps everything offline and reproducible. To pull a
specific epoch from Horizons and cache it:

```bash
.venv/bin/ceam-sim horizons --epoch "2025-01-01 00:00" --out cache/horizons_2025.yaml
.venv/bin/ceam-sim --set bodies.source=file --set bodies.file=cache/horizons_2025.yaml
```

---

## Outputs

```
outputs/<run.name>/
├── config.resolved.yaml   # exact resolved config (reproducibility)
├── archive.bin            # REBOUND SimulationArchive (bit-reproducible restart/analysis)
├── timeseries.npz         # all time series (energy, AM, AMD, MEGNO, elements, positions)
├── summary.json           # conservation, chaos, per-planet statistics
├── statistics.csv         # per-planet element statistics table
├── report.md              # human-readable report with embedded figures
└── figures/
    ├── orbits.png             # heliocentric orbit traces (inner / outer)
    ├── conservation.png       # |ΔE/E| and |ΔL/L| vs time
    ├── eccentricity.png       # e(t) per planet
    ├── inclination.png        # i(t) per planet
    ├── semimajor.png          # (a−a₀)/a₀ in ppm (bounded for symplectic)
    ├── body_panel.png         # a, e, i, ϖ for Mercury (the headline body)
    ├── megno.png              # MEGNO ⟨Y⟩ vs time with regular/chaotic references
    ├── amd.png                # angular-momentum deficit (secular stability)
    ├── secular_spectrum.png   # FFT of e·e^{iϖ} with Laskar g-frequency markers
    ├── hk_phase.png           # secular phase portrait (h=e cosϖ, k=e sinϖ)
    └── dashboard.png          # one-page overview
```

Reload a run for further analysis:

```python
import numpy as np
d = np.load("outputs/solar_default/timeseries.npz")
t = d["times_years"]; e_merc = d["elem__Mercury__e"]
import rebound
sa = rebound.Simulationarchive("outputs/solar_default/archive.bin")
```

---

## Experiment recipes

These are the scoped experiments from the project brief (Deliverables §3.4),
expressed directly as configuration. They are the natural next phase.

```bash
# --- GR on/off (the cheapest, most decisive chaos-sensitivity experiment) ---
.venv/bin/ceam-sim --gr gr_potential --t-end 5e6 --name gr_on
.venv/bin/ceam-sim --gr none        --t-end 5e6 --name gr_off

# --- Experiment B: Jupiter-mass MEGNO scan (one grid point shown) ---
.venv/bin/ceam-sim --set bodies.overrides.Jupiter.mass_scale=1.5 \
                   --t-end 5e6 --name jupiter_x1.5

# --- Experiment C: (e_Jup, e_Sat) plane (one grid point shown) ---
.venv/bin/ceam-sim --set bodies.overrides.Jupiter.e=0.06 \
                   --set bodies.overrides.Saturn.e=0.08 --name ejup_esat

# --- Experiment D: planet-removed cross-check ---
.venv/bin/ceam-sim --set bodies.overrides.Saturn.include=false --name no_saturn
```

A full sweep is then a loop over `--set` values calling `run_pipeline` per point,
each in its own process (`multiprocessing.Pool`) — the embarrassingly-parallel
model the brief recommends. Each worker must build its own `Simulation`.

---

## Scientific notes & honest caveats

- **Integrator**: WHFast with an 11th-order symplectic corrector (WHCKL-class),
  Jacobi coordinates, `safe_mode=0` (the framework `synchronize()`s before
  reading state). Energy is bounded (not secular) and angular momentum is
  conserved to ~machine precision — the standard validation the brief asks for.
- **GR**: `gr_potential` on the Sun (Nobili & Roxburgh 1986) preserves WHFast's
  symplecticity and gives Mercury's 43″/century precession. This is decisive for
  the science: without GR, Mercury's destabilisation probability over 5 Gyr rises
  from ~1% to >60% (Laskar 2008).
- **MEGNO + GR caveat**: REBOUNDx forces are **not** propagated through REBOUND's
  variational equations, so MEGNO measures the *Newtonian* variational growth
  around the GR-corrected reference trajectory. This is fine for qualitative chaos
  detection; for a rigorous GR Lyapunov exponent use a two-shadow-trajectory
  estimate (a planned addition) or compute MEGNO with GR off. The lazy/WHCKL
  kernel is also incompatible with variational particles, so MEGNO runs use the
  standard kernel automatically.
- **Run length for chaos**: the inner-system Lyapunov time is ~5 Myr, so MEGNO
  only diverges from ⟨Y⟩≈2 after *tens* of Myr. The 1-Myr default exercises the
  full pipeline and gives clean secular dynamics, but is *too short to diagnose
  chaos* — the report flags this. For chaos, run ≳ 20–50 Myr.
- **Secular spectrum resolution** is 1/T, so resolving g₁ (5.59″/yr) from g₅
  (4.26″/yr) needs multi-Myr runs; the figure annotates its own resolution.

### Key references (in `context/`)

Laskar 2008 (Icarus 196) · Batygin & Laughlin 2008 · Ito & Tanikawa 2002 ·
Brown & Rein 2020, 2023 · Mogavero & Laskar 2022 · Mogavero, Hoang & Laskar 2023.
See `context/Deliverables.md` for the full project brief.

---

## Project structure

```
solarsystem/
├── constants.py     # units, J2000 body table, reference secular frequencies
├── config.py        # dataclass config model + YAML/dotted-override loader
├── builder.py       # build rebound.Simulation (+ REBOUNDx GR/J2, MEGNO)
├── horizons.py      # optional JPL Horizons initial conditions
├── diagnostics.py   # energy/AM/AMD + orbital-element extraction
├── runner.py        # chunked symplectic integration + SimulationArchive
├── analysis.py      # statistics, conservation checks, secular spectra
├── visualize.py     # all figures
├── report.py        # Markdown report generator
├── pipeline.py      # build → run → analyse → visualise → report
└── cli.py           # `ceam-sim` command-line interface
config/default.yaml  # the single place to tweak parameters
tests/test_smoke.py  # fast sanity checks
```
