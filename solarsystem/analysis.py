"""Post-processing: summary statistics, conservation checks, chaos indicators
and secular (Fourier) spectra computed from a :class:`RunResult`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from . import constants
from .config import Config
from .runner import RunResult

ARCSEC_PER_CYCLE = 360.0 * 3600.0  # 1 cycle/yr -> arcsec/yr
RAD2DEG = 180.0 / math.pi


# --------------------------------------------------------------------------- #
# Per-planet element statistics
# --------------------------------------------------------------------------- #
def per_planet_stats(result: RunResult) -> list[dict[str, Any]]:
    rows = []
    for name in result.planet_names:
        el = result.elements[name]
        row: dict[str, Any] = {"planet": name, "mass_Msun": result.meta["masses"][name]}
        for fld, conv, unit in (("a", 1.0, "AU"),
                                ("e", 1.0, ""),
                                ("inc", RAD2DEG, "deg")):
            arr = el[fld] * conv
            row[f"{fld}_init{('_'+unit) if unit else ''}"] = float(arr[0])
            row[f"{fld}_mean"] = float(np.mean(arr))
            row[f"{fld}_std"] = float(np.std(arr))
            row[f"{fld}_min"] = float(np.min(arr))
            row[f"{fld}_max"] = float(np.max(arr))
            row[f"{fld}_range"] = float(np.max(arr) - np.min(arr))
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# Conservation & AMD
# --------------------------------------------------------------------------- #
def conservation_summary(result: RunResult, cfg: Config) -> dict[str, Any]:
    max_dE = float(np.max(result.dE_rel))
    max_dL = float(np.max(result.dL_rel))
    amd0 = float(result.amd[0]) if np.isfinite(result.amd[0]) else math.nan
    amd_max = float(np.nanmax(result.amd)) if np.any(np.isfinite(result.amd)) else math.nan
    amd_final = float(result.amd[-1]) if np.isfinite(result.amd[-1]) else math.nan
    return {
        "max_dE_rel": max_dE,
        "max_dL_rel": max_dL,
        "energy_tol": cfg.diagnostics.energy_tol,
        "angmom_tol": cfg.diagnostics.angmom_tol,
        "energy_ok": bool(max_dE <= cfg.diagnostics.energy_tol),
        "angmom_ok": bool(max_dL <= cfg.diagnostics.angmom_tol),
        "amd_initial": amd0,
        "amd_final": amd_final,
        "amd_max": amd_max,
        "amd_growth": float(amd_max / amd0) if amd0 and np.isfinite(amd0) else math.nan,
    }


# --------------------------------------------------------------------------- #
# Chaos indicators
# --------------------------------------------------------------------------- #
def chaos_summary(result: RunResult, cfg: Config) -> dict[str, Any]:
    if result.megno is None:
        return {"megno_enabled": False}
    megno = result.megno
    times = result.times_years
    megno_final = float(megno[-1])

    # Independent Lyapunov-time estimate from the late-time MEGNO slope:
    # for a chaotic orbit  <Y>(t) -> (lambda_max/2) t, so lambda = 2 dY/dt.
    t_lyap_slope = math.inf
    half = times.size // 2
    if times.size >= 8 and times[-1] > times[half]:
        slope, _ = np.polyfit(times[half:], megno[half:], 1)
        if slope > 0:
            lam = 2.0 * slope                # 1/yr
            if lam > 0:
                t_lyap_slope = 1.0 / lam

    classification = "regular" if megno_final < 4.0 else "chaotic"
    # Honest caveat: chaos only shows once t >> Lyapunov time (~5 Myr inner SS).
    converged = times[-1] >= 2.0 * constants.REFERENCE_LYAPUNOV_TIME_YEARS
    return {
        "megno_enabled": True,
        "megno_final": megno_final,
        "lyap_time_years_running": result.meta.get("lyap_time_years_final"),
        "lyap_time_years_slope": (None if not math.isfinite(t_lyap_slope)
                                  else float(t_lyap_slope)),
        "classification": classification,
        "diagnosis_converged": bool(converged),
        "reference_lyap_time_years": constants.REFERENCE_LYAPUNOV_TIME_YEARS,
    }


# --------------------------------------------------------------------------- #
# Secular spectra (FFT of the complex eccentricity / inclination vectors)
# --------------------------------------------------------------------------- #
def _spectrum(times_years: np.ndarray, complex_signal: np.ndarray):
    n = complex_signal.size
    if n < 16:
        return None
    # Resample onto a strictly uniform grid to kill integer-step jitter.
    t_uni = np.linspace(times_years[0], times_years[-1], n)
    re = np.interp(t_uni, times_years, complex_signal.real)
    im = np.interp(t_uni, times_years, complex_signal.imag)
    z = (re + 1j * im)
    z = z - z.mean()
    z = z * np.hanning(n)
    dt = (t_uni[-1] - t_uni[0]) / (n - 1)
    spec = np.fft.fftshift(np.fft.fft(z))
    freq = np.fft.fftshift(np.fft.fftfreq(n, d=dt)) * ARCSEC_PER_CYCLE  # arcsec/yr
    power = np.abs(spec) ** 2
    return freq, power


def secular_spectra(result: RunResult, cfg: Config) -> dict[str, Any]:
    if cfg.time.spacing != "linear":
        return {"available": False, "reason": "spacing != linear"}
    times = result.times_years
    span = float(times[-1] - times[0]) if times.size > 1 else 0.0
    if span <= 0 or times.size < 16:
        return {"available": False, "reason": "too few samples / zero span"}
    out: dict[str, Any] = {"available": True, "g": {}, "s": {},
                           "resolution_arcsec_yr": ARCSEC_PER_CYCLE / span}
    for name in result.planet_names:
        el = result.elements[name]
        # g-spectrum: complex eccentricity z = e * exp(i*pomega)  -> peaks near g_i
        zg = el["e"] * np.exp(1j * el["pomega"])
        sg = _spectrum(times, zg)
        if sg is not None:
            freq, power = sg
            out["g"][name] = {"freq": freq, "power": power,
                              "peak": _dominant_peak(freq, power, 0.1, 60.0)}
        # s-spectrum: complex inclination zeta = sin(i/2) * exp(i*Omega) -> s_i
        zs = np.sin(el["inc"] / 2.0) * np.exp(1j * el["Omega"])
        ss = _spectrum(times, zs)
        if ss is not None:
            freq, power = ss
            out["s"][name] = {"freq": freq, "power": power,
                              "peak": _dominant_peak(freq, power, -60.0, -0.1)}
    return out


def _dominant_peak(freq, power, lo, hi):
    mask = (freq >= lo) & (freq <= hi)
    if not np.any(mask):
        return None
    sub_f, sub_p = freq[mask], power[mask]
    return float(sub_f[int(np.argmax(sub_p))])


# --------------------------------------------------------------------------- #
# Top-level summary + persistence
# --------------------------------------------------------------------------- #
def summarize(result: RunResult, cfg: Config) -> dict[str, Any]:
    return {
        "meta": result.meta,
        "conservation": conservation_summary(result, cfg),
        "chaos": chaos_summary(result, cfg),
        "planets": per_planet_stats(result),
    }


def save_statistics(summary: dict, outdir: str | Path) -> None:
    outdir = Path(outdir)
    (outdir / "summary.json").write_text(_json(summary))
    # Flat per-planet CSV.
    rows = summary["planets"]
    if rows:
        cols = list(rows[0].keys())
        lines = [",".join(cols)]
        for r in rows:
            lines.append(",".join(_csv_cell(r[c]) for c in cols))
        (outdir / "statistics.csv").write_text("\n".join(lines) + "\n")


def _csv_cell(v) -> str:
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def _json(obj) -> str:
    def default(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    return json.dumps(obj, indent=2, default=default)
