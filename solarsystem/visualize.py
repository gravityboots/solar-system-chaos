"""Figure generation.

Every figure is a standalone function ``fig_*(result, cfg, ...) -> Path`` plus a
:func:`make_all` driver.  Plots use a consistent dark "space" theme (or light)
and the per-planet colours defined in :mod:`solarsystem.constants`.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import constants  # noqa: E402
from .config import Config  # noqa: E402
from .runner import RunResult  # noqa: E402

RAD2DEG = 180.0 / math.pi


# --------------------------------------------------------------------------- #
# Styling helpers
# --------------------------------------------------------------------------- #
def _theme(cfg: Config) -> dict:
    if cfg.plots.style == "light":
        return dict(bg="#ffffff", panel="#f4f4f8", fg="#1a1a1a",
                    grid="#cccccc", muted="#666666")
    return dict(bg="#0b0e17", panel="#141a28", fg="#e8ecf5",
                grid="#2a3550", muted="#8a93a8")


def _apply(theme: dict) -> None:
    plt.rcParams.update({
        "figure.facecolor": theme["bg"], "axes.facecolor": theme["panel"],
        "savefig.facecolor": theme["bg"], "text.color": theme["fg"],
        "axes.labelcolor": theme["fg"], "axes.edgecolor": theme["grid"],
        "xtick.color": theme["muted"], "ytick.color": theme["muted"],
        "grid.color": theme["grid"], "axes.titlecolor": theme["fg"],
        "axes.grid": True, "grid.alpha": 0.35, "font.size": 10,
        "axes.titlesize": 11, "figure.titlesize": 13, "legend.framealpha": 0.0,
    })


def _color(name: str) -> str:
    return constants.DEFAULT_BODIES.get(name, {}).get("color", "#cccccc")


def _time_axis(times_years: np.ndarray):
    span = times_years[-1] - times_years[0]
    if span >= 2e6:
        return times_years / 1e6, "Time [Myr]"
    if span >= 2e3:
        return times_years / 1e3, "Time [kyr]"
    return times_years, "Time [yr]"


def _decimate(n: int, max_points: int) -> slice:
    return slice(None, None, max(1, n // max_points))


def _save(fig, path: Path, cfg: Config) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=cfg.plots.dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _present(result: RunResult, group: list[str]) -> list[str]:
    return [n for n in result.planet_names if n in group]


def _ellipse_xy(a, e, inc, Omega, omega, n=512):
    """Ecliptic-projected (x, y) of an osculating ellipse from its elements.

    Reconstructed analytically by sweeping true anomaly, so it is a clean curve
    regardless of the (coarse) output cadence — connecting stored position
    samples directly would alias badly over long runs.
    """
    f = np.linspace(0.0, 2.0 * np.pi, n)
    r = a * (1.0 - e * e) / (1.0 + e * np.cos(f))
    u = omega + f
    cosO, sinO, cosi = np.cos(Omega), np.sin(Omega), np.cos(inc)
    x = r * (cosO * np.cos(u) - sinO * np.sin(u) * cosi)
    y = r * (sinO * np.cos(u) + cosO * np.sin(u) * cosi)
    return x, y


def _draw_orbits(ax, result, group, theme, label=True):
    ax.scatter([0], [0], c=constants.SUN_COLOR, s=110, marker="*", zorder=6,
               edgecolors="none", label="Sun" if label else None)
    for name in group:
        el = result.elements[name]
        x0, y0 = _ellipse_xy(el["a"][0], el["e"][0], el["inc"][0],
                             el["Omega"][0], el["omega"][0])
        ax.plot(x0, y0, lw=1.1, color=_color(name), alpha=0.95,
                label=name if label else None)
        if el["a"].size > 1:  # final-epoch orbit (dashed) shows precession/drift
            x1, y1 = _ellipse_xy(el["a"][-1], el["e"][-1], el["inc"][-1],
                                 el["Omega"][-1], el["omega"][-1])
            ax.plot(x1, y1, lw=0.9, color=_color(name), alpha=0.4, ls="--")
    ax.set_aspect("equal")


# --------------------------------------------------------------------------- #
# Individual figures
# --------------------------------------------------------------------------- #
def fig_orbits(result: RunResult, cfg: Config, path: Path) -> Path:
    theme = _theme(cfg); _apply(theme)
    inner = _present(result, constants.INNER_PLANETS)
    outer = _present(result, constants.OUTER_PLANETS)
    groups = [(g, t) for g, t in ((inner, "Inner planets"), (outer, "Outer planets")) if g]
    fig, axes = plt.subplots(1, len(groups), figsize=(6 * len(groups), 6), squeeze=False)
    for ax, (group, title) in zip(axes[0], groups):
        ax.set_title(title)
        _draw_orbits(ax, result, group, theme)
        ax.set_xlabel("x [AU]"); ax.set_ylabel("y [AU]")
        ax.legend(loc="upper right", fontsize=8, ncol=2)
    fig.suptitle("Osculating orbits (ecliptic projection) — "
                 "solid = initial epoch, dashed = final epoch")
    return _save(fig, path, cfg)


def fig_conservation(result: RunResult, cfg: Config, path: Path) -> Path:
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].semilogy(t, np.clip(result.dE_rel, 1e-18, None), color="#5ad1c2", lw=1.2)
    axes[0].axhline(cfg.diagnostics.energy_tol, color="#ff6b6b", ls="--", lw=1,
                    label=f"tol = {cfg.diagnostics.energy_tol:g}")
    axes[0].set_ylabel(r"$|\Delta E / E_0|$"); axes[0].legend(loc="best", fontsize=8)
    axes[0].set_title("Relative energy error")
    axes[1].semilogy(t, np.clip(result.dL_rel, 1e-18, None), color="#c2a3ff", lw=1.2)
    axes[1].axhline(cfg.diagnostics.angmom_tol, color="#ff6b6b", ls="--", lw=1,
                    label=f"tol = {cfg.diagnostics.angmom_tol:g}")
    axes[1].set_ylabel(r"$|\Delta L / L_0|$"); axes[1].set_xlabel(tlabel)
    axes[1].legend(loc="best", fontsize=8); axes[1].set_title("Relative angular-momentum error")
    fig.suptitle("Conservation diagnostics (symplectic integrator)")
    return _save(fig, path, cfg)


def _element_vs_time(result, cfg, path, field, conv, ylabel, title):
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    inner = _present(result, constants.INNER_PLANETS)
    outer = _present(result, constants.OUTER_PLANETS)
    groups = [g for g in (inner, outer) if g]
    fig, axes = plt.subplots(len(groups), 1, figsize=(11, 3.2 * len(groups)),
                             squeeze=False)
    for ax, group, gname in zip(axes[:, 0], groups,
                                ["Inner planets", "Outer planets"]):
        for name in group:
            ax.plot(t, result.elements[name][field] * conv, color=_color(name),
                    lw=1.0, label=name)
        ax.set_ylabel(ylabel); ax.legend(loc="upper right", fontsize=8, ncol=4)
        ax.set_title(gname)
    axes[-1, 0].set_xlabel(tlabel)
    fig.suptitle(title)
    return _save(fig, path, cfg)


def fig_eccentricity(result, cfg, path):
    return _element_vs_time(result, cfg, path, "e", 1.0, "eccentricity",
                            "Eccentricity evolution")


def fig_inclination(result, cfg, path):
    return _element_vs_time(result, cfg, path, "inc", RAD2DEG, "inclination [deg]",
                            "Inclination evolution")


def fig_semimajor(result: RunResult, cfg: Config, path: Path) -> Path:
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    fig, ax = plt.subplots(figsize=(11, 5))
    for name in result.planet_names:
        a = result.elements[name]["a"]
        ax.plot(t, (a - a[0]) / a[0] * 1e6, color=_color(name), lw=1.0, label=name)
    ax.set_xlabel(tlabel); ax.set_ylabel(r"$(a - a_0)/a_0$  [ppm]")
    ax.legend(loc="upper right", fontsize=8, ncol=4)
    ax.set_title("Semi-major-axis variation (symplectic: should stay bounded)")
    return _save(fig, path, cfg)


def fig_body_panel(result: RunResult, cfg: Config, path: Path,
                   name: str | None = None) -> Path:
    theme = _theme(cfg); _apply(theme)
    if name is None:
        name = "Mercury" if "Mercury" in result.planet_names else result.planet_names[0]
    el = result.elements[name]
    t, tlabel = _time_axis(result.times_years)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    specs = [("a", 1.0, "a [AU]"), ("e", 1.0, "eccentricity"),
             ("inc", RAD2DEG, "inclination [deg]"),
             ("pomega", RAD2DEG, r"$\varpi$ [deg]")]
    for ax, (fld, conv, lab) in zip(axes.flat, specs):
        y = el[fld] * conv
        if fld == "pomega":
            y = np.mod(y, 360.0)
            ax.scatter(t, y, s=2, color=_color(name))
        else:
            ax.plot(t, y, color=_color(name), lw=1.0)
        ax.set_xlabel(tlabel); ax.set_ylabel(lab)
    fig.suptitle(f"{name}: orbital-element evolution")
    return _save(fig, path, cfg)


def fig_megno(result: RunResult, cfg: Config, path: Path) -> Path:
    if result.megno is None:
        return path
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(t, result.megno, color="#ffb454", lw=1.3, label="MEGNO  ⟨Y⟩")
    ax.axhline(2.0, color="#5ad1c2", ls="--", lw=1, label="regular (⟨Y⟩→2)")
    ax.axhline(4.0, color="#ff6b6b", ls=":", lw=1, label="chaotic threshold (~4)")
    ax.set_xlabel(tlabel); ax.set_ylabel("MEGNO")
    final = result.megno[-1]
    ax.set_title(f"MEGNO chaos indicator  (final ⟨Y⟩ = {final:.2f})")
    ax.legend(loc="best", fontsize=8)
    return _save(fig, path, cfg)


def fig_amd(result: RunResult, cfg: Config, path: Path) -> Path:
    if not np.any(np.isfinite(result.amd)):
        return path
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    amd0 = result.amd[0]
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(t, result.amd / amd0, color="#7ee787", lw=1.2)
    ax.axhline(1.0, color=theme["muted"], ls="--", lw=0.8)
    ax.set_xlabel(tlabel); ax.set_ylabel("AMD / AMD$_0$")
    ax.set_title("Angular-momentum deficit (secular-stability measure)")
    return _save(fig, path, cfg)


def fig_secular_spectrum(result: RunResult, cfg: Config, spectra: dict,
                         path: Path) -> Path:
    if not spectra.get("available") or not spectra.get("g"):
        return path
    theme = _theme(cfg); _apply(theme)
    inner = _present(result, constants.INNER_PLANETS)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for name in inner:
        sp = spectra["g"].get(name)
        if sp is None:
            continue
        freq, power = sp["freq"], sp["power"]
        mask = (freq > 0) & (freq < 60)
        p = power[mask]
        if p.size == 0:
            continue
        if p.max() > 0:
            p = p / p.max()
        ax.semilogy(freq[mask], np.clip(p, 1e-6, None), lw=1.1,
                    color=_color(name), label=name)
    for key, g in constants.SECULAR_G.items():
        ax.axvline(g, color=theme["muted"], ls=":", lw=0.8, alpha=0.7)
        ax.text(g, 1.15, key, rotation=90, fontsize=7, color=theme["muted"],
                ha="center", va="bottom")
    ax.set_xlim(0, 35); ax.set_ylim(1e-6, 2)
    ax.set_xlabel("frequency [arcsec / yr]"); ax.set_ylabel("normalised power")
    res = spectra.get("resolution_arcsec_yr", float("nan"))
    ax.set_title(f"Secular spectrum of e·e$^{{i\\varpi}}$  (FFT resolution ≈ {res:.2f}″/yr)")
    ax.legend(loc="upper right", fontsize=8)
    return _save(fig, path, cfg)


def fig_hk_phase(result: RunResult, cfg: Config, path: Path) -> Path:
    theme = _theme(cfg); _apply(theme)
    inner = _present(result, constants.INNER_PLANETS)
    fig, ax = plt.subplots(figsize=(7.5, 7))
    for name in inner:
        el = result.elements[name]
        h = el["e"] * np.cos(el["pomega"])
        k = el["e"] * np.sin(el["pomega"])
        ax.plot(h, k, lw=0.6, color=_color(name), alpha=0.85, label=name)
    ax.axhline(0, color=theme["grid"], lw=0.6); ax.axvline(0, color=theme["grid"], lw=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$h = e\cos\varpi$"); ax.set_ylabel(r"$k = e\sin\varpi$")
    ax.set_title("Secular phase portrait (eccentricity vector)")
    ax.legend(loc="upper right", fontsize=8)
    return _save(fig, path, cfg)


def fig_dashboard(result: RunResult, cfg: Config, spectra: dict, path: Path) -> Path:
    theme = _theme(cfg); _apply(theme)
    t, tlabel = _time_axis(result.times_years)
    inner = _present(result, constants.INNER_PLANETS)
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(3, 2, hspace=0.32, wspace=0.22)

    # (0,0) inner orbits (osculating ellipses, initial vs final epoch)
    ax = fig.add_subplot(gs[0, 0])
    _draw_orbits(ax, result, inner, theme, label=False)
    ax.set_title("Inner orbits (solid=init, dashed=final)"); ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")

    # (0,1) eccentricity
    ax = fig.add_subplot(gs[0, 1])
    for name in inner:
        ax.plot(t, result.elements[name]["e"], lw=0.9, color=_color(name), label=name)
    ax.set_title("Eccentricity"); ax.set_xlabel(tlabel); ax.set_ylabel("e")
    ax.legend(fontsize=7, ncol=2)

    # (1,0) conservation
    ax = fig.add_subplot(gs[1, 0])
    ax.semilogy(t, np.clip(result.dE_rel, 1e-18, None), color="#5ad1c2", lw=1, label="|ΔE/E|")
    ax.semilogy(t, np.clip(result.dL_rel, 1e-18, None), color="#c2a3ff", lw=1, label="|ΔL/L|")
    ax.set_title("Conservation"); ax.set_xlabel(tlabel); ax.legend(fontsize=7)

    # (1,1) MEGNO
    ax = fig.add_subplot(gs[1, 1])
    if result.megno is not None:
        ax.plot(t, result.megno, color="#ffb454", lw=1.1)
        ax.axhline(2.0, color="#5ad1c2", ls="--", lw=0.8)
        ax.set_title(f"MEGNO (final {result.megno[-1]:.2f})")
    else:
        ax.set_title("MEGNO (disabled)")
    ax.set_xlabel(tlabel); ax.set_ylabel("⟨Y⟩")

    # (2,0) inclination
    ax = fig.add_subplot(gs[2, 0])
    for name in inner:
        ax.plot(t, result.elements[name]["inc"] * RAD2DEG, lw=0.9, color=_color(name))
    ax.set_title("Inclination"); ax.set_xlabel(tlabel); ax.set_ylabel("i [deg]")

    # (2,1) secular spectrum
    ax = fig.add_subplot(gs[2, 1])
    if spectra.get("available") and spectra.get("g"):
        for name in inner:
            sp = spectra["g"].get(name)
            if sp is None:
                continue
            freq, power = sp["freq"], sp["power"]
            mask = (freq > 0) & (freq < 35)
            p = power[mask]
            if p.size == 0:
                continue
            if p.max() > 0:
                p = p / p.max()
            ax.semilogy(freq[mask], np.clip(p, 1e-6, None), lw=0.9, color=_color(name))
        for g in constants.SECULAR_G.values():
            ax.axvline(g, color=theme["muted"], ls=":", lw=0.7, alpha=0.6)
        ax.set_ylim(1e-5, 2)
    ax.set_title("Secular spectrum (g modes)"); ax.set_xlabel("arcsec/yr")

    status = result.meta.get("status", "")
    fig.suptitle(f"{cfg.run.name} — {result.meta['t_reached_years']:.3g} yr, "
                 f"GR={result.meta['gr_mode']}, dt={result.meta['dt_days']} d  [{status}]")
    return _save(fig, path, cfg)


def fig_animation(result: RunResult, cfg: Config, path: Path) -> Path | None:
    """Optional inner-system orbit animation (GIF via Pillow)."""
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except Exception:
        return None
    theme = _theme(cfg); _apply(theme)
    inner = _present(result, constants.INNER_PLANETS)
    if not inner:
        return None
    n = result.positions[inner[0]].shape[0]
    frames = min(cfg.plots.animate_frames, n)
    idx = np.linspace(0, n - 1, frames).astype(int)
    lim = max(np.max(np.abs(result.positions[m][:, :2])) for m in inner) * 1.1
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.scatter([0], [0], c=constants.SUN_COLOR, s=120, marker="*", zorder=5)
    dots = {m: ax.plot([], [], "o", color=_color(m), ms=6)[0] for m in inner}
    trails = {m: ax.plot([], [], "-", color=_color(m), lw=0.6, alpha=0.5)[0] for m in inner}

    def update(fi):
        i = idx[fi]
        for m in inner:
            xy = result.positions[m]
            dots[m].set_data([xy[i, 0]], [xy[i, 1]])
            trails[m].set_data(xy[:i + 1, 0], xy[:i + 1, 1])
        return list(dots.values()) + list(trails.values())

    anim = FuncAnimation(fig, update, frames=frames, blit=True)
    anim.save(str(path), writer=PillowWriter(fps=24))
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def make_all(result: RunResult, cfg: Config, spectra: dict, figdir: Path) -> list[tuple[str, Path]]:
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    made: list[tuple[str, Path]] = []

    def add(title, fname, fn, *args):
        out = fn(result, cfg, *args, figdir / fname)
        if out is not None and Path(out).exists():
            made.append((title, Path(out)))

    add("Orbital traces", "orbits.png", fig_orbits)
    add("Conservation", "conservation.png", fig_conservation)
    add("Eccentricity", "eccentricity.png", fig_eccentricity)
    add("Inclination", "inclination.png", fig_inclination)
    add("Semi-major axis", "semimajor.png", fig_semimajor)
    add("Body panel", "body_panel.png", fig_body_panel)
    add("AMD", "amd.png", fig_amd)
    add("h–k phase portrait", "hk_phase.png", fig_hk_phase)
    if result.megno is not None:
        add("MEGNO", "megno.png", fig_megno)
    add("Secular spectrum", "secular_spectrum.png", fig_secular_spectrum, spectra)
    add("Dashboard", "dashboard.png", fig_dashboard, spectra)

    if cfg.plots.animate:
        out = fig_animation(result, cfg, figdir / "orbit_animation.gif")
        if out is not None:
            made.append(("Orbit animation", Path(out)))
    return made
