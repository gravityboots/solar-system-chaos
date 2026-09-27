"""End-to-end pipeline: build → integrate → analyse → visualise → report."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import analysis, builder, report, runner, visualize
from .config import Config


@dataclass
class PipelineOutput:
    outdir: Path
    result: runner.RunResult
    summary: dict
    spectra: dict
    figures: list
    report_path: Path | None


def run_pipeline(cfg: Config, progress: bool = True,
                 make_plots: bool | None = None) -> PipelineOutput:
    outdir = Path(cfg.run.output_root) / cfg.run.name
    outdir.mkdir(parents=True, exist_ok=True)
    cfg.save_yaml(outdir / "config.resolved.yaml")

    # 1. Build the simulation (bodies + integrator + physics + MEGNO).
    ctx = builder.build_simulation(cfg)

    # 2. Integrate, recording diagnostics + SimulationArchive.
    result = runner.run(ctx, outdir, progress=progress)
    result.save_npz(outdir / "timeseries.npz")

    # 3. Analyse: statistics, conservation, chaos, secular spectra.
    summary = analysis.summarize(result, cfg)
    spectra = analysis.secular_spectra(result, cfg)
    analysis.save_statistics(summary, outdir)

    # 4. Visualise.
    figures: list = []
    do_plots = cfg.plots.enabled if make_plots is None else make_plots
    if do_plots:
        figures = visualize.make_all(result, cfg, spectra, outdir / "figures")

    # 5. Report.
    report_path = report.write_report(result, cfg, summary, spectra, figures, outdir)

    return PipelineOutput(outdir=outdir, result=result, summary=summary,
                          spectra=spectra, figures=figures, report_path=report_path)


def format_console_summary(out: PipelineOutput) -> str:
    s = out.summary
    cons, chaos, meta = s["conservation"], s["chaos"], s["meta"]
    lines = [
        f"  run            : {meta['run_name']}  ({meta['status']})",
        f"  integrated to  : {meta['t_reached_years']:.4g} yr "
        f"({meta['n_samples']} samples, dt={meta['dt_days']} d, GR={meta['gr_mode']})",
        f"  max |dE/E|     : {cons['max_dE_rel']:.3e}  "
        f"({'ok' if cons['energy_ok'] else 'OVER tol'})",
        f"  max |dL/L|     : {cons['max_dL_rel']:.3e}  "
        f"({'ok' if cons['angmom_ok'] else 'OVER tol'})",
    ]
    if chaos.get("megno_enabled"):
        lines.append(f"  MEGNO (final)  : {chaos['megno_final']:.3f} "
                     f"-> {chaos['classification']}")
    lines.append(f"  outputs        : {out.outdir}")
    if out.report_path:
        lines.append(f"  report         : {out.report_path}")
    return "\n".join(lines)
