"""CEAM Solar-System chaos & stability framework (REBOUND + REBOUNDx).

A configuration-driven framework for symplectic N-body integration of the
Solar System with general-relativistic corrections, chaos diagnostics (MEGNO /
Lyapunov), standard statistical outputs and a full set of figures.

Typical use::

    from solarsystem import load_config, run_pipeline
    cfg = load_config("config/default.yaml", overrides=["time.t_end_years=1e5"])
    out = run_pipeline(cfg)

See ``README.md`` for the parameter reference and experiment recipes.
"""

from __future__ import annotations

from .config import Config, load_config
from .builder import build_simulation, SimContext
from .runner import run, RunResult
from .pipeline import run_pipeline, PipelineOutput, format_console_summary

__all__ = [
    "Config", "load_config", "build_simulation", "SimContext",
    "run", "RunResult", "run_pipeline", "PipelineOutput",
    "format_console_summary",
]

__version__ = "0.1.0"
