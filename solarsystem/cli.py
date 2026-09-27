"""Command-line interface: ``ceam-sim``.

Examples
--------
  ceam-sim                                   # default 1-Myr run -> outputs/solar_default
  ceam-sim --quick --name smoke              # fast pipeline check
  ceam-sim --t-end 5e6 --name mercury_5Myr   # longer run
  ceam-sim --set bodies.overrides.Jupiter.mass_scale=2.0 --name jup2x
  ceam-sim --gr none --name no_gr            # GR off (chaos-sensitivity experiment)
  ceam-sim horizons --epoch "2025-01-01 00:00" --out cache/horizons_2025.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

DEFAULT_CONFIG = "config/default.yaml"


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", default=None,
                   help=f"YAML config (default: {DEFAULT_CONFIG} if present)")
    p.add_argument("--set", dest="overrides", action="append", default=[],
                   metavar="a.b.c=value", help="override any config key (repeatable)")
    p.add_argument("--name", default=None, help="run name (output subdirectory)")
    p.add_argument("--t-end", type=float, default=None, help="end time in years")
    p.add_argument("--dt", type=float, default=None, help="timestep in days")
    p.add_argument("--gr", default=None,
                   choices=["none", "gr_potential", "gr", "gr_full"], help="GR mode")
    p.add_argument("--no-megno", action="store_true", help="disable MEGNO")
    p.add_argument("--no-plots", action="store_true", help="skip figure generation")
    p.add_argument("--no-progress", action="store_true", help="hide progress bar")
    p.add_argument("--animate", action="store_true", help="render orbit animation GIF")
    p.add_argument("--quick", action="store_true",
                   help="short low-res run for a fast pipeline check")


def _overrides_from_args(args) -> list[str]:
    ov = list(args.overrides)
    if args.quick:
        ov += ["time.t_end_years=20000", "time.n_outputs=500"]
    if args.name is not None:
        ov.append(f"run.name={args.name}")
    if args.t_end is not None:
        ov.append(f"time.t_end_years={args.t_end}")
    if args.dt is not None:
        ov.append(f"integrator.dt_days={args.dt}")
    if args.gr is not None:
        ov.append(f"physics.gr.mode={args.gr}")
    if args.no_megno:
        ov.append("chaos.megno.enabled=false")
    if args.no_plots:
        ov.append("plots.enabled=false")
    if args.animate:
        ov.append("plots.animate=true")
    return ov


def _resolve_config_path(arg) -> str | None:
    if arg:
        return arg
    return DEFAULT_CONFIG if Path(DEFAULT_CONFIG).exists() else None


def cmd_run(args) -> int:
    from .config import load_config
    from .pipeline import run_pipeline, format_console_summary

    cfg = load_config(_resolve_config_path(args.config), _overrides_from_args(args))
    print(f"[ceam-sim] starting run '{cfg.run.name}' "
          f"(t_end={cfg.time.t_end_years:g} yr, dt={cfg.integrator.dt_days} d, "
          f"GR={cfg.physics.gr.mode})")
    out = run_pipeline(cfg, progress=not args.no_progress)
    print("\n[ceam-sim] done.\n" + format_console_summary(out))
    return 0


def cmd_horizons(args) -> int:
    from . import constants
    from .horizons import fetch_body_table, save_body_table

    names = args.bodies or constants.PLANET_ORDER
    print(f"[ceam-sim] querying JPL Horizons for {names} @ {args.epoch} ...")
    table = fetch_body_table(names, epoch=args.epoch, earth_moon=args.earth_moon)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    save_body_table(table, args.out)
    print(f"[ceam-sim] wrote body table -> {args.out}")
    print("           use it with:  --set bodies.source=file "
          f"--set bodies.file={args.out}")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="ceam-sim",
        description="CEAM Solar-System chaos & stability framework (REBOUND + REBOUNDx). "
                    "With no subcommand, runs the simulation pipeline.")
    sub = parser.add_subparsers(dest="command")
    p_hor = sub.add_parser("horizons", help="fetch initial conditions from JPL Horizons")
    p_hor.add_argument("--epoch", default="2000-01-01 12:00")
    p_hor.add_argument("--out", default="cache/horizons_table.yaml")
    p_hor.add_argument("--bodies", nargs="*", default=None)
    p_hor.add_argument("--earth-moon", default="barycenter",
                       choices=["barycenter", "geocenter"])
    _add_common(parser)  # run options at the top level (so bare + --help show them)

    if argv and argv[0] == "run":     # "run" is optional sugar for the default action
        argv = argv[1:]
    if argv and argv[0] == "horizons":
        return cmd_horizons(parser.parse_args(argv))

    args = parser.parse_args(argv)
    if getattr(args, "command", None) == "horizons":
        return cmd_horizons(args)
    return cmd_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
