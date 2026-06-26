"""Configuration model for a simulation run.

Every scientifically meaningful knob lives here so that experiments are *data*,
not code: vary Jupiter's mass, freeze a planet, toggle GR, change the timestep,
etc., all from ``config/default.yaml`` or ``--set a.b.c=value`` overrides.

The model is a set of nested dataclasses (so the available knobs are
self-documenting), built from a plain dict that is produced by deep-merging the
default YAML, an optional user YAML, and dotted command-line overrides.
"""

from __future__ import annotations

import copy
import dataclasses
import re
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

import yaml

# Matches a pure numeric literal, including sign-less exponents like "1.0e6"
# that PyYAML's resolver leaves as strings.
_NUMERIC_RE = re.compile(r"^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$")


# --------------------------------------------------------------------------- #
# Dataclasses describing the configuration tree
# --------------------------------------------------------------------------- #
@dataclass
class WHFastConfig:
    safe_mode: int = 0          # 0 = fast (must synchronize before reading state)
    corrector: int = 11         # 11th-order symplectic corrector (WHCKL)
    corrector2: int = 0
    kernel: str = "lazy"        # "default" | "lazy" (modified kernel -> WHCKL)


@dataclass
class IntegratorConfig:
    name: str = "whfast"        # whfast | ias15 | mercurius | trace | saba | whckl
    dt_days: float = 4.0        # production timestep; 1.0 for validation
    whfast: WHFastConfig = field(default_factory=WHFastConfig)
    exact_finish_time: int = 0  # 0 keeps whole-dt steps (uniform sampling for FFT)


@dataclass
class TimeConfig:
    t_end_years: float = 1.0e6
    n_outputs: int = 2000       # number of recorded snapshots (time series length)
    spacing: str = "linear"     # "linear" (needed for FFT) | "log"


@dataclass
class BodyOverride:
    """Per-body tweaks. ``None`` means "use the table / Horizons value"."""
    include: bool = True
    mass_scale: float = 1.0     # multiply the nominal mass (e.g. 2.0 -> 2 M_Jup)
    mass: float | None = None   # absolute mass in M_sun (overrides scale)
    a: float | None = None
    e: float | None = None
    inc: float | None = None    # degrees
    Omega: float | None = None  # degrees
    pomega: float | None = None # degrees
    l: float | None = None      # degrees (mean longitude)


@dataclass
class BodiesConfig:
    source: str = "builtin"     # "builtin" | "horizons" | "file"
    epoch: str = "2000-01-01 12:00"  # used for source=horizons
    file: str | None = None     # used for source=file (a cached body table YAML)
    include: list[str] = field(default_factory=lambda: [
        "Mercury", "Venus", "Earth", "Mars",
        "Jupiter", "Saturn", "Uranus", "Neptune"])
    overrides: dict[str, BodyOverride] = field(default_factory=dict)


@dataclass
class GRConfig:
    mode: str = "gr_potential"  # "none" | "gr_potential" | "gr" | "gr_full"
    sources: list[str] = field(default_factory=lambda: ["Sun"])
    c: float | None = None      # speed of light; None -> reboundx constant


@dataclass
class J2Config:
    enabled: bool = False
    J2: float = 2.2e-7
    R_sun_AU: float = 0.00465047


@dataclass
class PhysicsConfig:
    gr: GRConfig = field(default_factory=GRConfig)
    j2: J2Config = field(default_factory=J2Config)


@dataclass
class MegnoConfig:
    enabled: bool = True
    seed: int = 42


@dataclass
class ChaosConfig:
    megno: MegnoConfig = field(default_factory=MegnoConfig)


@dataclass
class FreezeConfig:
    """Experimental: hold listed bodies on their initial Keplerian orbit.

    Implemented via a heartbeat that resets the frozen bodies' states each
    step (simple, non-symplectic; see Deliverables.md "Frozen outer planets").
    Leave empty for a fully self-consistent integration.
    """
    bodies: list[str] = field(default_factory=list)


@dataclass
class ArchiveConfig:
    enabled: bool = True
    interval_years: float | None = None  # None -> match the output cadence


@dataclass
class DiagnosticsConfig:
    energy_tol: float = 1.0e-7       # pass/fail on max |dE/E|
    angmom_tol: float = 1.0e-10      # pass/fail on max |dL/L|
    compute_amd: bool = True


@dataclass
class PlotsConfig:
    enabled: bool = True
    dpi: int = 140
    style: str = "dark"             # "dark" | "light"
    max_orbit_points: int = 4000    # decimation for orbit-trace plots
    animate: bool = False
    animate_frames: int = 240


@dataclass
class RunConfig:
    name: str = "solar_default"
    output_root: str = "outputs"
    seed: int = 42


@dataclass
class Config:
    run: RunConfig = field(default_factory=RunConfig)
    integrator: IntegratorConfig = field(default_factory=IntegratorConfig)
    time: TimeConfig = field(default_factory=TimeConfig)
    bodies: BodiesConfig = field(default_factory=BodiesConfig)
    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    chaos: ChaosConfig = field(default_factory=ChaosConfig)
    freeze: FreezeConfig = field(default_factory=FreezeConfig)
    archive: ArchiveConfig = field(default_factory=ArchiveConfig)
    diagnostics: DiagnosticsConfig = field(default_factory=DiagnosticsConfig)
    plots: PlotsConfig = field(default_factory=PlotsConfig)

    # -- (de)serialisation ------------------------------------------------- #
    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Config":
        return _build_dataclass(cls, data or {})

    def save_yaml(self, path: str | Path) -> None:
        Path(path).write_text(yaml.safe_dump(self.to_dict(), sort_keys=False))


# --------------------------------------------------------------------------- #
# Helpers: dataclass<->dict, deep merge, dotted overrides
# --------------------------------------------------------------------------- #
def _build_dataclass(cls, data: Any):
    """Recursively construct dataclass ``cls`` from ``data`` (dict)."""
    if not is_dataclass(cls):
        return data
    kwargs = {}
    type_hints = {f.name: f.type for f in fields(cls)}
    for f in fields(cls):
        if f.name not in data:
            continue
        val = data[f.name]
        ftype = type_hints[f.name]
        # Nested dataclass
        if is_dataclass(_strip_optional(ftype)) and isinstance(val, dict):
            kwargs[f.name] = _build_dataclass(_strip_optional(ftype), val)
        # dict[str, BodyOverride]
        elif f.name == "overrides" and isinstance(val, dict):
            kwargs[f.name] = {k: _build_dataclass(BodyOverride, v)
                              for k, v in val.items()}
        else:
            kwargs[f.name] = val
    return cls(**kwargs)


def _strip_optional(ftype):
    """Best-effort resolve of a possibly-stringified / Optional annotation."""
    name_map = {
        "WHFastConfig": WHFastConfig, "IntegratorConfig": IntegratorConfig,
        "TimeConfig": TimeConfig, "BodiesConfig": BodiesConfig,
        "BodyOverride": BodyOverride, "GRConfig": GRConfig, "J2Config": J2Config,
        "PhysicsConfig": PhysicsConfig, "MegnoConfig": MegnoConfig,
        "ChaosConfig": ChaosConfig, "FreezeConfig": FreezeConfig,
        "ArchiveConfig": ArchiveConfig, "DiagnosticsConfig": DiagnosticsConfig,
        "PlotsConfig": PlotsConfig, "RunConfig": RunConfig,
    }
    if isinstance(ftype, str):
        return name_map.get(ftype, ftype)
    return ftype


def _normalize_numbers(obj: Any) -> Any:
    """Recursively coerce numeric-looking strings (e.g. ``"1.0e6"``) to numbers."""
    if isinstance(obj, dict):
        return {k: _normalize_numbers(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_normalize_numbers(v) for v in obj]
    if isinstance(obj, str) and _NUMERIC_RE.match(obj.strip()):
        s = obj.strip()
        try:
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return obj
    return obj


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge ``override`` into a copy of ``base``."""
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _coerce_scalar(text: str) -> Any:
    """Parse a CLI string into bool/int/float/None/list/dict.

    Handles scientific notation like ``2e4`` that YAML 1.1 leaves as a string.
    """
    t = text.strip()
    low = t.lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("null", "~", ""):   # NB: "none" is a valid value (e.g. gr.mode)
        return None
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)        # catches 2e4, 1e-3, .5, etc.
    except ValueError:
        pass
    try:
        return yaml.safe_load(t)   # lists/dicts/quoted strings
    except yaml.YAMLError:
        return t


def apply_dotted(data: dict, dotted: str) -> dict:
    """Apply a single ``a.b.c=value`` override to a nested dict (in place-ish)."""
    if "=" not in dotted:
        raise ValueError(f"Override '{dotted}' must be of the form key.path=value")
    key, _, raw = dotted.partition("=")
    keys = key.strip().split(".")
    value = _coerce_scalar(raw.strip())
    node = data
    for k in keys[:-1]:
        node = node.setdefault(k, {})
        if not isinstance(node, dict):
            raise ValueError(f"Cannot descend into non-mapping at '{k}' in '{dotted}'")
    node[keys[-1]] = value
    return data


def load_config(config_path: str | Path | None = None,
                overrides: list[str] | None = None,
                base: dict | None = None) -> Config:
    """Build a :class:`Config` from defaults + optional YAML + dotted overrides.

    Order of precedence (low -> high): dataclass defaults, ``base`` dict,
    ``config_path`` YAML, then each ``a.b.c=value`` in ``overrides``.
    """
    data = Config().to_dict()
    if base:
        data = deep_merge(data, base)
    if config_path:
        loaded = yaml.safe_load(Path(config_path).read_text()) or {}
        data = deep_merge(data, loaded)
    for ov in (overrides or []):
        apply_dotted(data, ov)
    data = _normalize_numbers(data)
    return Config.from_dict(data)
