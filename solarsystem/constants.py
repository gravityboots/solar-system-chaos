"""Physical constants and default Solar-System body data.

All simulation work is done in REBOUND's ``('AU', 'yr2pi', 'Msun')`` unit system:
G = 1, lengths in AU, masses in solar masses, and time in years/(2*pi) so that
one Julian year = 2*pi time units.  In these units the speed of light is
``reboundx.constants.C = 10065.32`` AU/(yr/2pi); see :data:`SPEED_OF_LIGHT`.

The default body table holds masses (M_sun) and J2000 heliocentric ecliptic
osculating elements (the Standish "approximate positions" set), which lets the
framework build the whole system *offline* and lets every element be perturbed
from config.  Use :mod:`solarsystem.horizons` to refresh these from JPL Horizons.
"""

from __future__ import annotations

import math

# --- Units --------------------------------------------------------------------
UNITS = ("AU", "yr2pi", "Msun")
TWO_PI = 2.0 * math.pi
DAYS_PER_YEAR = 365.25
#: Speed of light in AU / (yr/2pi).  Matches reboundx.constants.C exactly.
SPEED_OF_LIGHT = 10065.32012038219
#: Solar radius in AU (for the optional J2 oblateness term).
R_SUN_AU = 0.00465047
#: Sun's gravitational oblateness coefficient (Mecheri et al. 2004).
J2_SUN = 2.2e-7
ARCSEC_PER_RAD = 180.0 / math.pi * 3600.0


def days_to_code_time(dt_days: float) -> float:
    """Convert a timestep in days to REBOUND yr2pi time units."""
    return dt_days / DAYS_PER_YEAR * TWO_PI


def years_to_code_time(years: float) -> float:
    """Convert a duration in Julian years to REBOUND yr2pi time units."""
    return years * TWO_PI


def code_time_to_years(t: float) -> float:
    """Convert REBOUND yr2pi time units back to Julian years."""
    return t / TWO_PI


# --- Default body table -------------------------------------------------------
# Masses are 1/(JPL reciprocal mass).  Earth is the Earth-Moon barycentre, the
# standard choice for inner-planet chaos work (Laskar & Gastineau 2009).
# Elements: a [AU], e, inc [deg], Omega [deg] (long. ascending node),
#           pomega [deg] (long. of perihelion), l [deg] (mean longitude),
# in the J2000 ecliptic frame.  REBOUND consumes (a, e, inc, Omega, pomega, l).

SUN_MASS = 1.0

# name -> dict
DEFAULT_BODIES: dict[str, dict] = {
    "Mercury": dict(mass=1.0 / 6023600.0, a=0.38709927, e=0.20563593, inc=7.00497902,
                    Omega=48.33076593, pomega=77.45779628, l=252.25032350,
                    group="inner", color="#b1b1b1"),
    "Venus":   dict(mass=1.0 / 408523.71, a=0.72333566, e=0.00677672, inc=3.39467605,
                    Omega=76.67984255, pomega=131.60246718, l=181.97909950,
                    group="inner", color="#e6c27a"),
    "Earth":   dict(mass=1.0 / 328900.5614, a=1.00000261, e=0.01671123, inc=-0.00001531,
                    Omega=0.0, pomega=102.93768193, l=100.46457166,
                    group="inner", color="#3b82c4"),
    "Mars":    dict(mass=1.0 / 3098708.0, a=1.52371034, e=0.09339410, inc=1.84969142,
                    Omega=49.55953891, pomega=-23.94362959, l=-4.55343205,
                    group="inner", color="#c1440e"),
    "Jupiter": dict(mass=1.0 / 1047.3486, a=5.20288700, e=0.04838624, inc=1.30439695,
                    Omega=100.47390909, pomega=14.72847983, l=34.39644051,
                    group="outer", color="#d8a06b"),
    "Saturn":  dict(mass=1.0 / 3497.898, a=9.53667594, e=0.05386179, inc=2.48599187,
                    Omega=113.66242448, pomega=92.59887831, l=49.95424423,
                    group="outer", color="#e3d9a0"),
    "Uranus":  dict(mass=1.0 / 22902.98, a=19.18916464, e=0.04725744, inc=0.77263783,
                    Omega=74.01692503, pomega=170.95427630, l=313.23810451,
                    group="outer", color="#9fd8e0"),
    "Neptune": dict(mass=1.0 / 19412.24, a=30.06992276, e=0.00859048, inc=1.77004347,
                    Omega=131.78422574, pomega=44.96476227, l=-55.12002969,
                    group="outer", color="#4f6fd8"),
}

#: Canonical planet ordering used everywhere (Sun is implicit index 0).
PLANET_ORDER = ["Mercury", "Venus", "Earth", "Mars",
                "Jupiter", "Saturn", "Uranus", "Neptune"]

INNER_PLANETS = [n for n in PLANET_ORDER if DEFAULT_BODIES[n]["group"] == "inner"]
OUTER_PLANETS = [n for n in PLANET_ORDER if DEFAULT_BODIES[n]["group"] == "outer"]

SUN_COLOR = "#ffd24d"

# --- Reference secular frequencies (Laskar et al. 2004), arcsec / yr ----------
# g_i are perihelion precession eigenmodes; s_i are nodal regression eigenmodes.
# Peaks in the FFT of e*exp(i*pomega) should fall near g_i; peaks in
# sin(inc/2)*exp(i*Omega) near s_i.  Used to annotate the secular spectra.
SECULAR_G = {  # arcsec/yr
    "g1": 5.59, "g2": 7.452, "g3": 17.368, "g4": 17.916,
    "g5": 4.257, "g6": 28.245, "g7": 3.087, "g8": 0.673,
}
SECULAR_S = {  # arcsec/yr (s1 ~ 0, s5 ~ 0 by symmetry / definition)
    "s1": -5.61, "s2": -7.06, "s3": -18.848, "s4": -17.751,
    "s5": 0.0, "s6": -26.347, "s7": -2.993, "s8": -0.692,
}

#: The headline chaos-driving secular resonance (Laskar 1989/2008).
RESONANCE_G1_MINUS_G5 = SECULAR_G["g1"] - SECULAR_G["g5"]  # ~1.33 arcsec/yr

#: Reference inner-system Lyapunov time (Laskar) for sanity-checking, in years.
REFERENCE_LYAPUNOV_TIME_YEARS = 5.0e6
