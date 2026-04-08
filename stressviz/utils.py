# stressviz/utils.py

import numpy as np
import math
try:
    import pandas as pd
except Exception:
    pd = None

def num_or_none(x):
    if x is None:
        return None
    if pd is not None and isinstance(x, (pd.Timestamp, type(getattr(pd, "NaT", None)))):
        return None
    try:
        v = float(x)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except Exception:
        try:
            v = float(str(x))
            return None if (math.isnan(v) or math.isinf(v)) else v
        except Exception:
            return None

def deg2rad(deg):
    return np.radians(deg)

def rad2deg(rad):
    return np.degrees(rad)

def normalize_longitude(lon):
    """Wrap longitude to [-180, 180] or [0, 360] as needed."""
    return (lon + 360) % 360

# --- SatStress import shim (handles SatStress.satstress layout) ---
import importlib
import pkgutil

def get_satstress_handles():
    """
    Return (pkg, StressCalc, Diurnal) regardless of how the repo exposes them.
    Tries SatStress.*, then satstress.*, and scans submodules (incl. satstress).
    """
    # 1) Load package (capitalized first, then lowercase)
    pkg = None
    for name in ("SatStress", "satstress"):
        try:
            pkg = importlib.import_module(name)
            break
        except Exception:
            continue
    if pkg is None:
        raise ImportError("Could not import SatStress/satstress (is it on sys.path?)")

    StressCalc = getattr(pkg, "StressCalc", None)
    Diurnal    = getattr(pkg, "Diurnal", None)

    # 2) Try common submodules, including the IMPORTANT 'satstress' module
    submods_to_try = [
        f"{pkg.__name__}.satstress",     # <-- many forks put classes here
        f"{pkg.__name__}.gridcalc",
        f"{pkg.__name__}.stresscalc",
        f"{pkg.__name__}.grid_calc",
        f"{pkg.__name__}.diurnal",
        f"{pkg.__name__}.stressdef",
        f"{pkg.__name__}.eccentric",
    ]
    for modname in submods_to_try:
        try:
            m = importlib.import_module(modname)
        except Exception:
            continue
        StressCalc = StressCalc or getattr(m, "StressCalc", None)
        Diurnal    = Diurnal    or getattr(m, "Diurnal", None)
        if StressCalc and Diurnal:
            break

    # 3) Still missing? Walk all submodules and look for the classes
    if (StressCalc is None or Diurnal is None) and hasattr(pkg, "__path__"):
        for finder, name, ispkg in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            try:
                m = importlib.import_module(name)
            except Exception:
                continue
            if StressCalc is None and hasattr(m, "StressCalc"):
                StressCalc = getattr(m, "StressCalc")
            if Diurnal is None and hasattr(m, "Diurnal"):
                Diurnal = getattr(m, "Diurnal")
            if StressCalc and Diurnal:
                break

    if StressCalc is None or Diurnal is None:
        missing = []
        if StressCalc is None: missing.append("StressCalc")
        if Diurnal is None:    missing.append("Diurnal")
        raise ImportError(f"Found {pkg.__name__} but missing: {', '.join(missing)}")

    return pkg, StressCalc, Diurnal


# --- Resolve specific stress-definition classes (Diurnal/NSR/Secular) and build calc ---

def _find_class_in_pkg(pkg, class_name: str):
    """Best-effort search for a class inside the SatStress package tree."""
    import importlib, pkgutil
    tried = []

    # common candidates first
    candidate_mods = [
        pkg.__name__,                         # top-level
        f"{pkg.__name__}.satstress",          # common fork layout
        f"{pkg.__name__}.stress",             # some forks
        f"{pkg.__name__}.stressdef",
        f"{pkg.__name__}.diurnal",
        f"{pkg.__name__}.nsr",
        f"{pkg.__name__}.secular",
        f"{pkg.__name__}.eccentric",
    ]

    for modname in candidate_mods:
        try:
            m = importlib.import_module(modname)
            if hasattr(m, class_name):
                return getattr(m, class_name)
        except Exception as e:
            tried.append(f"{modname}: {e}")

    # fall back: scan full package
    if hasattr(pkg, "__path__"):
        for _, modname, _ in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            try:
                m = importlib.import_module(modname)
                if hasattr(m, class_name):
                    return getattr(m, class_name)
            except Exception as e:
                tried.append(f"{modname}: {e}")

    raise ImportError(
        f"Could not find class '{class_name}' in package '{pkg.__name__}'. "
        "Tried common modules and walked package."
    )


def resolve_stresscalc(kind: str = "diurnal"):
    """
    Return the SatStress *stress-definition class* for the requested kind.
    Example:
        CalcClass = resolve_stresscalc("diurnal")  # -> Diurnal class
        calc = CalcClass(satellite)                # if that class takes a satellite
    """
    kind_l = (kind or "").lower()
    pkg, _StressCalc, _Diurnal = get_satstress_handles()

    if kind_l in ("diurnal", "eccentric"):
        # Prefer the already-located Diurnal if available
        try:
            return _find_class_in_pkg(pkg, "Diurnal")
        except Exception:
            return _Diurnal

    if kind_l in ("nsr", "non-synchronous", "nonsynchronous", "non_synchronous"):
        return _find_class_in_pkg(pkg, "NSR")

    if kind_l in ("secular",):
        return _find_class_in_pkg(pkg, "Secular")

    raise ValueError(f"Unknown stress kind '{kind}'. Expected 'diurnal', 'nsr', or 'secular'.")


def build_stresscalc(kind: str, satellite):
    """
    Convenience: construct a StressCalc configured with the requested stress-definition.
    Returns a StressCalc instance, trying common constructor signatures.

    Usage:
        sc = build_stresscalc("diurnal", sat)
        # sc now has sat bound and includes the Diurnal calc
    """
    pkg, StressCalc, _ = get_satstress_handles()
    StressClass = resolve_stresscalc(kind)

    # mirror _build_stresscalc but with the resolved class
    for attempt in (
        lambda: StressCalc([StressClass], satellite=satellite),
        lambda: StressCalc([StressClass]),
        lambda: StressCalc(StressClass, satellite=satellite),
        lambda: StressCalc(StressClass),
    ):
        try:
            sc = attempt()
            # Bind satellite if constructor didn’t
            if not getattr(sc, "satellite", None):
                sc.satellite = satellite
            return sc
        except Exception:
            continue

    # last resort: construct bare then patch in
    try:
        sc = StressCalc([StressClass])
        if not getattr(sc, "satellite", None):
            sc.satellite = satellite
        return sc
    except Exception as e:
        raise TypeError(f"Could not construct StressCalc for kind '{kind}': {e}")

def true_to_mean_anomaly_deg(nu_deg: float, e: float) -> float:
    """Convert true anomaly (deg) -> mean anomaly (deg) for 0 <= e < 1; returns [0,360)."""
    if nu_deg is None:
        return None
    nu = math.radians(nu_deg % 360.0)
    denom = 1.0 + e * math.cos(nu)
    if abs(denom) < 1e-15:
        return float(nu_deg % 360.0)
    cosE = (e + math.cos(nu)) / denom
    sinE = (math.sqrt(max(0.0, 1 - e*e)) * math.sin(nu)) / denom
    E = math.atan2(sinE, cosE)
    M = E - e * math.sin(E)
    return (math.degrees(M) % 360.0)