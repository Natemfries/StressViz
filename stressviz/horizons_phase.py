# stressviz/horizons_phase.py
"""
Horizons-backed ephemeris helpers for StressViz.

Primary use case:
  - Given a UTC ISO timestamp, return Europa's true anomaly (nu) in its orbit
    about Jupiter, along with eccentricity (and optionally mean anomaly).

Notes:
  - Horizons elements() prefers epochs in TDB. We convert UTC -> JD(TDB) and pass
    a scalar epoch to avoid TLIST formatting issues.
  - For "Europa around Jupiter", use location="500@599" (Jupiter center).
    If you specifically want the Jupiter system barycenter, use "500@5".
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Union, Iterable, List, Tuple

import numpy as np

try:
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None  # type: ignore

from astropy.time import Time
from astroquery.jplhorizons import Horizons


# ---------- public dataclass ----------

@dataclass(frozen=True)
class BodyPhase:
    utc_iso: str
    target: str
    center: str
    nu_deg: float
    e: float
    M_deg: Optional[float] = None
    colnames: Optional[str] = None

    @property
    def nu_rad(self) -> float:
        return float(np.deg2rad(self.nu_deg) % (2.0 * np.pi))

    @property
    def M_rad(self) -> Optional[float]:
        if self.M_deg is None:
            return None
        return float(np.deg2rad(self.M_deg) % (2.0 * np.pi))


# Backward-compatible alias
EuropaPhase = BodyPhase


HORIZONS_BODY_IDS = {
    "EUROPA": "502",
    "GANYMEDE": "503",
    "ENCELADUS": "602",
}


HORIZONS_CENTER_IDS = {
    "JUPITER": "599",
    "SATURN": "699",
}


def _horizons_id_for_target(target: str) -> str:
    t = str(target).strip().upper()
    return HORIZONS_BODY_IDS.get(t, t)


def _location_for_center(center: str) -> str:
    c = center.strip().upper()

    if c in {"JUPITER", "599"}:
        return "500@599"

    if c in {"SATURN", "699"}:
        return "500@699"

    if c in {"JSB", "JUPITER_BARYCENTER", "BARYCENTER", "5"}:
        return "500@5"

    raise ValueError(
        f"Unknown center={center!r}. "
        "Use 'JUPITER' (599), 'SATURN' (699), or 'JSB' (5)."
    )


# ---------- internal helpers ----------

def _normalize_utc_iso(utc: Union[str, "pd.Timestamp"]) -> str:
    if pd is not None and isinstance(utc, pd.Timestamp):
        # keep it simple; your project may already have an _iso_z() helper
        s = utc.isoformat()
        return s.replace("+00:00", "Z")
    return str(utc).strip()


def _utc_to_jd_tdb_scalar(utc: Union[str, "pd.Timestamp"]) -> Tuple[str, float]:
    """
    Convert UTC time -> scalar JD in TDB, suitable for Horizons(elements).
    """
    utc_iso = _normalize_utc_iso(utc)
    # Let astropy parse. This handles "YYYY-mm-ddTHH:MM:SSZ" and similar.
    t = Time(utc_iso, scale="utc")
    jd_tdb = float(t.tdb.jd)
    return utc_iso, jd_tdb


def _pick_float(row, cols, *names) -> Optional[float]:
    for n in names:
        if n in cols:
            try:
                return float(row[n])
            except Exception:
                return None
    return None


# ---------- public API (single epoch) ----------

@lru_cache(maxsize=4096)
def orbital_phase_from_horizons(
    utc_iso: str,
    *,
    target: str = "EUROPA",
    center: str = "JUPITER",
    prefer_M_from_horizons: bool = False,
    debug_payload: bool = False,
) -> BodyPhase:
    """
    Return target body's true anomaly about the requested center at a UTC time.

    Examples
    --------
    target="EUROPA", center="JUPITER"       -> id="502", location="500@599"
    target="ENCELADUS", center="SATURN"    -> id="602", location="500@699"
    """
    utc_norm, jd_tdb = _utc_to_jd_tdb_scalar(utc_iso)
    location = _location_for_center(center)
    target_id = _horizons_id_for_target(target)

    obj = Horizons(
        id=target_id,
        location=location,
        epochs=jd_tdb,
        id_type=None,
    )

    if debug_payload:
        payload = obj.elements(get_query_payload=True)
        print("[DEBUG] Horizons payload:", payload)

    el = obj.elements()
    row = el[0]
    cols = set(el.colnames)

    nu = _pick_float(row, cols, "nu", "TA", "true_anom")
    e  = _pick_float(row, cols, "e", "EC", "ecc")
    Mh = _pick_float(row, cols, "M", "MA", "mean_anom")

    if nu is None or e is None or (not np.isfinite(nu)) or (not np.isfinite(e)):
        raise RuntimeError(f"Missing nu/e from Horizons. Columns={','.join(el.colnames)}")

    nu_deg = float(nu) % 360.0
    e_val = float(e)

    M_deg: Optional[float] = None
    if prefer_M_from_horizons and (Mh is not None) and np.isfinite(Mh):
        M_deg = float(Mh) % 360.0

    return BodyPhase(
        utc_iso=utc_norm,
        target=str(target).strip().upper(),
        center=str(center).strip().upper(),
        nu_deg=nu_deg,
        e=e_val,
        M_deg=M_deg,
        colnames=",".join(el.colnames),
    )


# ---------- optional convenience (batch epochs) ----------

def europa_nu_series_from_horizons(
    utc_isos: Iterable[str],
    *,
    center: str = "JUPITER",
    debug_payload: bool = False,
) -> List[EuropaPhase]:
    """
    Batch query Horizons for many epochs at once.

    Returns a list of EuropaPhase (nu_deg, e). M_deg is not included here by default
    because it varies by Horizons column availability; you can add it similarly if needed.
    """
    utc_list = [str(s).strip() for s in utc_isos]
    if len(utc_list) == 0:
        return []

    # Astropy can parse a list of ISO strings directly.
    t = Time(utc_list, scale="utc")
    jd_tdb_list = [float(x) for x in t.tdb.jd]

    location = _location_for_center(center)
    obj = Horizons(
        id="502",
        location=location,
        epochs=jd_tdb_list,
        id_type=None,
    )

    if debug_payload:
        payload = obj.elements(get_query_payload=True)
        print("[DEBUG] Horizons payload:", payload)

    el = obj.elements()
    cols = set(el.colnames)

    # vector columns
    def col_as_float(name_candidates: Tuple[str, ...]) -> Optional[np.ndarray]:
        for n in name_candidates:
            if n in cols:
                try:
                    return np.array(el[n], dtype=float)
                except Exception:
                    return None
        return None

    nu_arr = col_as_float(("nu", "TA", "true_anom"))
    e_arr  = col_as_float(("e", "EC", "ecc"))

    if nu_arr is None or e_arr is None:
        raise RuntimeError(f"Missing nu/e from Horizons. Columns={','.join(el.colnames)}")

    phases: List[EuropaPhase] = []
    for utc_iso, nu, e in zip(utc_list, nu_arr, e_arr):
        if (not np.isfinite(nu)) or (not np.isfinite(e)):
            continue
        phases.append(
            EuropaPhase(
                utc_iso=_normalize_utc_iso(utc_iso),
                nu_deg=float(nu) % 360.0,
                e=float(e),
                M_deg=None,
                colnames=",".join(el.colnames),
            )
        )

    return phases
