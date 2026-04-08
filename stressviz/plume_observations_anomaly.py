#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from dateutil import parser as dtparser
from dateutil import tz

from astropy.time import Time
from astroquery.jplhorizons import Horizons


# ----------------------------
# Config: match encounters_io.py style
# ----------------------------
DATA_DIR = str(Path(__file__).resolve().parents[1] / "Data")
DEFAULT_PATH = os.path.join(DATA_DIR, "plume_observations.txt")


# ----------------------------
# Orbital math
# ----------------------------
def true_to_mean_anomaly_deg(nu_deg: float, e: float) -> float:
    """Convert true anomaly ν (deg) to mean anomaly M (deg) for an elliptic orbit."""
    nu = math.radians(float(nu_deg) % 360.0)
    e = float(e)

    if abs(e) < 1e-12:
        return float(nu_deg) % 360.0

    t = math.tan(nu / 2.0)
    fac = math.sqrt((1.0 - e) / (1.0 + e))
    E = 2.0 * math.atan2(fac * t, 1.0)
    M = E - e * math.sin(E)
    return math.degrees(M) % 360.0


# ----------------------------
# Time parsing
# ----------------------------
def _parse_start_datetime_utc(row: pd.Series) -> Optional[pd.Timestamp]:
    """
    Build start UTC datetime.
    Pref:
      1) Start (ISO)
      2) Start Date + Start Time (UTC)
    """
    iso = str(row.get("Start (ISO)", "")).strip()
    date_s = str(row.get("Start Date (mm/dd/yyyy)", "")).strip()
    time_s = str(row.get("Start Time (UTC)", "")).strip()

    if iso and iso.lower() != "nan":
        if not iso.endswith("T") and len(iso) >= 10:
            try:
                dt = dtparser.isoparse(iso)
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=tz.UTC)
                else:
                    dt = dt.astimezone(tz.UTC)
                return pd.Timestamp(dt)
            except Exception:
                pass

    if not date_s or date_s.lower() == "nan":
        return None
    if not time_s or time_s.lower() == "nan":
        return None

    try:
        dt = dtparser.parse(f"{date_s} {time_s}")
        dt = dt.replace(tzinfo=tz.UTC)  # treat as UTC
        return pd.Timestamp(dt)
    except Exception:
        return None


def _parse_exposure_s(row: pd.Series) -> Optional[float]:
    raw = str(row.get("Total Exposure (s)", "")).strip()
    if not raw or raw.lower() == "nan":
        return None
    try:
        return float(raw)
    except Exception:
        return None


def _mid_exposure_utc(row: pd.Series) -> Optional[pd.Timestamp]:
    start = _parse_start_datetime_utc(row)
    exp_s = _parse_exposure_s(row)
    if start is None or exp_s is None:
        return None
    return start + pd.to_timedelta(exp_s / 2.0, unit="s")


def _iso_z(ts: pd.Timestamp) -> str:
    dt = ts.to_pydatetime()
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=tz.UTC)
    else:
        dt = dt.astimezone(tz.UTC)
    return dt.isoformat().replace("+00:00", "Z")


# ----------------------------
# Horizons query
# ----------------------------
def _horizons_elements_europa(
    mid_utc: pd.Timestamp, debug: bool = False
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[str]]:
    """
    Query Horizons for Europa osculating elements.

    Key details:
      - For elements(), use TDB epochs
      - Pass scalar epoch (JD TDB) to avoid TLIST formatting issues
    """
    dt_utc = mid_utc.to_pydatetime()
    t = Time(dt_utc, scale="utc")
    jd_tdb = float(t.tdb.jd)

    obj = Horizons(
        id="502",            # Europa
        location="500@599",    # Jupiter system barycenter
        epochs=jd_tdb,       # scalar, not [jd]
        id_type=None,
    )

    if debug:
        payload = obj.elements(get_query_payload=True)
        print("[DEBUG] Horizons payload:", payload)

    el = obj.elements()
    row = el[0]
    cols = set(el.colnames)

    def pick(*names):
        for n in names:
            if n in cols:
                try:
                    return float(row[n])
                except Exception:
                    return None
        return None

    nu = pick("nu", "TA", "true_anom")
    M  = pick("M", "MA", "mean_anom")
    e  = pick("e", "EC", "ecc")

    return nu, M, e, ",".join(el.colnames)


# ----------------------------
# Script main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Precompute Europa ν + M for plume observations and update the TSV in place."
    )
    ap.add_argument(
        "--path",
        default=DEFAULT_PATH,
        help="Path to plume_observations.txt (default: Data/plume_observations.txt)",
    )
    ap.add_argument("--force", action="store_true", help="Recompute even if ν/M already present")
    ap.add_argument("--use-ma-from-horizons", action="store_true", help="Prefer Horizons-provided M when available")
    ap.add_argument("--dry-run", action="store_true", help="Do everything but do not write output")
    ap.add_argument("--debug-one", action="store_true", help="On first failure, print Horizons payload")
    args = ap.parse_args()

    path = os.path.abspath(os.path.expanduser(args.path))
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    print("[precompute] Using file:", path)

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    out_cols = ["Mid (ISO)", "true_anom_deg", "mean_anom_deg", "eccentricity", "phase_src"]
    for c in out_cols:
        if c not in df.columns:
            df[c] = ""

    n_total = len(df)
    n_ok = n_skip = n_fail = 0
    debug_used = False

    for i, row in df.iterrows():
        obs = str(row.get("Observer", "")).strip() or "UNKNOWN"

        if not args.force:
            have_nu = str(row.get("true_anom_deg", "")).strip()
            have_m  = str(row.get("mean_anom_deg", "")).strip()
            if have_nu and have_m:
                n_skip += 1
                continue

        mid = _mid_exposure_utc(row)
        if mid is None:
            n_skip += 1
            continue

        mid_iso = _iso_z(mid)
        df.at[i, "Mid (ISO)"] = mid_iso

        try:
            debug = bool(args.debug_one and (not debug_used))
            nu, M_h, e, colnames = _horizons_elements_europa(mid, debug=debug)
            if debug:
                debug_used = True

            if e is not None and np.isfinite(e):
                df.at[i, "eccentricity"] = f"{float(e):.12g}"

            if args.use_ma_from_horizons and (M_h is not None) and np.isfinite(M_h):
                M = float(M_h) % 360.0
                df.at[i, "mean_anom_deg"] = f"{M:.6f}"
                df.at[i, "phase_src"] = "horizons:M"
                if nu is not None and np.isfinite(nu):
                    df.at[i, "true_anom_deg"] = f"{float(nu)%360.0:.6f}"
                n_ok += 1
                continue

            if nu is None or e is None or (not np.isfinite(nu)) or (not np.isfinite(e)):
                raise RuntimeError(f"Missing nu/e from Horizons. Columns={colnames}")

            nu = float(nu) % 360.0
            e = float(e)
            M = true_to_mean_anomaly_deg(nu, e)

            df.at[i, "true_anom_deg"] = f"{nu:.6f}"
            df.at[i, "mean_anom_deg"] = f"{float(M)%360.0:.6f}"
            df.at[i, "phase_src"] = "horizons:nu->M"
            n_ok += 1

            if n_ok <= 5:
                print(f"[ok] row={i} {obs} mid={mid_iso} -> nu={nu:.6f} e={e:.6g} M={df.at[i,'mean_anom_deg']}")

        except Exception as e:
            n_fail += 1
            print(f"[FAIL] row={i} {obs} mid={mid_iso}: {e!r}")

    if n_ok == 0:
        print("[precompute] No successful Horizons results; NOT overwriting the file.")
        return

    print(f"[precompute] total={n_total} ok={n_ok} skip={n_skip} fail={n_fail}")

    if args.dry_run:
        print("[precompute] dry-run: not writing output")
        return

    tmp_fd, tmp_path = tempfile.mkstemp(
        dir=os.path.dirname(path),
        prefix="plume_observations_",
        suffix=".tmp",
        text=True,
    )
    os.close(tmp_fd)
    df.to_csv(tmp_path, sep="\t", index=False)
    shutil.move(tmp_path, path)

    print("[precompute] Updated file in place:", path)


if __name__ == "__main__":
    main()

