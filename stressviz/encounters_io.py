# stressviz/encounters_io.py
import os
import re
import pandas as pd
from pathlib import Path
from typing import Dict, Callable, Any, List
from .encounters import EncounterRegistry
from .utils import true_to_mean_anomaly_deg as _true2mean

DATA_DIR = str(Path(__file__).resolve().parents[1] / "Data")
ENCOUNTERS_FILE = "Europa_Encounters_21F31_V7_LP01_ver2.txt"

# zero-based indices: 0(body_id), 1(enc_no), 3(ET sec), 11(lon E), 12(lat pc), 14(true anomaly)
_USECOLS = [0, 1, 3, 11, 12, 14]
_COLNAMES = ["body_id", "enc_no", "tca_et_sec", "lon_e_deg", "lat_pc_deg", "true_anom_deg"]

def load_europa_min(path: str | None = None) -> pd.DataFrame:
    if path is None:
        path = os.path.join(DATA_DIR, ENCOUNTERS_FILE)

    df = pd.read_csv(
        path,
        comment="%",             # skip header comments
        sep=",",
        skipinitialspace=True,
        header=0,                # first non-% line is header
        usecols=_USECOLS,
        dtype=str,               # avoid auto Timestamps
        parse_dates=False,
        engine="python",
    )
    df.columns = _COLNAMES

    # Europa only
    df = df[df["body_id"].str.strip() == "502"].copy()

    # Preserve zero-padded label and numeric sort key
    df["enc_label"]  = df["enc_no"].astype(str).str.strip()
    df["enc_no_num"] = pd.to_numeric(df["enc_no"], errors="coerce")

    # Explicit numeric coercion
    for c in ("tca_et_sec", "lon_e_deg", "lat_pc_deg", "true_anom_deg"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Optional pretty string for UI
    try:
        from .utils import format_et_seconds_to_display
        df["tca_tdb_str"] = df["tca_et_sec"].map(format_et_seconds_to_display)
    except Exception:
        df["tca_tdb_str"] = ""

    # Sorted by encounter number
    df = df.sort_values("enc_no_num", kind="stable", ignore_index=True)
    return df[["enc_label","enc_no_num","tca_et_sec","true_anom_deg","lat_pc_deg","lon_e_deg","tca_tdb_str"]]

def build_plume_rows_from_excel(xlsx_path: str, horizons_lookup: Callable[[str], dict]):
    df = pd.read_excel(xlsx_path)

    rows = []
    for _, row in df.iterrows():
        observer = str(row["Observer"]).strip()

        # Use your existing Start ISO column or recompute it
        utc_iso = str(row.get("Start ISO") or "").strip()
        if not utc_iso:
            continue

        # Get Horizons anomaly values
        h = horizons_lookup(utc_iso)

        # Build encounter rows the way EncounterRegistry expects
        safe_time = utc_iso.replace("-", "").replace(":", "").replace("T", "_")
        enc_id = f"PLUME_{observer}_{safe_time}"

        rows.append({
            "encounter_id": enc_id,
            "encounter_tag": "plume",
            "et_tdb_sec": h["et_tdb_sec"],
            "utc_iso": utc_iso,
            "true_anom_deg": h.get("true_anom_deg"),
            "mean_anom_deg": h.get("mean_anom_deg"),
            "lat_deg": None,
            "lon_deg": None,
            "period_hours": None,
        })
    return rows


ISO_UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")

def build_plume_rows_from_txt(txt_path: str) -> list[dict]:
    """
    Read plume observations from a tab-delimited text file and return rows
    shaped the way EncounterRegistry expects.

    This does NOT query Horizons/SPICE. It normalizes input and creates rows.

    NEW: propagates any precomputed phase columns if present in the TSV:
      - Mid (ISO)
      - true_anom_deg
      - mean_anom_deg
      - eccentricity
      - phase_src
      - period_hours (optional)
      - detection
    """
    import re
    import pandas as pd

    def _to_float_or_none(x):
        s = str(x).strip()
        if not s or s.lower() == "nan":
            return None
        try:
            return float(s)
        except Exception:
            return None

    def _to_str_or_none(x):
        s = str(x).strip()
        if not s or s.lower() == "nan":
            return None
        return s

    df = pd.read_csv(txt_path, sep="\t", dtype=str).fillna("")

    # Normalize the canonical columns we rely on
    df = df.rename(columns={
        "Start (ISO)": "utc_iso",
        "Observer": "observer",
        "Total Exposure (s)": "exposure_s",
    })

    for col in ("utc_iso", "observer"):
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in {txt_path}")

    # Optional columns (only present after precompute)
    has_mid = "Mid (ISO)" in df.columns
    has_ta  = "true_anom_deg" in df.columns
    has_ma  = "mean_anom_deg" in df.columns
    has_ecc = "eccentricity" in df.columns
    has_src = "phase_src" in df.columns
    has_P   = "period_hours" in df.columns
    has_det = "detection" in df.columns

    rows: list[dict] = []

    for _, r in df.iterrows():
        utc_iso = str(r.get("utc_iso", "")).strip()
        observer = str(r.get("observer", "UNKNOWN")).strip() or "UNKNOWN"

        # Optional exposure (float or None)
        exposure_s = _to_float_or_none(r.get("exposure_s", ""))

        # Pull precomputed values if present
        mid_iso = _to_str_or_none(r.get("Mid (ISO)", "")) if has_mid else None
        ta_deg  = _to_float_or_none(r.get("true_anom_deg", "")) if has_ta else None
        ma_deg  = _to_float_or_none(r.get("mean_anom_deg", "")) if has_ma else None
        ecc     = _to_float_or_none(r.get("eccentricity", "")) if has_ecc else None
        src     = _to_str_or_none(r.get("phase_src", "")) if has_src else None
        P_h     = _to_float_or_none(r.get("period_hours", "")) if has_P else None
        det     = _to_str_or_none(r.get("detection", "")) if has_det else None

        # Normalize detection flag for use in ID
        det_norm = (det or "").strip().lower()
        if det_norm == "y":
            det_tag = "Y"
        elif det_norm == "n":
            det_tag = "N"
        else:
            det_tag = "U"   # unknown / missing

        # Stable encounter_id: observer + start time + detection
        safe_time = utc_iso.replace("-", "").replace(":", "").replace("T", "_").replace("Z", "Z")
        safe_obs = re.sub(r"[^A-Za-z0-9]+", "", observer)
        enc_id = f"PLUME_{safe_obs}_{safe_time}_{det_tag}"

        # Normalize anomalies to [0,360) if present
        if ta_deg is not None:
            try:
                ta_deg = float(ta_deg) % 360.0
            except Exception:
                ta_deg = None
        if ma_deg is not None:
            try:
                ma_deg = float(ma_deg) % 360.0
            except Exception:
                ma_deg = None

        rows.append({
            "encounter_id": enc_id,
            "encounter_tag": "plume",

            # from file
            "utc_iso": utc_iso,
            "observer": observer,
            "exposure_s": exposure_s,

            # optional precomputed (if present in file)
            "mid_utc_iso": mid_iso,
            "true_anom_deg": ta_deg,
            "mean_anom_deg": ma_deg,
            "eccentricity": ecc,
            "phase_src": src,
            "period_hours": P_h,
            "detection": det_tag,

            # keep keys EncounterRegistry might expect (even if None)
            "et_tdb_sec": None,
            "lat_deg": None,
            "lon_deg": None,
            "detected": None,
        })

    return rows



def load_plume_observations_txt(txt_path: str | None = None) -> list[dict]:
    """
    Load plume observation encounters from the default text file and return
    rows suitable for EncounterRegistry.
    """
    if txt_path is None:
        txt_path = os.path.join(DATA_DIR, "plume_observations.txt")  
    if not os.path.exists(txt_path):
        raise FileNotFoundError(f"Plume observation file not found: {txt_path}")
    return build_plume_rows_from_txt(txt_path)

_SPLIT_RE = re.compile(r"[,\t]+|\s{2,}")  # commas, tabs, or 2+ spaces

def load_juice_europa_flybys_txt(path: str | None = None) -> list[dict]:
    """
    Read JUICE_Europa_flybys.txt.

    Columns (header optional):
      Observer, Encounter, Start(ISO), true anom deg

    Returns list of dict rows with normalized keys:
      observer, encounter, utc_iso, true_anom_deg
    """
    if path is None:
        path = os.path.join(DATA_DIR, "JUICE_Europa_flybys.txt")

    rows: list[dict] = []
    with open(path, "r", encoding="utf-8") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            low = line.lower()
            # header-ish line
            if ("observer" in low) and ("encounter" in low) and ("start" in low) and ("anom" in low):
                continue

            parts = [p.strip() for p in _SPLIT_RE.split(line) if p.strip()]
            if len(parts) < 4:
                raise ValueError(
                    f"{path}:{line_no}: expected 4 columns (Observer, Encounter, Start(ISO), true anom deg); "
                    f"got {len(parts)}: {parts}"
                )

            observer, encounter, start_iso, nu = parts[0], parts[1], parts[2], parts[3]

            rows.append({
                "observer": observer,
                "encounter": encounter,
                "utc_iso": start_iso,
                "true_anom_deg": nu,  # keep as str; normalized later
            })

    return rows

