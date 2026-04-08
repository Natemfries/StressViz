# stressviz/encounters.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Callable, Any
import math
import datetime as _dt

@dataclass
class Encounter:
    enc_id: str
    tag: Optional[str]
    et_tdb_sec: Optional[float]
    utc_iso: Optional[str]
    lat_deg: Optional[float]
    lon_deg: Optional[float]
    nu_deg: Optional[float]
    M_deg: Optional[float]
    period_hours: Optional[float]

class EncounterRegistry:
    def __init__(
        self,
        get_eccentricity: Optional[Callable[[], float]] = None,
        default_period_hours: float = 85.228,
        utc_default: Optional[_dt.datetime] = _dt.datetime(2000, 1, 1, tzinfo=_dt.timezone.utc),
        true2mean: Optional[Callable[[float, float], Optional[float]]] = None,
    ):
        self._get_e = get_eccentricity or (lambda: 0.0)
        self._P_h   = float(default_period_hours)
        self._utc_default = utc_default
        self._true2mean = true2mean or (lambda nu, e: None)
        self.by_id: Dict[str, Encounter] = {}

    def build_registry_from_raw(self, rows: Iterable[dict]) -> None:
        """
        rows: iterable of dicts with keys like
          encounter_id, encounter_tag, et_tdb_sec, utc_iso, true_anom_deg, mean_anom_deg,
          lat_deg, lon_deg, period_hours
        """
        e = 0.0
        try:
            e = float(self._get_e())
        except Exception:
            pass

        out: Dict[str, Encounter] = dict(self.by_id)
        for r in rows:
            enc_id = str(r.get("encounter_id") or r.get("enc_id") or "").strip()
            if not enc_id:
                continue

            et  = _f(r.get("et_tdb_sec"))
            nu  = _f(r.get("true_anom_deg") or r.get("nu_deg"))
            M   = _f(r.get("mean_anom_deg"))
            if M is None and nu is not None:
                try:
                    mm = self._true2mean(float(nu), float(e))
                    M = None if mm is None else float(mm) % 360.0
                except Exception:
                    M = None

            P_h = _f(r.get("period_hours"))
            if not (isinstance(P_h, float) and P_h > 0):
                P_h = self._P_h

            out[enc_id] = Encounter(
                enc_id=enc_id,
                tag=(r.get("encounter_tag") or None),
                et_tdb_sec=et,
                utc_iso=(r.get("utc_iso") or None),
                lat_deg=_f(r.get("lat_deg") or r.get("lat")),
                lon_deg=_f(r.get("lon_deg") or r.get("lon")),
                nu_deg=nu,
                M_deg=M if M is not None else 0.0,   # never None; safe default
                period_hours=P_h,
            )
        self.by_id = out

    def list_ids(self):
        return sorted(self.by_id.keys(), key=str)

    def get(self, enc_id: str) -> Optional[dict]:
        rec = self.by_id.get(enc_id)
        if not rec:
            return None
        return {
            "encounter_id": rec.enc_id,
            "encounter_tag": rec.tag,
            "et_tdb_sec": rec.et_tdb_sec,
            "utc_iso": rec.utc_iso,
            "lat_deg": rec.lat_deg,
            "lon_deg": rec.lon_deg,
            "true_anom_deg": rec.nu_deg,
            "mean_anom_deg": rec.M_deg,
            "period_hours": rec.period_hours,
        }

    def to_scalar_meta(self) -> Dict[str, dict]:
        """
        Format expected by ScalarPlotPanel.register_encounters:
        {enc_id: {'utc_ca': dt, 'M_ca_deg': float, 'period_hours': float}}
        """
        meta = {}
        for enc_id, rec in self.by_id.items():
            utc_ca = _parse_iso(rec.utc_iso) or self._utc_default
            meta[enc_id] = {
                "utc_ca": utc_ca,
                "M_ca_deg": float(rec.M_deg if rec.M_deg is not None else 0.0),
                "period_hours": float(rec.period_hours if rec.period_hours else self._P_h),
            }
        return meta
    
    def _ensure_plume_state_for_selected(self, plume_id: str):
        """
        Resolve ET/true/M for ONE selected plume observation.
        Must cache results somewhere that later plotting code actually reads.
        """
        meta = self._enc_meta_by_id.get(plume_id, {})
        if meta.get("mean_anom_deg") is not None:
            return  # already cached

        utc = meta.get("utc_iso") or meta.get("utc_ca")
        if not utc:
            return
        
        et_sec, true_anom_deg, period_hours = self._query_horizons_for_plume_utc(utc)

        meta["et_tdb_sec"] = et_sec
        meta["true_anom_deg"] = true_anom_deg
        meta["period_hours"] = period_hours

        if true_anom_deg is not None:
            meta["mean_anom_deg"] = float(_true2mean(true_anom_deg)) % 360.0

        self._enc_meta_by_id[plume_id] = meta

        if meta.get("mean_anom_deg") is not None:
            self._selected_encounter_M_deg = meta["mean_anom_deg"]



def _parse_iso(s: Optional[str]):
    if not s:
        return None
    s = s.strip()
    if s.endswith("Z"):
        s = s[:-1]
    try:
        return _dt.datetime.fromisoformat(s).replace(tzinfo=_dt.timezone.utc)
    except Exception:
        return None

def _f(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        sx = str(x).strip()
        if not sx or sx.lower() == "nan":
            return None
        return float(sx)
    except Exception:
        return None

