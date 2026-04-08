# stressviz/sat_file_parser.py
import re
from dataclasses import dataclass, field
from typing import List, Optional, Dict

# ---------- Canonical schema (SatStress GUI compatible) ----------
@dataclass
class Layer:
    name: str
    density_kgm3: Optional[float] = None          # ρ
    young_pa: Optional[float] = None              # E
    poisson: Optional[float] = None               # ν
    thickness_m: Optional[float] = None           # H (meters)
    viscosity_pas: Optional[float] = None         # η

@dataclass
class Structure:
    system_id: str = "JupiterEuropa"
    planet_mass_kg: Optional[float] = None
    orbit_eccentricity: Optional[float] = None
    semi_major_axis_m: Optional[float] = None
    nsr_period_yrs: Optional[float] = None
    layers: List[Layer] = field(default_factory=list)

class SatFileParseError(Exception):
    pass

# ---------- Helpers ----------
ALIASES = {"ICE1": "ICE_UPPER", "ICE2": "ICE_LOWER", "ROCK": "CORE"}
ORDER = ["ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE"]

def _norm_name(s: str) -> str:
    return ALIASES.get(s.strip().upper(), s.strip().upper())

def _rank(name: str) -> int:
    try: return ORDER.index(_norm_name(name))
    except ValueError: return 10_000

def _is_number(s: str) -> bool:
    try:
        float(str(s).replace(",", ""))
        return True
    except Exception:
        return False

def _f(s) -> Optional[float]:
    try:
        return float(str(s).replace(",", ""))
    except Exception:
        return None

# ---------- Parser ----------
def parse_sat_file(path: str) -> Structure:
    """
    Parse a SatStress-style .sat file (as displayed by SatStress GUI):
      Header keys (case-insensitive): System ID, Planet Mass, Orbit Eccentricity,
      Orbit Semimajor Axis [m], NSR Period [yrs].
      Layer rows: LayerID, Density [kg/m^3], Young's [Pa], Poisson's ratio,
                  Thickness [m], Viscosity [Pa s]
    Comments starting with # or ! are ignored. Blank lines are skipped.
    Column header line for the layer table is optional.
    """
    try:
        with open(path, "r") as f:
            raw = f.readlines()
    except OSError as e:
        raise SatFileParseError(f"Could not read {path}: {e}")

    # strip comments and blanks
    lines: List[str] = []
    for ln in raw:
        ln = re.split(r"[#!]", ln, maxsplit=1)[0].strip()
        if ln:
            lines.append(ln)

    st = Structure()

    # --- headers ---
    # accept "key = val", "key: val", or "key  val"
    header_re = re.compile(r"^\s*([A-Za-z][A-Za-z0-9 _\[\]/']*?)\s*(?:[:=]\s*|\s{2,}|\s)(.+?)\s*$")
    for ln in lines:
        m = header_re.match(ln)
        if not m:
            continue
        k_raw, v_raw = m.group(1).strip(), m.group(2).strip()
        k = k_raw.lower().replace("'", "").replace("  ", " ")
        # map to canonical fields
        if "system" in k:
            st.system_id = v_raw
        elif "planet" in k and "mass" in k:
            st.planet_mass_kg = _f(v_raw)
        elif "eccentr" in k:
            st.orbit_eccentricity = _f(v_raw)
        elif ("semi" in k and "axis" in k) or ("a" == k.strip()):
            st.semi_major_axis_m = _f(v_raw)
        elif "nsr" in k and ("period" in k or "yrs" in k):
            st.nsr_period_yrs = _f(v_raw)

    # --- locate (optional) layer header row to map columns by name ---
    col_names: Optional[List[str]] = None
    for ln in lines:
        if re.search(r"density|young|poisson|thick|viscos", ln, re.I):
            # tokenize by commas or whitespace
            cand = [t.strip().lower() for t in re.split(r"[,\s]+", ln) if t.strip()]
            # must contain at least 3 known column-ish words
            score = sum(int(any(w in c for w in ("density","young","poisson","thick","viscos","layer")))
                        for c in cand)
            if score >= 3:
                col_names = cand
                break

    # --- layer rows ---
    # lines that look like: [NAME] v v v v v
    layers: List[Layer] = []
    for ln in lines:
        toks = re.split(r"[,\s]+", ln.strip())
        if not toks:
            continue
        has_name = not _is_number(toks[0])
        nums_str = toks[1:] if has_name else toks
        if len(nums_str) < 3 or not all(_is_number(x) for x in nums_str[:3]):
            continue  # not a numeric row
        name = _norm_name(toks[0] if has_name else f"Layer{len(layers)+1}")
        nums = [_f(x) for x in nums_str]

        if col_names:
            # map by names (best-effort)
            m: Dict[str, Optional[float]] = {}
            for i, n in enumerate(nums):
                if i >= len(col_names):
                    break
                m[col_names[i]] = n
            rho = m.get("density") or m.get("density[kg/m3]") or m.get("density_kg/m3")
            E   = (m.get("young") or m.get("youngs") or m.get("young's") or
                   m.get("youngsmodulus") or m.get("young'smodulus") or m.get("young's") or
                   m.get("young's") or m.get("young's") )
            nu  = (m.get("poisson") or m.get("poisson's") or m.get("poissons") or
                   m.get("poisson_ratio") or m.get("poisson'sratio"))
            Hm  = m.get("thickness") or m.get("thickness[m]") or m.get("thickness_m")
            eta = m.get("viscosity") or m.get("viscosity[pa s]") or m.get("viscosity[pa*s]") or m.get("viscosity_pas")
            # accept km too
            Hkm = m.get("thickness[km]") or m.get("thickness_km")
            if Hm is None and Hkm is not None:
                Hm = 1000.0 * Hkm
        else:
            # default SatStress order by position:
            # Density, Young, Poisson, Thickness(m), Viscosity
            rho = nums[0] if len(nums) > 0 else None
            E   = nums[1] if len(nums) > 1 else None
            nu  = nums[2] if len(nums) > 2 else None
            Hm  = nums[3] if len(nums) > 3 else None
            eta = nums[4] if len(nums) > 4 else None

        # finalize layer
        L = Layer(
            name=name,
            density_kgm3=rho,
            young_pa=E,
            poisson=nu,
            thickness_m=Hm,
            viscosity_pas=eta,
        )
        layers.append(L)

    if not layers:
        raise SatFileParseError("No layer rows detected in .sat.")

    # sort to SatStress GUI order (outer→inner)
    layers.sort(key=lambda L: _rank(L.name))
    st.layers = layers
    return st

