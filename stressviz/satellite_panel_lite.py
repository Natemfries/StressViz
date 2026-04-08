# stressviz/satellite_panel_lite.py
import wx

from .satstress_interface import load_satellite_from_file, SatFileLoadError

# ----------------------------- canonical order & utils -----------------------------

ORDER = ["ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE"]
ALIASES = {"ICE1": "ICE_UPPER", "ICE2": "ICE_LOWER", "ROCK": "CORE"}


def _norm_id(s):
    if s is None:
        return ""
    name = str(s).strip().upper()
    return ALIASES.get(name, name)


def _rank(name: str) -> int:
    name = _norm_id(name)
    try:
        return ORDER.index(name)
    except ValueError:
        return 10_000  # unknowns at the end (stable order among themselves)


def _coerce_float(v):
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        try:
            s = str(v).replace(",", "").strip()
            if not s:
                return None
            if s.lower() in ("inf", "infinity"):
                return None
            return float(s)
        except Exception:
            return None


def _fmt_sci(v):
    """SatStress-like scientific formatting (e.g., 3.30E-01, 9.28E+09)."""
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.3E}"
    except Exception:
        return str(v)


def _get_field(obj, cand_keys):
    """Get a field from dict-like or object-like using any candidate key."""
    if isinstance(obj, dict):
        for k in cand_keys:
            if k in obj:
                return obj[k]
        low = {str(k).lower(): v for k, v in obj.items()}
        for k in cand_keys:
            lk = str(k).lower()
            if lk in low:
                return low[lk]
        return None
    for k in cand_keys:
        if hasattr(obj, k):
            return getattr(obj, k)
    names = dir(obj)
    lowmap = {n.lower(): n for n in names}
    for k in cand_keys:
        lk = str(k).lower()
        if lk in lowmap:
            return getattr(obj, lowmap[lk])
    return None


def _layer_to_satstress_dict(L):
    """
    Normalize any SatStress layer object/dict to SatStress-like DTO keys:
      { density, young, poisson, thickness_m, viscosity }
    Thickness coerced to meters if provided in km.
    """
    density = _get_field(L, ["DENSITY", "density", "rho", "RHO"])
    young = _get_field(L, ["YOUNGS_MODULUS", "young", "youngs", "E", "YOUNG"])
    poisson = _get_field(L, ["POISSONS_RATIO", "poisson", "nu", "POISSON", "NU"])
    thick_m = _get_field(L, ["THICKNESS", "thickness", "THICKNESS_M", "H_m"])
    if thick_m is None:
        thick_km = _get_field(L, ["THICKNESS_KM", "thickness_km", "H_km"])
        if thick_km is not None:
            thick_m = _coerce_float(thick_km) * 1000.0
    viscosity = _get_field(L, ["VISCOSITY", "viscosity", "eta", "ETA"])

    return {
        "density": _coerce_float(density),
        "young": _coerce_float(young),
        "poisson": _coerce_float(poisson),
        "thickness_m": _coerce_float(thick_m),
        "viscosity": _coerce_float(viscosity),
    }


def _sat_to_dto(sat):
    """
    Normalize a loaded SatStress satellite (or our DTO) to:
    {
      system_id, planet_mass_kg, ecc, a_m, nsr_years,
      layers: [ {density, young, poisson, thickness_m, viscosity} ... ]  (sorted outer→inner)
    }
    """
    if sat is None:
        return None

    if isinstance(sat, dict):
        system_id = sat.get("system_id") or sat.get("SYSTEM_ID")
        planet_mass_kg = sat.get("planet_mass_kg") or sat.get("PLANET_MASS")
        ecc = sat.get("ecc") or sat.get("ORBIT_ECCENTRICITY")
        a_m = (
            sat.get("a_m")
            or sat.get("ORBIT_SEMIMAJOR_AXIS")
            or sat.get("SEMI_MAJOR_AXIS_M")
        )
        if a_m is None:
            a_km = sat.get("SEMI_MAJOR_AXIS_KM")
            if a_km is not None:
                a_m = _coerce_float(a_km) * 1000.0
        nsr_years = sat.get("nsr_years") or sat.get("NSR_PERIOD")
        layers_in = sat.get("layers", [])
    else:
        system_id = _get_field(sat, ["SYSTEM_ID", "system_id", "SystemID", "ID", "name", "Name"])
        planet_mass_kg = _get_field(
            sat,
            [
                "PLANET_MASS",
                "planet_mass",
                "planet_mass_kg",
                "planetMass",
                "Mplanet",
                "M_p",
                "mass",
                "primary_mass",
            ],
        )
        ecc = _get_field(sat, ["ORBIT_ECCENTRICITY", "ecc", "Ecc", "e", "eccentricity", "Eccentricity"])
        a_m = _get_field(
            sat,
            [
                "ORBIT_SEMIMAJOR_AXIS",
                "a_m",
                "a",
                "semi_major_axis_m",
                "semimajor_axis_m",
                "semiMajorAxis_m",
            ],
        )
        if a_m is None:
            a_km = _get_field(sat, ["SEMI_MAJOR_AXIS_KM", "semimajor_axis_km", "a_km"])
            if a_km is not None:
                a_m = _coerce_float(a_km) * 1000.0
        nsr_years = _get_field(sat, ["NSR_PERIOD", "nsr_years", "NSR", "nsr"])
        layers_in = _get_field(sat, ["layers", "Layers", "satlayers", "sat_layers"]) or []

    tmp = []
    for item in layers_in:
        name = _norm_id(getattr(item, "name", getattr(item, "ID", getattr(item, "layer_id", ""))))
        tmp.append((name, _layer_to_satstress_dict(item)))

    if tmp:
        tmp.sort(key=lambda t: _rank(t[0]))

    norm_layers = [d for _, d in tmp]

    return {
        "system_id": system_id if (system_id is None or isinstance(system_id, str)) else str(system_id),
        "planet_mass_kg": _coerce_float(planet_mass_kg),
        "ecc": _coerce_float(ecc),
        "a_m": _coerce_float(a_m),
        "nsr_years": nsr_years,  # may be "infinity" string
        "layers": norm_layers,
    }

EUROPA_DEFAULT_DTO = {
    "system_id": "JupiterEuropa",
    "planet_mass_kg": 1.8987e27,
    "ecc": 0.0094,
    "a_m": 6.709e8,
    # NSR_PERIOD from the file is 3.1556926e12 s → ~9.999786422288134e4 years
    "nsr_years": 99997.86422288134,

    "layers": [
        {  # ICE_UPPER
            "layer_id": "ICE_UPPER",
            "density": 917.0,
            "young": 9.28e9,
            "poisson": 0.33,
            "thickness_m": 12000.0,
            "viscosity": 1.0e22,
        },
        {  # ICE_LOWER
            "layer_id": "ICE_LOWER",
            "density": 917.0,
            "young": 9.28e9,
            "poisson": 0.331,
            "thickness_m": 8000.0,
            "viscosity": 1.0e17,
        },
        {  # OCEAN
            "layer_id": "OCEAN",
            "density": 1000.0,
            "young": 0.0,
            "poisson": 0.5,
            "thickness_m": 100000.0,
            "viscosity": 0.0,
        },
        {  # CORE
            "layer_id": "CORE",
            "density": 3847.6,
            "young": 1.0e11,
            "poisson": 0.25,
            "thickness_m": 1.391e6,
            "viscosity": 0.0,
        },
    ],
}




# --------------------------- the panel -------------------------------


class SatellitePanelLite(wx.Panel):
    """
    SatStress-style satellite tab (no dependency on satstressgui):
      - Load from file button (uses your satstress_interface)
      - Save to file button (present but disabled)
      - Global params: System ID, Planet Mass [kg], Orbit Eccentricity,
                       Orbit Semimajor Axis [m], NSR Period [yrs]
      - Fixed 4-layer grid of individual text fields in canonical order:
        ICE_UPPER, ICE_LOWER, OCEAN, CORE
      - Emits on_satellite_changed(obj_or_dto) on load/apply
    """

    def __init__(self, parent, on_satellite_changed=None):
        super().__init__(parent)
        print("[SatellitePanelLite] __init__ starting")
        self.on_satellite_changed = on_satellite_changed
        self.satellite = None  # object (from file) or dict DTO (manual)
        self.last_loaded_path = None

        root = wx.BoxSizer(wx.VERTICAL)

        # ---------- Preset row (above Load/Save) ----------
        preset_row = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_europa = wx.Button(self, label=u"Europa Preset")
        preset_row.Add(self.btn_europa, 0, wx.ALL, self.FromDIP(3))
        root.Add(preset_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.FromDIP(3))


        # ---------- Top row: Load / Save ----------
        top_row = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_load_file = wx.Button(self, label=u"Load from file")
        self.btn_save_file = wx.Button(self, label=u"Save to file")
        self.btn_save_file.Enable(False)  # present but disabled for now
        top_row.Add(self.btn_load_file, 1, wx.ALL | wx.EXPAND, 3)
        top_row.Add(self.btn_save_file, 1, wx.ALL | wx.EXPAND, 3)
        root.Add(top_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 3)

        # ---------- Two columns under the buttons ----------
        row = wx.BoxSizer(wx.HORIZONTAL)

        # spacing + paddings (DIP-safe)
        GUTTER    = self.FromDIP(6)   # space between left/right
        BOX_PAD   = self.FromDIP(6)    # outer padding inside the static box
        INNER_PAD = self.FromDIP(2)    # tiny pad to avoid text touching border

        # make all panels share the same background (kills that “white ring” look)
        bg = self.GetBackgroundColour()

        # ===== LEFT: compact form inside a rounded StaticBox =====
        INPUT_W = self.FromDIP(120)
        LEFT_W  = self.FromDIP(300)

        left_box   = wx.StaticBox(self, label="")
        left_sbs   = wx.StaticBoxSizer(left_box, wx.VERTICAL)
        left_inner = wx.Panel(left_box)
        left_box.SetBackgroundColour(bg)
        left_inner.SetBackgroundColour(bg)

        form = wx.FlexGridSizer(0, 2, self.FromDIP(6), self.FromDIP(8))  # label | field

        def L(lbl): 
            st = wx.StaticText(left_inner, label=lbl)
            return st

        self.txt_system_id      = wx.TextCtrl(left_inner, style=wx.TE_PROCESS_ENTER)
        self.txt_planet_mass_kg = wx.TextCtrl(left_inner, style=wx.TE_RIGHT)
        self.txt_ecc            = wx.TextCtrl(left_inner, style=wx.TE_RIGHT)
        self.txt_a_m            = wx.TextCtrl(left_inner, style=wx.TE_RIGHT)
        self.txt_nsr_years      = wx.TextCtrl(left_inner)

        for ctrl in (self.txt_system_id, self.txt_planet_mass_kg,
                     self.txt_ecc, self.txt_a_m, self.txt_nsr_years):
            ctrl.SetMinSize((INPUT_W, -1))


        def add_row(label, ctrl):
            form.Add(L(label + ":"), 0, wx.ALIGN_CENTER_VERTICAL)
            form.Add(ctrl,           0, wx.ALIGN_CENTER_VERTICAL)

        add_row("System ID",                self.txt_system_id)
        add_row("Planet Mass [kg]",         self.txt_planet_mass_kg)
        add_row("Orbit Eccentricity",       self.txt_ecc)
        add_row("Orbit Semimajor Axis [m]", self.txt_a_m)
        add_row("NSR Period [yrs]",         self.txt_nsr_years)

        left_inner.SetSizer(form)
        left_sbs.Add(left_inner, 0, wx.EXPAND | wx.ALL, INNER_PAD)
        left_sbs.SetMinSize((LEFT_W, -1))
        row.Add(left_sbs, 0, wx.ALL, self.FromDIP(6))

        row.AddSpacer(GUTTER)

        # ===== RIGHT: text-field matrix inside a rounded StaticBox =====
        props  = [("rho","Density [kg/m3]"), ("E","Young's Modulus [Pa]"),
                  ("nu","Poisson's Ratio"), ("H","Thickness [m]"), ("eta","Viscosity [Pa s]")]
        layers = tuple(globals().get("ORDER", ("ICE_UPPER","ICE_LOWER","OCEAN","CORE")))

        right_box   = wx.StaticBox(self, label="")
        right_sbs   = wx.StaticBoxSizer(right_box, wx.VERTICAL)
        right_inner = wx.Panel(right_box)
        right_box.SetBackgroundColour(bg)
        right_inner.SetBackgroundColour(bg)

        # 1 header row + N layer rows; 1 label col + len(props) data cols
        rg = wx.FlexGridSizer(rows=1+len(layers), cols=1+len(props),
                              vgap=self.FromDIP(4), hgap=self.FromDIP(6))
        # equal width for all data columns
        for c in range(1, 1+len(props)):
            rg.AddGrowableCol(c, 1)

        def _header(text):
            st = wx.StaticText(right_inner, label=text)
            f = st.GetFont(); f.MakeBold(); st.SetFont(f)
            return st

        # header row
        rg.Add(_header("Layer ID"), 0, wx.ALIGN_CENTER_VERTICAL)
        for _, hdr in props:
            rg.Add(_header(hdr), 0, wx.ALIGN_CENTER_VERTICAL)

        # data rows
        self.layer_fields = {}
        CELL_MIN_W = self.FromDIP(120)

        for layer in layers:
            rg.Add(wx.StaticText(right_inner, label=layer), 0, wx.ALIGN_CENTER_VERTICAL)
            for key, _hdr in props:
                tc = wx.TextCtrl(right_inner, style=wx.TE_RIGHT)
                tc.SetMinSize((CELL_MIN_W, -1))
                self.layer_fields[(layer, key)] = tc
                rg.Add(tc, 1, wx.EXPAND)

        right_inner.SetSizer(rg)
        right_sbs.Add(right_inner, 1, wx.EXPAND | wx.ALL, INNER_PAD)
        row.Add(right_sbs, 1, wx.ALL | wx.EXPAND, self.FromDIP(6))

        # attach to root
        root.Add(row, 1, wx.EXPAND)

        self.SetSizer(root)
        self.Layout()

        # ---------- Events ----------
        self.btn_load_file.Bind(wx.EVT_BUTTON, self._on_click_load_file)
        self.btn_europa.Bind(wx.EVT_BUTTON, self._on_click_europa_defaults)

        # self.btn_save_file.Bind(wx.EVT_BUTTON, self._on_click_save_file)  # when implemented

    # ---------------------- event handlers ----------------------

    def _on_click_load_file(self, evt):
        with wx.FileDialog(
            self,
            message="Load from satellite file",
            wildcard="Satellite files (*.sat;*.satellite;*.sat.bak)|*.sat;*.satellite;*.sat.bak|All files|*.*",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
        ) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            path = dlg.GetPath()
            self.last_loaded_path = path

        try:
            sat = load_satellite_from_file(path)
            self.satellite = sat
            self._mirror_satellite_into_form(sat)
            self._notify()
        except SatFileLoadError as e:
            wx.MessageBox(str(e), "SatStress import error", wx.OK | wx.ICON_ERROR)
        except Exception as e:
            wx.MessageBox(f"Failed to load .sat:\n{e}", "Parse error", wx.OK | wx.ICON_ERROR)

    def _on_apply(self, evt):
        dto = self._dto_from_form()
        self.satellite = dto
        self._notify()

    def _notify(self):
        if self.on_satellite_changed:
            self.on_satellite_changed(self.satellite)

    def get_system_id(self) -> str:
        return self.txt_system_id.GetValue().strip()

    # --------------------- UI <-> DTO helpers --------------------

    def _dto_from_form(self):
        def getf(ctrl):
            s = ctrl.GetValue().strip()
            if not s:
                return None
            if s.lower() == "infinity":
                return "infinity"
            try:
                return float(s.replace(",", ""))
            except Exception:
                return s  # leave as-is

        dto = {
            "system_id": self.txt_system_id.GetValue().strip() or None,
            "planet_mass_kg": getf(self.txt_planet_mass_kg),
            "ecc": getf(self.txt_ecc),
            "a_m": getf(self.txt_a_m),
            "nsr_years": self.txt_nsr_years.GetValue().strip() or None,
            "layers": [],
        }

        # Build layers from individual text boxes in canonical order
        def val(layer, key):
            ctrl = self.layer_fields.get((layer, key))
            if not ctrl:
                return None
            s = ctrl.GetValue().strip()
            if not s:
                return None
            try:
                return float(s.replace(",", ""))
            except Exception:
                return s

        for lname in ORDER:
            dto["layers"].append(
                {
                    "layer_id": lname,
                    "density": val(lname, "rho"),
                    "young": val(lname, "E"),
                    "poisson": val(lname, "nu"),
                    "thickness_m": val(lname, "H"),
                    "viscosity": val(lname, "eta"),
                }
            )
        return dto
    
    def _apply_dto_to_ui(self, dto: dict):
        """Push a normalized DTO into the current controls (left globals + right layer fields)."""

        def setv(ctrl, val, sci=False, allow_inf=False):
            if val is None or val == "":
                s = ""
            elif allow_inf and isinstance(val, str) and val.strip().lower() == "infinity":
                s = "infinity"
            else:
                try:
                    f = float(val)
                    s = _fmt_sci(f) if sci else str(val)
                except Exception:
                    s = str(val)
            ctrl.SetValue(s)

        # Left side (system_id as-is; numbers in sci)
        setv(self.txt_system_id,      dto.get("system_id"), sci=False)
        setv(self.txt_planet_mass_kg, dto.get("planet_mass_kg"), sci=True)
        setv(self.txt_ecc,            dto.get("ecc"),            sci=True)
        setv(self.txt_a_m,            dto.get("a_m"),            sci=True)
        setv(self.txt_nsr_years,      dto.get("nsr_years"),      sci=True, allow_inf=True)

        # Right side (per-layer) — everything numeric in sci
        want = ("ICE_UPPER","ICE_LOWER","OCEAN","CORE")
        by_name = {}
        for L in (dto.get("layers") or []):
            nm = str(L.get("layer_id") or L.get("name") or "").strip().upper()
            if nm:
                by_name[nm] = L

        for i, lname in enumerate(want):
            row = by_name.get(lname)
            if row is None:
                layers = dto.get("layers") or []
                row = layers[i] if i < len(layers) else {}

            mapping = {
                "rho": row.get("density"),
                "E":   row.get("young"),
                "nu":  row.get("poisson"),
                "H":   row.get("thickness_m"),
                "eta": row.get("viscosity"),
            }
            for key, val in mapping.items():
                ctrl = self.layer_fields.get((lname, key))
                if ctrl:
                    try:
                        # format as sci unless blank
                        s = "" if val in (None, "") else _fmt_sci(float(val))
                    except Exception:
                        s = "" if val in (None, "") else str(val)
                    ctrl.SetValue(s)


    def _on_click_europa_defaults(self, _evt):
        """Fill UI with the Europa preset and notify parent."""
        dto = {
            "system_id": EUROPA_DEFAULT_DTO["system_id"],
            "planet_mass_kg": EUROPA_DEFAULT_DTO["planet_mass_kg"],
            "ecc": EUROPA_DEFAULT_DTO["ecc"],
            "a_m": EUROPA_DEFAULT_DTO["a_m"],
            "nsr_years": EUROPA_DEFAULT_DTO["nsr_years"],
            "layers": [dict(L) for L in EUROPA_DEFAULT_DTO.get("layers", [])],
        }
        self._apply_dto_to_ui(dto)
        self.satellite = dto
        self._notify()


    def _mirror_satellite_into_form(self, sat):
        """
        Fill the left globals and the right per-layer fields (self.layer_fields).
        Designed for dict files with keys:

          SYSTEM_ID, PLANET_MASS, ORBIT_ECCENTRICITY, ORBIT_SEMIMAJOR_AXIS, NSR_PERIOD

        Layers via indexed keys (0/1/2/3):
          LAYER_ID_i, DENSITY_i, YOUNGS_MODULUS_i, POISSONS_RATIO_i, THICKNESS_i, VISCOSITY_i

        Also supports object-style sat["layers"] if present.
        """
        import re

        PROPS = ("rho", "E", "nu", "H", "eta")
        SEC_PER_YEAR = 365.25 * 86400.0

        def _fmt(v):
            if v in (None, ""):
                return ""
            try:
                f = float(v)
            except Exception:
                return str(v)
            if f == 0:
                return "0"
            return f"{f:.6g}" if (abs(f) >= 1e6 or abs(f) < 1e-3) else f"{f:g}"

        def _num(v):
            if v is None:
                return None
            if isinstance(v, (int, float)):
                return float(v)
            s = str(v).replace(",", " ").strip()
            m = re.search(r"([-+]?\d+(\.\d+)?([eE][-+]?\d+)?)", s)
            return float(m.group(1)) if m else None

        def _get(obj, *names, default=None):
            for n in names:
                if hasattr(obj, n):
                    return getattr(obj, n)
                if isinstance(obj, dict):
                    if n in obj:
                        return obj[n]
                    ln = str(n).lower()
                    for k in obj.keys():
                        if str(k).lower() == ln:
                            return obj[k]
            return default

        # ---------- LEFT: globals ----------
        sys_id = (sat.get("SYSTEM_ID") if isinstance(sat, dict) else None) or _get(
            sat, "name", "ID", "system_id", "system", "body_name"
        )
        self.txt_system_id.SetValue(_fmt(sys_id))

        planet_mass = (sat.get("PLANET_MASS") if isinstance(sat, dict) else None) or _get(
            sat, "planet_mass", "planetMass", "M_primary", "primary_mass", "mass_primary", "M_p", "Mp"
        )
        self.txt_planet_mass_kg.SetValue(_fmt(planet_mass))

        ecc = (sat.get("ORBIT_ECCENTRICITY") if isinstance(sat, dict) else None) or _get(
            sat, "ecc", "eccentricity", "orbit_eccentricity", "Eccentricity", "e"
        )
        self.txt_ecc.SetValue(_fmt(ecc))

        a_val = None
        if isinstance(sat, dict):
            if "ORBIT_SEMIMAJOR_AXIS" in sat:
                a_val = _num(sat["ORBIT_SEMIMAJOR_AXIS"])
            if a_val is None:
                for k in (
                    "Orbit Semimajor Axis [m]",
                    "Semi-major axis [m]",
                    "Semimajor axis [m]",
                    "Semimajor Axis (m)",
                    "a [m]",
                ):
                    if k in sat:
                        a_val = _num(sat[k])
                        break
            if a_val is None:
                for k, v in sat.items():
                    kl = str(k).lower()
                    if "semi" in kl and "major" in kl and "axis" in kl:
                        a_val = _num(v)
                        if a_val is not None and ("km" in kl or " km" in str(v).lower()):
                            a_val *= 1000.0
                        break
        if a_val is None:
            a_val = _num(
                _get(
                    sat,
                    "a",
                    "a_m",
                    "semi_major_axis",
                    "semimajor_axis",
                    "semiMajorAxis",
                    "a_SI",
                    "a_meters",
                )
            )
        self.txt_a_m.SetValue(_fmt(a_val))

        nsr_years = None
        if isinstance(sat, dict) and "NSR_PERIOD" in sat:
            sec = _num(sat["NSR_PERIOD"])
            if sec and sec > 0:
                nsr_years = sec / SEC_PER_YEAR
        if nsr_years is None:
            nsr_years = _get(sat, "nsr_years", "nsrPeriodYears", "nsr_period_years")
            if nsr_years is None:
                sec = _num(_get(sat, "nsr_seconds", "nsrPeriod", "nsr_period", "NSR"))
                if sec and sec > 0:
                    nsr_years = sec / SEC_PER_YEAR
        self.txt_nsr_years.SetValue(_fmt(nsr_years))

        # ---------- RIGHT: canonical layer map ----------
        layer_map = {ln: {} for ln in ORDER}

        # (A) object-style list (optional)
        layers_obj = _get(sat, "layers", "layer_list", "Layers", default=[]) or []
        for L in layers_obj:
            nm = str(_get(L, "name", "ID", "label", "layer_name", default="")).strip().upper()
            nm = ALIASES.get(nm, nm)
            if not nm or nm not in layer_map:
                continue
            layer_map[nm].update(
                {
                    "rho": _get(L, "rho", "density", "Density"),
                    "E": _get(L, "E", "young", "young_modulus"),
                    "nu": _get(L, "nu", "poisson", "poisson_ratio"),
                    "H": _get(L, "H", "thickness", "Thickness"),
                    "eta": _get(L, "eta", "viscosity", "Viscosity"),
                }
            )

        # (B) indexed keys (allow 0-based): ..._0, ..._1, ...
        if isinstance(sat, dict):
            pat = re.compile(
                r"^(LAYER_ID|DENSITY|YOUNGS_MODULUS|POISSONS_RATIO|THICKNESS|VISCOSITY)[_\s]*(\d+)$", re.I
            )
            idx_rec = {}
            for k, v in sat.items():
                m = pat.match(str(k).strip())
                if not m:
                    continue
                base = m.group(1).upper()
                i = int(m.group(2))
                rec = idx_rec.setdefault(i, {})
                rec[base] = v

            by_index = []
            for i in sorted(idx_rec.keys()):
                rec = idx_rec[i]
                raw_name = rec.get("LAYER_ID")
                name = (str(raw_name).strip().upper() if raw_name else "")
                name = ALIASES.get(name, name)

                row = {
                    "rho": rec.get("DENSITY"),
                    "E": rec.get("YOUNGS_MODULUS"),
                    "nu": rec.get("POISSONS_RATIO"),
                    "H": rec.get("THICKNESS"),
                    "eta": rec.get("VISCOSITY"),
                }

                if name and name in layer_map:
                    for k2, v2 in row.items():
                        if v2 is not None:
                            layer_map[name][k2] = v2
                else:
                    by_index.append(row)

            # Fill any still-empty canonical rows by leftover index order
            for pos, lname in enumerate(ORDER):
                if not layer_map[lname] and pos < len(by_index):
                    layer_map[lname] = by_index[pos]

        # ---------- Push to UI ----------
        for lname in ORDER:
            rec = layer_map.get(lname, {})
            for key in ("rho", "E", "nu", "H", "eta"):
                ctrl = self.layer_fields.get((lname, key))
                if ctrl:
                    ctrl.SetValue(_fmt(rec.get(key)))

    # Optional: summary refresher if you later re-add a matrix view
    def _refresh_matrix_from_form(self):
        dto = _sat_to_dto(self.satellite) or self._dto_from_form()

        rows = [
            ("System ID", dto.get("system_id")),
            ("Planet Mass [kg]", dto.get("planet_mass_kg")),
            ("Orbit Eccentricity", dto.get("ecc")),
            ("Orbit Semimajor Axis [m]", dto.get("a_m")),
            ("NSR Period [yrs]", dto.get("nsr_years")),
            ("Layers", len(dto.get("layers", []))),
        ]
        for i, L in enumerate(dto.get("layers", [])):
            rows += [
                (f"Layer {i+1} Density [kg/m3]", _fmt_sci(L.get("density"))),
                (f"Layer {i+1} Young [Pa]", _fmt_sci(L.get("young"))),
                (f"Layer {i+1} Poisson", _fmt_sci(L.get("poisson"))),
                (f"Layer {i+1} Thickness [m]", _fmt_sci(L.get("thickness_m"))),
                (f"Layer {i+1} Viscosity [Pa s]", _fmt_sci(L.get("viscosity"))),
            ]

        g = getattr(self, "matrix_grid", None)
        if g is None:
            return
        if g.GetNumberRows() > 0:
            g.DeleteRows(0, g.GetNumberRows())
        if rows:
            g.AppendRows(len(rows))
            for r, (k, v) in enumerate(rows):
                g.SetCellValue(r, 0, str(k))
                g.SetCellValue(r, 1, "" if v is None else str(v))
            g.AutoSizeColumns()

    def get_loaded_path(self):
        return self.last_loaded_path
    
    # --- Build a DTO from the Satellite panel (robust, no dependency on grid) ---
def _get_current_dto(self):
    sp = getattr(self, "sat_panel", None)
    if sp is None:
        return None

    # If SatellitePanelLite already has get_dto(), use it.
    if hasattr(sp, "get_dto"):
        try:
            return sp.get_dto()
        except Exception:
            pass

    def _num_from(ctrl):
        if not ctrl:
            return None
        try:
            s = ctrl.GetValue().strip()
        except Exception:
            return None
        if not s or s.lower() == "infinity":
            return None
        try:
            return float(s.replace(",", ""))
        except Exception:
            return s  # leave strings as-is

    dto = {
        "system_id": (sp.txt_system_id.GetValue().strip() if hasattr(sp, "txt_system_id") else None) or "Preset",
        "planet_mass_kg": _num_from(getattr(sp, "txt_planet_mass_kg", None)),
        "ecc":            _num_from(getattr(sp, "txt_ecc", None)),
        "a_m":            _num_from(getattr(sp, "txt_a_m", None)),
        "nsr_years":      _num_from(getattr(sp, "txt_nsr_years", None)),
        "layers": [],
    }

    order  = ("ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE")
    keymap = {"rho":"density", "E":"young", "nu":"poisson", "H":"thickness_m", "eta":"viscosity"}
    lf = getattr(sp, "layer_fields", {}) or {}

    for lname in order:
        row = {"layer_id": lname}
        for k_ui, k_dto in keymap.items():
            ctrl = lf.get((lname, k_ui))
            row[k_dto] = _num_from(ctrl)
        dto["layers"].append(row)

    return dto


    # --- Materialize a temp .sat file from DTO and load it as a SatStress Satellite ---
    def _dto_to_temp_satellite(self, dto):
        import tempfile, os
        from .satstress_interface import load_satellite_from_file

        def W(lines, key, val):
            if val is None or val == "":
                return
            lines.append(f"{key} = {val}")

        lines = []
        W(lines, "SYSTEM_ID", dto.get("system_id", "Preset"))
        W(lines, "PLANET_MASS", dto.get("planet_mass_kg"))
        W(lines, "ORBIT_ECCENTRICITY", dto.get("ecc"))
        W(lines, "ORBIT_SEMIMAJOR_AXIS", dto.get("a_m"))

        ny = dto.get("nsr_years")
        if ny is not None:
            try:
                W(lines, "NSR_PERIOD", float(ny) * 365.25 * 86400.0)  # years → seconds
            except Exception:
                pass

        order = ("ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE")
        by    = { (L.get("layer_id") or "").upper(): L for L in (dto.get("layers") or []) }
        for i, lname in enumerate(order):
            L = by.get(lname, {})
            W(lines, f"LAYER_ID_{i}", lname)
            W(lines, f"DENSITY_{i}",        L.get("density"))
            W(lines, f"YOUNGS_MODULUS_{i}", L.get("young"))
            W(lines, f"POISSONS_RATIO_{i}", L.get("poisson"))
            W(lines, f"THICKNESS_{i}",      L.get("thickness_m"))
            W(lines, f"VISCOSITY_{i}",      L.get("viscosity"))

        txt = "\n".join(lines) + "\n"

        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sat", mode="w")
        try:
            tmp.write(txt); tmp.close()
            sat = load_satellite_from_file(tmp.name)
        finally:
            try: os.unlink(tmp.name)
            except Exception: pass

        return sat


    def _dto_from_current_fields(self):
        """
        Build a minimal DTO from the current SatellitePanelLite UI.
        Use this if self.satellite is a dict preset or if fields were edited.
        """
        dto = {}
        try:
            # Left globals
            dto["system_id"]      = self.sat_panel.txt_system_id.GetValue().strip() or "Preset"
            def _num(ctrl):
                s = ctrl.GetValue().strip()
                if not s or s.lower() == "infinity":
                    return None
                try:
                    return float(s.replace(",", ""))
                except Exception:
                    return None

            dto["planet_mass_kg"] = _num(self.sat_panel.txt_planet_mass_kg)
            dto["ecc"]            = _num(self.sat_panel.txt_ecc)
            dto["a_m"]            = _num(self.sat_panel.txt_a_m)

            s_nsr = self.sat_panel.txt_nsr_years.GetValue().strip()
            if s_nsr and s_nsr.lower() != "infinity":
                try:
                    dto["nsr_years"] = float(s_nsr.replace(",", ""))
                except Exception:
                    dto["nsr_years"] = None
            else:
                dto["nsr_years"] = None  # treat “infinity” / blank as not-set

            # Right layers: pull from the text boxes we created
            order = ("ICE_UPPER","ICE_LOWER","OCEAN","CORE")
            props = ("rho","E","nu","H","eta")
            dto["layers"] = []
            for lname in order:
                row = {"layer_id": lname}
                for p in props:
                    ctrl = self.sat_panel.layer_fields.get((lname, p))
                    val = None
                    if ctrl:
                        s = ctrl.GetValue().strip()
                        if s:
                            try:
                                val = float(s.replace(",", ""))
                            except Exception:
                                val = None
                    # normalize keys to what SatStress expects
                    keymap = {"rho":"density", "E":"young", "nu":"poisson", "H":"thickness_m", "eta":"viscosity"}
                    row[keymap[p]] = val
                dto["layers"].append(row)
        except Exception:
            # fallback – don’t crash the caller
            pass
        return dto






