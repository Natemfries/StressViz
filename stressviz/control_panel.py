# control_panel.py
from __future__ import annotations

import os
import sys
import math
import wx
import re
import numpy as np
import pandas as pd
import traceback
import datetime as _dt
from datetime import timezone
from typing import Optional
from bisect import bisect_left
from wx.lib.scrolledpanel import ScrolledPanel
from .encounters import EncounterRegistry
from .encounters_io import load_europa_min, load_plume_observations_txt
from .scalar_plot_panel import get_or_create_scalar_popup
from .ui_style import instr_label, style_staticbox, apply_font
from .utils import true_to_mean_anomaly_deg as _true2mean
from stressviz.onboarding import show_getting_started


# Optional: cached Love I/O (ok if missing)
try:
    from .love_io import read_love_from_sat  # , write_love_to_sat
except Exception:
    read_love_from_sat = None

# Make sure the SatStress folder (which contains the 'satstress' package) is on sys.path
SATSTRESS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "SatStress"))
if SATSTRESS_DIR not in sys.path:
    sys.path.insert(0, SATSTRESS_DIR)

# Safer import than pulling in the legacy GUI
from .satellite_panel_lite import SatellitePanelLite
from .scalar_plot_panel import ScalarPlotPanel


def _make_point_panel(parent, **kw):
    """Create the actual PointStressPanel (no fallbacks to other panels)."""
    import inspect
    from .point_panel import PointStressPanel  # <- explicit

    # Pass only kwargs the constructor accepts (avoids unexpected kw errors)
    params = inspect.signature(PointStressPanel.__init__).parameters
    allowed = {k: v for k, v in kw.items() if k in params}
    return PointStressPanel(parent, **allowed)


# ------------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------------

def _import_diurnal():
    # Import the Diurnal class from the SatStress package
    from SatStress.satstress.satstress import Diurnal
    return Diurnal

def _import_stresscalc():
    """
    Import SatStress StressCalc with a couple of safe fallbacks.
    """
    try:
        # most forks (same place as Diurnal)
        from SatStress.satstress.satstress import StressCalc
        return StressCalc
    except Exception:
        pass
    try:
        # some forks split it out
        from SatStress.satstress.stresscalc import StressCalc
        return StressCalc
    except Exception:
        pass
    try:
        # lowercase package name fallback
        from satstress.satstress import StressCalc
        return StressCalc
    except Exception:
        pass
    try:
        from satstress.stresscalc import StressCalc
        return StressCalc
    except Exception as e:
        raise ImportError("Could not import SatStress StressCalc") from e


def _try_satstress_formatter(z):
    """Try to use SatStress's own complex formatter if present."""
    try:
        # If the top-level module exposes a formatter
        from SatStress import satstressgui as _sg
        for name in ("cplxstr", "cplxStr", "format_complex", "fmt_complex"):
            if hasattr(_sg, name):
                return getattr(_sg, name)(z)
    except Exception:
        pass
    try:
        # Or a utils submodule inside the package (if it exists in your fork)
        from SatStress.satstress import utils as _su
        for name in ("cplxstr", "cplxStr", "format_complex", "fmt_complex"):
            if hasattr(_su, name):
                return getattr(_su, name)(z)
    except Exception:
        pass
    return None  # not available

def _fmt_complex_like_satstress(z, sig=6):
    """Fallback formatting similar to typical SatStress GUI style."""
    try:
        real = float(getattr(z, "real", z))
        imag = float(getattr(z, "imag", 0.0))
    except Exception:
        if isinstance(z, (tuple, list)) and len(z) >= 2:
            real, imag = float(z[0]), float(z[1])
        elif isinstance(z, dict):
            real = float(z.get("real", z.get("r", 0.0)))
            imag = float(z.get("imag", z.get("i", 0.0)))
        else:
            return str(z)
    sign = "+" if imag >= 0 else "−"
    return f"{real:.{sig}g} {sign} {abs(imag):.{sig}g}i"

def _mag_phase_tooltip(z):
    try:
        real = float(getattr(z, "real", z))
        imag = float(getattr(z, "imag", 0.0))
    except Exception:
        if isinstance(z, (tuple, list)) and len(z) >= 2:
            real, imag = float(z[0]), float(z[1])
        elif isinstance(z, dict):
            real = float(z.get("real", z.get("r", 0.0)))
            imag = float(z.get("imag", z.get("i", 0.0)))
        else:
            return None
    mag = math.hypot(real, imag)
    ph = math.degrees(math.atan2(imag, real))
    return f"|z|={mag:.6g}, phase={ph:.3f}°"

def _norm_col(name): return re.sub(r"[^a-z0-9]+","_", str(name).strip().lower()).strip("_")

def _parse_utc(s):
    try:
        dt = pd.to_datetime(s, utc=True)
        return dt.to_pydatetime() if isinstance(dt, pd.Timestamp) else None
    except Exception:
        try:
            return datetime.fromisoformat(str(s).replace("Z","+00:00")).astimezone(timezone.utc)
        except Exception:
            return None

def _fmt_iso_z(dt):
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") if isinstance(dt, datetime) else ""

ORDER = ["ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE"]
ALIASES = {"ICE1": "ICE_UPPER", "ICE2": "ICE_LOWER", "ROCK": "CORE"}

def _lname(x):
    return str(getattr(x, "name", getattr(x, "ID", ""))).strip().upper()

def _rank(name: str) -> int:
    name = ALIASES.get(name, name)
    try:
        return ORDER.index(name)
    except ValueError:
        return 10_000  # unknowns keep their relative order (stable sort)

def normalize_satellite_layer_order(sat):
    try:
        layers = list(getattr(sat, "layers", []))
        # stable sort: known names first in canonical order, others after
        sat.layers = sorted(layers, key=lambda L: _rank(_lname(L)))
    except Exception:
        pass
    return sat

# --- Horizons config for Europa around Jupiter ---
_HORIZONS_TARGET = "502"   # Europa
_HORIZONS_CENTER = "599"   # Jupiter barycenter


# ------------------------------------------------------------------------------------
# Control Panel
# ------------------------------------------------------------------------------------

class AnalysisControlPanel(ScrolledPanel):
    def __init__(self, parent, orbital_panel=None, stress_panel=None):
        super().__init__(parent, style=wx.TAB_TRAVERSAL)
        self.orbital_panel = orbital_panel
        self.stress_panel = stress_panel
        self.structure_model = None

        # SatStress Satellite object if loaded/materialized; otherwise None
        self.satellite = None

        # --- Child panels / widgets ------------------------------------------------
        self.sat_panel = SatellitePanelLite(self, on_satellite_changed=self._on_satellite_changed)

        # Love numbers area
        self.love_h_input = wx.TextCtrl(self, value="", style=wx.TE_READONLY)
        self.love_k_input = wx.TextCtrl(self, value="", style=wx.TE_READONLY)
        self.love_l_input = wx.TextCtrl(self, value="", style=wx.TE_READONLY)
        for ctrl in (self.love_h_input, self.love_k_input, self.love_l_input):
            ctrl.SetMinSize((180, -1))

        self.compute_love_btn = wx.Button(self, label="Compute Love Numbers")
        self.compute_love_btn.Enable(False)
        self.compute_love_btn.Bind(wx.EVT_BUTTON, self.on_compute_love)

        # Point panel
        self.point_panel = _make_point_panel(
            self,
            get_satellite=lambda: self.satellite,
            get_diurnal_cls=_import_diurnal,
            resolve_encounter=getattr(self, '_resolve_encounter', None),
            resolve_true_anomaly=getattr(self, '_resolve_true_anomaly', None),
            resolve_true_anomaly_from_et=getattr(self, '_resolve_true_anomaly_from_et', None),
            resolve_nu_from_time=getattr(self, '_resolve_nu_from_time', None),
            on_compute=self._on_point_compute,
        )
        self.scalar_frame = None
        self.scalar_panel = None

        get_e = None
        sp = getattr(self, "scalar_panel_ref", None)
        if sp and hasattr(sp, "_get_eccentricity"):
            get_e = sp._get_eccentricity

        self.enc_registry = EncounterRegistry(
            get_eccentricity=get_e or (lambda: 0.0),
            default_period_hours=getattr(self, "period_hours", 85.228),
            true2mean=_true2mean,
        )

        self.encounters_by_id = {}

        # ---- Grid & orbit range controls (SatStress-style) ----
        # Defaults mirror your previous hardcoded arrays (91×181 grid, 13 ν steps)
        self.txt_lat_min = wx.TextCtrl(self, value="-90")
        self.txt_lat_max = wx.TextCtrl(self, value="90")
        self.sp_lat_n    = wx.SpinCtrl(self, min=2, max=2001, initial=10)

        self.txt_lon_min = wx.TextCtrl(self, value="-180")      # use "0" if you prefer
        self.txt_lon_max = wx.TextCtrl(self, value="180")    # use "360" if you prefer
        self.sp_lon_n    = wx.SpinCtrl(self, min=2, max=4001, initial=10)

        self.txt_nu_min  = wx.TextCtrl(self, value="0")
        self.txt_nu_max  = wx.TextCtrl(self, value="360")
        self.sp_nu_n     = wx.SpinCtrl(self, min=2, max=2001, initial=10)


        self.update_btn = wx.Button(self, label="Plot")
        self.update_btn.Enable(False)
        self.Bind(wx.EVT_BUTTON, self.on_open_map, self.update_btn)

        # Getting Started button (opens the onboarding dialog)
        self.btn_getting_started = wx.Button(self, label="Getting Started")
        self.btn_getting_started.Bind(wx.EVT_BUTTON, self._on_open_getting_started)


        # --- A handy Europa preset button -----------------------------------------
        #self.btn_europa_preset = wx.Button(self, label="Europa Preset")
        #self.btn_europa_preset.Bind(wx.EVT_BUTTON, self._on_click_europa_preset)

        # --- Layout ----------------------------------------------------------------
        root = self.layout()
        if root is not None:
            self.SetSizer(root)

        # Scrolling
        self.SetupScrolling(scroll_x=False, scroll_y=True, rate_x=5, rate_y=10)

        # Try autoloading default encounters
        self._auto_load_default_encounters()

        self.Layout()
        self.FitInside()
        self.Bind(wx.EVT_SIZE, self._on_resize)

    # ============================== Layout ========================================
    def _on_resize(self, evt):
        self.Layout()
        self.FitInside()
        evt.Skip()

    def layout(self):
        sizer = wx.BoxSizer(wx.VERTICAL)

        # Header
        hdr = wx.BoxSizer(wx.HORIZONTAL)
        hdr.Add(instr_label(self, "Load or input satellite parameters:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        hdr.AddStretchSpacer(1)
        hdr.Add(self.btn_getting_started, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        sizer.Add(hdr, 0, wx.EXPAND)


        
        # Satellite panel
        sizer.Add(self.sat_panel, 0, wx.EXPAND | wx.ALL, 5)

        # ───────────────── Love (left) + Point Inputs (right) ─────────────────
        row = wx.BoxSizer(wx.HORIZONTAL)

        # Love Numbers (boxed)
        love_box   = wx.StaticBox(self, label="Input or compute diurnal love numbers:")
        style_staticbox(love_box)
        love_sbs   = wx.StaticBoxSizer(love_box, wx.VERTICAL)
        love_inner = wx.Panel(love_box)

        for w in (self.love_h_input, self.love_k_input, self.love_l_input, self.compute_love_btn):
            if w.GetParent() is not love_inner:
                w.Reparent(love_inner)

        grid = wx.FlexGridSizer(0, 2, 4, 8)
        grid.Add(wx.StaticText(love_inner, label="h2:"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.love_h_input, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(wx.StaticText(love_inner, label="k2:"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.love_k_input, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(wx.StaticText(love_inner, label="l2:"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.love_l_input, 0, wx.ALIGN_CENTER_VERTICAL)

        love_col = wx.BoxSizer(wx.VERTICAL)
        love_col.Add(grid, 0, wx.EXPAND | wx.ALL, 6)              # inner padding inside the box
        love_col.Add(self.compute_love_btn, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)
        love_inner.SetSizer(love_col)

        love_sbs.Add(love_inner, 0, wx.EXPAND, 0)

        # Add Love box with OUTER margins (no right margin between the two)
        row.Add(love_sbs, 0, wx.EXPAND | wx.LEFT | wx.TOP | wx.BOTTOM, self.FromDIP(6))

        # EXACT gap between the two boxes — tweak here
        GAP = self.FromDIP(18)
        row.AddSpacer(GAP)

        # Point Inputs mounted next to Love Numbers
        inputs_host = wx.Panel(self)
        self.point_panel.mount_inputs_into(inputs_host)

        # Add the inputs host with OUTER margins (no left margin between the two)
        row.Add(inputs_host, 1, wx.EXPAND | wx.RIGHT | wx.TOP | wx.BOTTOM, self.FromDIP(6))

        sizer.Add(row, 0, wx.EXPAND)


        # ───────────────── Results full width below ─────────────────
        results_host = wx.Panel(self)
        # requires PointStressPanel.mount_results_into
        self.point_panel.mount_results_into(results_host)
        sizer.Add(results_host, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(6))

        # ── Grid / Orbit ranges (SatStress-style) ──
        range_box   = wx.StaticBox(self, label="Grid / Orbit ranges")
        range_sbs   = wx.StaticBoxSizer(range_box, wx.VERTICAL)
        range_inner = wx.Panel(range_box)

        # Make sure the controls belong to the inner panel managed by this sizer
        for w in (
            self.txt_lat_min, self.txt_lat_max, self.sp_lat_n,
            self.txt_lon_min, self.txt_lon_max, self.sp_lon_n,
            self.txt_nu_min,  self.txt_nu_max,  self.sp_nu_n,
        ):
            if w.GetParent() is not range_inner:
                w.Reparent(range_inner)


        grid = wx.FlexGridSizer(0, 4, 6, 10)
        grid.AddGrowableCol(1, 1)
        grid.AddGrowableCol(2, 1)

        # Header
        grid.Add(wx.StaticText(range_inner, label=""),                0)
        grid.Add(wx.StaticText(range_inner, label="Minimum value"),   0)
        grid.Add(wx.StaticText(range_inner, label="Maximum value"),   0)
        grid.Add(wx.StaticText(range_inner, label="Number of grid points"), 0)

        # Latitude
        grid.Add(wx.StaticText(range_inner, label="Latitude"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.txt_lat_min, 1, wx.EXPAND)
        grid.Add(self.txt_lat_max, 1, wx.EXPAND)
        grid.Add(self.sp_lat_n,    0)

        # Longitude
        grid.Add(wx.StaticText(range_inner, label="Longitude"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.txt_lon_min, 1, wx.EXPAND)
        grid.Add(self.txt_lon_max, 1, wx.EXPAND)
        grid.Add(self.sp_lon_n,    0)

        # Second header for ν line
        grid.Add(wx.StaticText(range_inner, label=""), 0)
        grid.Add(wx.StaticText(range_inner, label="Minimum"), 0)
        grid.Add(wx.StaticText(range_inner, label="Maximum"), 0)
        grid.Add(wx.StaticText(range_inner, label="Number of increments"), 0)

        # Orbital position (ν)
        grid.Add(wx.StaticText(range_inner, label="Orbital position (Periapse = 0) [°]"),
                0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.txt_nu_min,  1, wx.EXPAND)
        grid.Add(self.txt_nu_max,  1, wx.EXPAND)
        grid.Add(self.sp_nu_n,     0)

        range_inner.SetSizer(grid)
        range_sbs.Add(range_inner, 0, wx.EXPAND | wx.ALL, 6)
        sizer.Add(range_sbs, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 8)

        # Plot button 
        sizer.Add(self.update_btn, 0, wx.ALIGN_CENTER | wx.ALL, 10)


        self.SetSizer(sizer)
        return sizer


    def _on_open_getting_started(self, _evt):
        frame = wx.GetTopLevelParent(self) or self
        try:
            show_getting_started(frame)
        except Exception as e:
            wx.MessageBox(f"Couldn’t open Getting Started:\n{e}",
                          "StressViz", style=wx.OK | wx.ICON_ERROR)


    # ======================= Europa preset wiring =================================
    def _europa_preset_dto(self):
        """Exact numbers matching EuropaSample.sat, in SI; NSR given in seconds."""
        return {
            "system_id": "JupiterEuropa",
            "planet_mass_kg": 1.899e27,
            "ecc": 9.4e-3,
            "a_m": 6.709e8,
            "nsr_years": 3.15576e12 / (365.25 * 86400.0),  # 3.15576e12 s -> years
            "layers": [
                {"layer_id": "ICE_UPPER", "density": 917.0,  "young": 9.28e9,  "poisson": 0.33,  "thickness_m": 1.2e4,  "viscosity": 1.0e22},
                {"layer_id": "ICE_LOWER", "density": 917.0,  "young": 9.28e9,  "poisson": 0.331, "thickness_m": 8.0e3,  "viscosity": 1.0e17},
                {"layer_id": "OCEAN",     "density": 1000.0, "young": 0.0,     "poisson": 0.5,   "thickness_m": 1.0e5,  "viscosity": 0.0},
                {"layer_id": "CORE",      "density": 3847.6, "young": 1.0e11,  "poisson": 0.25,  "thickness_m": 1.391e6,"viscosity": 0.0},
            ],
        }

    def _fmtE(self, x):
        try:
            return f"{float(x):.6E}"
        except Exception:
            return str(x) if x is not None else ""

    def _mirror_dto_into_sat_panel(self, dto):
        """Populate the Satellite panel fields, using scientific notation."""
        sp = self.sat_panel
        if not sp:
            return
        sp.txt_system_id.SetValue(dto.get("system_id") or "")
        sp.txt_planet_mass_kg.SetValue(self._fmtE(dto.get("planet_mass_kg")))
        sp.txt_ecc.SetValue(self._fmtE(dto.get("ecc")))
        sp.txt_a_m.SetValue(self._fmtE(dto.get("a_m")))
        sp.txt_nsr_years.SetValue(self._fmtE(dto.get("nsr_years")))

        lf = getattr(sp, "layer_fields", {})
        for L in ("ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE"):
            rec = next((r for r in dto["layers"] if str(r.get("layer_id")).upper() == L), {})
            def put(key_ui, key_dto):
                ctrl = lf.get((L, key_ui))
                if ctrl:
                    ctrl.SetValue(self._fmtE(rec.get(key_dto)))
            put("rho", "density")
            put("E",   "young")
            put("nu",  "poisson")
            put("H",   "thickness_m")
            put("eta", "viscosity")

    def _on_click_europa_preset(self, _evt):
        """Fill the UI with Europa values AND materialize a SatStress Satellite object."""
        dto = self._europa_preset_dto()
        self._mirror_dto_into_sat_panel(dto)

        # Build/Load a real SatStress Satellite object from the DTO
        try:
            sat = self._dto_to_satellite(dto)
            self._on_satellite_changed(sat)     # propagate like a file load
        except Exception:
            # _dto_to_satellite already messaged if it failed; leave UI populated
            pass

    # ====================== Satellite change / DTO materialization =================
    def _on_satellite_changed(self, satellite):
        satellite = normalize_satellite_layer_order(satellite)
        self.satellite = satellite
        self.compute_love_btn.Enable(self.satellite is not None)

        # If a file path exists, try to read cached Love numbers
        path = getattr(self.sat_panel, "get_loaded_path", lambda: None)()
        if path and read_love_from_sat:
            try:
                found = read_love_from_sat(path)
            except Exception:
                found = None
            if found:
                h2, k2, l2 = found
                def fmt(z):
                    sign = "+" if z.imag >= 0 else "−"
                    return f"{z.real:.6g} {sign} {abs(z.imag):.3g}i"
                self.love_h_input.SetValue(fmt(h2))
                self.love_k_input.SetValue(fmt(k2))
                self.love_l_input.SetValue(fmt(l2))
            else:
                for box in (self.love_h_input, self.love_k_input, self.love_l_input):
                    box.SetValue("")

        # keep orbit panel synced
        if self.orbital_panel and hasattr(self.orbital_panel, "set_satellite"):
            try:
                self.orbital_panel.set_satellite(satellite)
                if hasattr(self.orbital_panel, "redraw"):
                    self.orbital_panel.redraw()
            except Exception:
                pass

    def _get_current_dto(self):
        """Harvest current UI fields into a DTO (works with the text-field matrix)."""
        sp = self.sat_panel

        def getf(ctrl):
            if not ctrl: return None
            s = ctrl.GetValue().strip()
            if not s: return None
            if s.lower() == "infinity": return "infinity"
            try: return float(s.replace(",", ""))
            except Exception: return s

        dto = {
            "system_id": sp.txt_system_id.GetValue().strip() if getattr(sp, "txt_system_id", None) else None,
            "planet_mass_kg": getf(getattr(sp, "txt_planet_mass_kg", None)),
            "ecc":            getf(getattr(sp, "txt_ecc", None)),
            "a_m":            getf(getattr(sp, "txt_a_m", None)),
            "nsr_years":      (sp.txt_nsr_years.GetValue().strip()
                               if getattr(sp, "txt_nsr_years", None) else None),
            "layers": [],
        }

        lf = getattr(sp, "layer_fields", {})
        order = ("ICE_UPPER","ICE_LOWER","OCEAN","CORE")
        def F(layer, key):
            ctrl = lf.get((layer, key))
            if not ctrl: return None
            s = ctrl.GetValue().strip()
            if not s: return None
            try: return float(s.replace(",", ""))
            except Exception: return s

        for L in order:
            dto["layers"].append({
                "layer_id":   L,
                "density":    F(L, "rho"),
                "young":      F(L, "E"),
                "poisson":    F(L, "nu"),
                "thickness_m":F(L, "H"),
                "viscosity":  F(L, "eta"),
            })
        return dto

    def _dto_to_satellite(self, dto):
        """
        Return a SatStress Satellite object by writing several .sat dialect variants
        to a temp folder and trying to load each via load_satellite_from_file().
        """
        import tempfile, pathlib
        from .satstress_interface import load_satellite_from_file

        # ---------- helpers ----------
        def fnum(v):
            if v in (None, "", "infinity"): return None
            try: return float(v)
            except Exception: return None

        def fmtE(v):
            v = fnum(v)
            if v is None: return None
            # Uppercase E, no leading plus
            return f"{v:.6E}".replace("+0", "+").replace("-0", "-")

        def nsr_line(mode: str):
            """
            mode: 'omit' | 'inf' | 'seconds'
            """
            if mode == "omit":
                return None
            if mode == "inf":
                return "NSR_PERIOD = infinity"
            # seconds from years if provided
            y = dto.get("nsr_years")
            if y in (None, "", "infinity"):
                return "NSR_PERIOD = infinity"
            try:
                seconds = float(y) * 365.25 * 86400.0
                return f"NSR_PERIOD = {fmtE(seconds)}"
            except Exception:
                return "NSR_PERIOD = infinity"

        def a_lines(axis_key_variant: str):
            """
            axis_key_variant: 'A' | 'A_M' (some loaders want explicit _M)
            """
            a = fmtE(dto.get("a_m"))
            if axis_key_variant == "A_M":
                return [f"ORBIT_SEMIMAJOR_AXIS_M = {a or ''}".rstrip()]
            return [f"ORBIT_SEMIMAJOR_AXIS = {a or ''}".rstrip()]

        # canonical layer order we show in UI (outer→inner)
        ORDER = ("ICE_UPPER", "ICE_LOWER", "OCEAN", "CORE")

        # collect by name for stable addressing
        byname = {str((L.get("layer_id") or "")).upper(): L for L in (dto.get("layers") or [])}

        def layer_rows(index_scheme: str, prop_keys: str, sep: str, order_variant: str):
            """
            index_scheme: '0' (0..3) | '1' (1..4) | 'sample' (CORE=0,OCEAN=1,ICE_LOWER=2,ICE_UPPER=3)
            prop_keys: 'YOUNGS' | 'YOUNG'  (key spelling)
            sep: '=' | ':'   (assignment token)
            order_variant: 'outer_in' | 'sample' (see index_scheme=sample mapping)
            """
            # index mapping
            if index_scheme == "sample":
                names = ("CORE", "OCEAN", "ICE_LOWER", "ICE_UPPER")  # typical legacy files
                idx_for = {nm: i for i, nm in enumerate(names)}
            else:
                names = ORDER if order_variant == "outer_in" else ORDER[::-1]
                base = 0 if index_scheme == "0" else 1
                idx_for = {nm: (i + base) for i, nm in enumerate(names)}

            # property spellings
            KEY_Y = "YOUNGS_MODULUS" if prop_keys == "YOUNGS" else "YOUNG_MODULUS"

            lines = ["NUM_LAYERS = 4", ""]
            for nm in names:
                L = byname.get(nm, {})
                i = idx_for[nm]
                def LN(k, v):
                    vv = fmtE(v) or ""
                    return f"{k} {sep} {vv}".rstrip()
                lines.append(f"LAYER_ID_{i} {sep} {nm}")
                lines.append(LN(f"DENSITY_{i}",        L.get("density")))
                lines.append(LN(f"{KEY_Y}_{i}",        L.get("young")))
                lines.append(LN(f"POISSONS_RATIO_{i}", L.get("poisson")))
                lines.append(LN(f"THICKNESS_{i}",      L.get("thickness_m")))
                lines.append(LN(f"VISCOSITY_{i}",      L.get("viscosity")))
            return lines

        # ---------- header variants ----------
        sys_id = (dto.get("system_id") or "Preset")
        mass   = fmtE(dto.get("planet_mass_kg"))
        ecc    = fmtE(dto.get("ecc"))

        def header(nsr_mode: str, axis_variant: str):
            hdr = [
                f"SYSTEM_ID = {sys_id}",
                f"PLANET_MASS = {mass or ''}".rstrip(),
                f"ORBIT_ECCENTRICITY = {ecc or ''}".rstrip(),
            ]
            hdr.extend(a_lines(axis_variant))
            nsr = nsr_line(nsr_mode)
            if nsr:
                hdr.append(nsr)
            hdr.append("")  # blank line
            return hdr

        # ---------- write & try ----------
        tmp_dir = pathlib.Path(tempfile.gettempdir()) / "StressViz_presets"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        def write_variant(tag, body_lines):
            p = tmp_dir / f"preset_{tag}.sat"
            p.write_text("\n".join(body_lines) + "\n")
            return str(p)

        variants = []

        # Very common combos first
        variants.append((
            "eq_u0_youngs_outer",
            header("inf", "A") + layer_rows(index_scheme="0", prop_keys="YOUNGS", sep="=", order_variant="outer_in")
        ))
        variants.append((
            "eq_u1_youngs_outer",
            header("inf", "A") + layer_rows(index_scheme="1", prop_keys="YOUNGS", sep="=", order_variant="outer_in")
        ))
        variants.append((
            "eq_u0_young_outer_A_M",
            header("inf", "A_M") + layer_rows(index_scheme="0", prop_keys="YOUNG", sep="=", order_variant="outer_in")
        ))
        # Sample-like (CORE=0..ICE_UPPER=3) with colon sep (seen in some forks)
        variants.append((
            "col_sample_youngs",
            header("inf", "A") + layer_rows(index_scheme="sample", prop_keys="YOUNGS", sep=":", order_variant="sample")
        ))
        # No NSR line at all (some parsers choke on 'infinity')
        variants.append((
            "eq_u0_no_nsr",
            header("omit", "A") + layer_rows(index_scheme="0", prop_keys="YOUNGS", sep="=", order_variant="outer_in")
        ))
        # NSR in seconds (if given)
        variants.append((
            "eq_u0_nsr_seconds",
            header("seconds", "A") + layer_rows(index_scheme="0", prop_keys="YOUNGS", sep="=", order_variant="outer_in")
        ))
        # Reversed layer order (inner→outer) just in case
        variants.append((
            "eq_u0_rev",
            header("inf", "A") + layer_rows(index_scheme="0", prop_keys="YOUNGS", sep="=", order_variant="inner_out")
        ))

        last_exc = None
        tried = []

        for tag, content in variants:
            path = write_variant(tag, content)
            tried.append(path)
            try:
                sat = load_satellite_from_file(path)
                print(f"[StressViz] Loaded preset via: {path}")
                return sat
            except Exception as e:
                last_exc = e
                print(f"[StressViz] Variant failed: {path} -> {e}")

        import wx
        wx.MessageBox(
            "Could not create a SatStress Satellite from the preset.\n\n"
            "Tried these temp files (you can open them via 'Load from file'):\n"
            + "\n".join(tried) + f"\n\nLast error: {last_exc}",
            "Preset .sat creation failed", wx.OK | wx.ICON_ERROR
        )
        raise last_exc or RuntimeError("All .sat variants failed")


    # =========================== Compute Love Numbers =============================
    def on_compute_love(self, _evt):
        """Compute Love numbers, materializing a SatStress Satellite if needed."""
        # Use already-loaded object if present; otherwise synthesize from UI/preset
        sat_obj = self.satellite
        if sat_obj is None or isinstance(sat_obj, dict):
            dto = self._get_current_dto()
            if not dto:
                wx.MessageBox("No satellite inputs available.", "Needs .sat",
                              wx.OK | wx.ICON_WARNING)
                return
            try:
                sat_obj = self._dto_to_satellite(dto)
                self.satellite = sat_obj
            except Exception:
                return  # error already shown

        # Compute via SatStress
        try:
            Diurnal = _import_diurnal()
            diurn = Diurnal(sat_obj)

            # Call whatever method exists
            for m in ("calcLove", "calculate_love", "compute_love", "calc_love"):
                if hasattr(diurn, m):
                    getattr(diurn, m)()
                    break

            love = getattr(diurn, "love", None)

            def _get(v, names):
                if v is None: return None
                if isinstance(v, dict):
                    for n in names:
                        if n in v: return v[n]
                for n in names:
                    if hasattr(v, n): return getattr(v, n)
                return None

            h2 = _get(love, ("h2","h","H2","H"))
            k2 = _get(love, ("k2","k","K2","K"))
            l2 = _get(love, ("l2","l","L2","L"))
            if any(x is None for x in (h2,k2,l2)) and isinstance(love,(list,tuple)) and len(love)>=3:
                h2,k2,l2 = love[0],love[1],love[2]

            def fmt_complex(z):
                try:
                    z = complex(z)
                    sign = "+" if z.imag >= 0 else "−"
                    return f"{z.real:.3g} {sign} {abs(z.imag):.3g}i"
                except Exception:
                    return str(z) if z is not None else ""

            self.love_h_input.SetValue(fmt_complex(h2))
            self.love_k_input.SetValue(fmt_complex(k2))
            self.love_l_input.SetValue(fmt_complex(l2))
            #wx.MessageBox("Love numbers computed.", "Success", wx.OK | wx.ICON_INFORMATION)
        except Exception as e:
            import traceback; traceback.print_exc()
            wx.MessageBox(f"Love number computation failed:\n{e}",
                          "Error", wx.OK | wx.ICON_ERROR)

    # ================= Encounters + resolvers ==================

    def load_encounters_excel(self, path, sheet_name=0, show_message=True):
        """Excel path (kept for compatibility) — read as STRINGS to avoid Timestamp issues.
        Populates self.encounters_by_id, Point panel choices, and pushes meta to ScalarPlotPanel."""

        try:
            df = pd.read_excel(path, sheet_name=sheet_name, dtype=str)
            if isinstance(df, dict):
                df = next(iter(df.values()))
        except Exception as e:
            if show_message:
                wx.MessageBox(f"Failed to read Excel:\n{e}", "Encounters", wx.OK | wx.ICON_ERROR)
            return False

        COL_ID  = "Encounter Tag"
        COL_ET  = "Time of C/A (TCA) (ET/TDB)"
        COL_NU  = "True Anomaly of Body (deg)"
        COL_LAT = "Latitude Planetocentric (deg)"
        COL_LON = "E. Longitude (deg)"

        missing = [c for c in (COL_ID, COL_ET, COL_NU, COL_LAT, COL_LON) if c not in df.columns]
        if missing:
            if show_message:
                wx.MessageBox("Excel is missing columns:\n  " + "\n  ".join(missing),
                            "Encounters", wx.OK | wx.ICON_ERROR)
            return False

        # Safe caster (avoids float(Timestamp))
        def _num_or_none(x):
            try:
                if x is None:
                    return None
                s = str(x).strip()
                if not s or s.lower() == "nan":
                    return None
                return float(s)
            except Exception:
                return None

        # Eccentricity
        try:
            e = float(self.get_eccentricity())
        except Exception:
            e = 0.0

        # Optional: ET→UTC string
        def _et_to_utc_iso_safe(et_sec):
            try:
                if et_sec is None:
                    return None
                from astropy.time import Time, TimeDelta
                t0 = Time('2000-01-01T12:00:00', scale='tdb')
                t  = t0 + TimeDelta(float(et_sec), format='sec')
                return t.utc.isot
            except Exception:
                return None

        self.encounters_by_id = {}
        et_nu_pairs = []

        slim = df[[COL_ID, COL_ET, COL_NU, COL_LAT, COL_LON]].copy()

        # Mission-aware prefix (default "Clipper")
        MISSION_PREFIX = getattr(self, "mission_prefix", "Clipper")

        for _, r in slim.iterrows():
            # original tag from the sheet (e.g., "E12" or "12")
            orig = (str(r[COL_ID]).strip() if pd.notna(r[COL_ID]) else "")
            if not orig:
                continue

            # display/primary ID with mission prefix
            enc = f"{MISSION_PREFIX}-{orig}"

            et  = _num_or_none(r[COL_ET])
            nu  = _num_or_none(r[COL_NU])
            lat = _num_or_none(r[COL_LAT])
            lon = _num_or_none(r[COL_LON])
            utc_iso = _et_to_utc_iso_safe(et)
            M = _true2mean(nu, e) if nu is not None else None

            self.encounters_by_id[enc] = {
                "encounter_id": enc,       # e.g., "Clipper-E12"
                "encounter_tag": orig,     # keep the raw tag too (nice for filenames, etc.)
                "et_tdb_sec": et,
                "utc_iso": utc_iso,
                "true_anom_deg": nu,
                "mean_anom_deg": M,
                "lat_deg": lat,
                "lon_deg": lon,
                # period_hours can be added elsewhere; Scalar push will default if missing
            }

            if et is not None and nu is not None:
                et_nu_pairs.append((et, nu))

        et_nu_pairs.sort(key=lambda x: x[0])
        self._et_list  = [t for t, _ in et_nu_pairs]
        self._nu_by_et = dict(et_nu_pairs)

        # ---- Update Point panel choices ----
        ids = sorted(self.encounters_by_id.keys())
        if hasattr(self.point_panel, "set_encounter_choices"):
            try:
                self.point_panel.set_encounter_choices(ids, select=(ids[0] if ids else None))
            except Exception:
                traceback.print_exc()

        # ---- Also push to ScalarPlotPanel so it gets M_ca + period (needed for time nudges) ----
        try:
            self._push_encounters_to_scalar_panel(select_id=(ids[0] if ids else None))
        except Exception:
            traceback.print_exc()

        return bool(ids)



    def _auto_load_default_encounters(self):
        """Prefer hardened TXT; fall back to Excel. Then append plume epochs.
        Populates Point panel choices and pushes to ScalarPlotPanel.
        """
        import os
        import traceback

        print("\n[autoload] ===== starting _auto_load_default_encounters =====")

        MISSION_PREFIX = getattr(self, "mission_prefix", "Clipper")
        P_H = getattr(self, "period_hours", 85.228)
        print(f"[autoload] mission_prefix={MISSION_PREFIX!r}, period_hours={P_H!r}")

        # Eccentricity for ν→M; try scalar panel, else 0.
        try:
            e = float(getattr(getattr(self, "scalar_panel_ref", None), "_get_eccentricity")())
            print(f"[autoload] eccentricity from scalar_panel_ref: {e}")
        except Exception:
            try:
                e = float(self.get_eccentricity())
                print(f"[autoload] eccentricity from self.get_eccentricity(): {e}")
            except Exception:
                e = 0.0
                print("[autoload] eccentricity fallback to 0.0")

        loaded = False

        # -------- TXT path (recommended) --------
        try:
            try:
                from .encounters_io import load_europa_min
                _loader = load_europa_min
                print("[autoload] using loader: .encounters_io.load_europa_min")
            except Exception:
                from .encounters import load_europa_min  # type: ignore
                _loader = load_europa_min
                print("[autoload] using loader: .encounters.load_europa_min")

            try_paths_txt = [
                os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "Europa_Encounters_21F31_V7_LP01_ver2.txt"),
                os.path.join(os.getcwd(), "data", "Europa_Encounters_21F31_V7_LP01_ver2.txt"),
                os.path.join(os.getcwd(), "Europa_Encounters_21F31_V7_LP01_ver2.txt"),
            ]
            print("[autoload] TXT candidate paths:")
            for p in try_paths_txt:
                print(f"    {p}  exists={os.path.isfile(p)}")

            for p in try_paths_txt:
                if os.path.isfile(p):
                    print(f"[autoload] found TXT file: {p}")
                    try:
                        df = _loader(p)
                        print(f"[autoload] TXT load succeeded: len(df)={len(df)}")
                        try:
                            print(f"[autoload] TXT columns={list(df.columns)}")
                        except Exception:
                            pass
                        try:
                            print("[autoload] TXT head:")
                            print(df.head())
                        except Exception:
                            pass

                        self.encounters_by_id = {}
                        et_nu_pairs = []

                        for i, (_, r) in enumerate(df.iterrows()):
                            raw = str(r["enc_label"]).strip()  # e.g., "12" or "00"
                            orig = f"E{raw}" if raw and not raw.startswith("E") else raw
                            enc = f"{MISSION_PREFIX}-{orig}"

                            et = r.get("tca_et_sec")
                            nu = r.get("true_anom_deg")
                            lat = r.get("lat_pc_deg")
                            lon = r.get("lon_e_deg")

                            try:
                                et = float(et) if et is not None else None
                                nu = float(nu) if nu is not None else None
                                lat = float(lat) if lat is not None else None
                                lon = float(lon) if lon is not None else None
                            except Exception:
                                pass

                            rec = {
                                "encounter_id": enc,
                                "encounter_tag": orig,
                                "et_tdb_sec": et,
                                "utc_iso": None,
                                "true_anom_deg": nu,
                                "mean_anom_deg": None,
                                "lat_deg": lat,
                                "lon_deg": lon,
                                "period_hours": P_H,
                            }

                            if nu is not None:
                                try:
                                    rec["mean_anom_deg"] = float(_true2mean(nu, e)) % 360.0
                                except Exception:
                                    print(f"[autoload] warning: _true2mean failed for {enc} with nu={nu}")
                                    traceback.print_exc()

                            self.encounters_by_id[enc] = rec
                            if et is not None and nu is not None:
                                et_nu_pairs.append((et, nu))

                            if i < 5:
                                print(f"[autoload] sample rec {i}: {enc} -> {rec}")

                        et_nu_pairs.sort(key=lambda x: x[0])
                        self._et_list = [t for t, _ in et_nu_pairs]
                        self._nu_by_et = dict(et_nu_pairs)

                        print(f"[autoload] built encounters_by_id count={len(self.encounters_by_id)}")
                        print(f"[autoload] first ids={list(self.encounters_by_id.keys())[:10]}")
                        print(f"[autoload] et_nu_pairs count={len(et_nu_pairs)}")

                        loaded = True
                    except Exception:
                        print("[autoload] TXT load/build failed:")
                        traceback.print_exc()
                    break
        except Exception:
            print("[autoload] outer TXT block failed:")
            traceback.print_exc()

        # -------- Excel fallback --------
        if not loaded:
            print("[autoload] TXT path did not load encounters, trying Excel fallback")
            try_paths_xlsx = [
                os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "Europa_Encounters.xlsx"),
                os.path.join(os.getcwd(), "data", "Europa_Encounters.xlsx"),
                os.path.join(os.getcwd(), "Europa_Encounters.xlsx"),
            ]
            print("[autoload] XLSX candidate paths:")
            for p in try_paths_xlsx:
                print(f"    {p}  exists={os.path.isfile(p)}")

            for p in try_paths_xlsx:
                if os.path.isfile(p):
                    print(f"[autoload] found XLSX file: {p}")
                    try:
                        self.load_encounters_excel(p, show_message=False)
                        print(f"[autoload] Excel load finished; raw encounters_by_id count={len(getattr(self, 'encounters_by_id', {}) or {})}")

                        if getattr(self, "encounters_by_id", None):
                            rebased = {}
                            for k, rec in self.encounters_by_id.items():
                                orig = k if not str(k).startswith(f"{MISSION_PREFIX}-") else str(k).split("-", 1)[1]
                                newk = f"{MISSION_PREFIX}-{orig}"
                                rec = rec.copy()
                                rec["encounter_id"] = newk
                                rec.setdefault("encounter_tag", orig)
                                rec.setdefault("period_hours", P_H)
                                if rec.get("mean_anom_deg") is None and rec.get("true_anom_deg") is not None:
                                    try:
                                        rec["mean_anom_deg"] = float(_true2mean(float(rec["true_anom_deg"]), e)) % 360.0
                                    except Exception:
                                        print(f"[autoload] warning: Excel _true2mean failed for {newk}")
                                        traceback.print_exc()
                                rebased[newk] = rec
                            self.encounters_by_id = rebased

                        print(f"[autoload] rebased encounters_by_id count={len(getattr(self, 'encounters_by_id', {}) or {})}")
                        print(f"[autoload] first ids={list(getattr(self, 'encounters_by_id', {}).keys())[:10]}")
                        loaded = True
                    except Exception:
                        print("[autoload] Excel fallback failed:")
                        traceback.print_exc()
                    break

        # -------- Add plume observation epochs --------
        # try:
        #     n_plume = self._try_autoload_plume_observations(e=e, P_H=P_H)
        #     print(f"[plume] helper added {n_plume} plume encounters")
        # except Exception as _plume_err:
        #     print(f"[plume] helper raised {_plume_err}")
        #     traceback.print_exc()

        # -------- Push choices to UI --------
        ids = sorted(getattr(self, "encounters_by_id", {}).keys())
        print(f"[autoload] final ids count={len(ids)}")
        print(f"[autoload] final ids sample={ids[:10]}")
        print(f"[autoload] point_panel={getattr(self, 'point_panel', None)}")

        if hasattr(self.point_panel, "set_encounter_choices"):
            try:
                print("[autoload] calling point_panel.set_encounter_choices(...)")
                self.point_panel.set_encounter_choices(ids, select=(ids[0] if ids else None))
            except Exception:
                print("[autoload] point_panel.set_encounter_choices failed:")
                traceback.print_exc()

        if hasattr(self, "encounter_choice") and self.encounter_choice:
            try:
                print("[autoload] updating encounter_choice widget")
                self.encounter_choice.SetItems(ids)
                self.encounter_choice.Enable(bool(ids))
                if ids:
                    self.encounter_choice.SetSelection(0)
            except Exception:
                print("[autoload] updating encounter_choice failed:")
                traceback.print_exc()

        # -------- Also push to ScalarPlotPanel --------
        try:
            print("[autoload] pushing encounters to scalar panel")
            self._push_encounters_to_scalar_panel(select_id=(ids[0] if ids else None))
        except Exception:
            print("[autoload] _push_encounters_to_scalar_panel failed:")
            traceback.print_exc()

        print("[autoload] ===== done _auto_load_default_encounters =====\n")


    def _push_encounters_to_scalar_panel(self, select_id=None):
        """Send current encounters to the ScalarPlotPanel dropdown (GUI thread-safe),
        with verbose debug + hard defaults so M_ca/period/utc are never None."""
        

        def _dbg(msg, *a):
            print("[StressViz/_push_encounters] " + (msg % a if a else msg))

        try:
            sp = getattr(self, "scalar_panel_ref", None)
            if not sp:
                self._pending_push_to_scalar = True
                self._pending_scalar_select_id = select_id
                _dbg("No scalar_panel_ref; deferring push until scalar panel is created.")
                return
            
            self._pending_push_to_scalar = False
            self._pending_scalar_select_id = None
            
            if not getattr(self, "encounters_by_id", None):
                _dbg("No encounters_by_id on controller; aborting.")
                return

            def _parse_iso_z(s):
                if not s:
                    return None
                s = str(s).strip()
                if s.endswith("Z"):
                    s = s[:-1]
                try:
                    return _dt.datetime.fromisoformat(s).replace(tzinfo=_dt.timezone.utc)
                except Exception:
                    return None

            # Eccentricity for ν→M fallback
            try:
                e = float(sp._get_eccentricity()) if hasattr(sp, "_get_eccentricity") else 0.0
            except Exception:
                e = 0.0

            # Robust default period (controller -> panel -> hardcoded)
            P_H_DEFAULT = (
                getattr(self, "period_hours", None)
                or getattr(sp, "period_hours", None)
                or 85.228  # Europa default (hours)
            )

            # Hard default CA time if missing
            UTC_DEFAULT = _dt.datetime(2000, 1, 1, tzinfo=_dt.timezone.utc)

            ids = sorted(self.encounters_by_id.keys(), key=str)
            _dbg("Preparing %d encounters; e=%.6g, default P_h=%.6g", len(ids), e, P_H_DEFAULT)

            meta_by_id = {}
            for eid in ids:
                rec = self.encounters_by_id.get(eid) or {}

                # ---- CA time ----
                utc_ca = rec.get("utc_dt") or _parse_iso_z(rec.get("utc_iso") or rec.get("utc")) or UTC_DEFAULT

                # ---- Mean anomaly at CA ----
                M_ca = rec.get("mean_anom_deg")
                if M_ca is None:
                    # 1) Try ν stored on the record (normal Clipper path)
                    nu = rec.get("true_anom_deg") or rec.get("nu_deg")
                    if nu is not None:
                        try:
                            M_ca = float(_true2mean(float(nu), float(e))) % 360.0
                        except Exception:
                            _dbg("eid=%s: ν→M conversion failed (nu=%r, e=%r)", eid, nu, e)
                            M_ca = None

                # 2) If we STILL don't have M, try resolving ν from the UTC time (plume path)
                if M_ca is None and utc_ca is not None:
                    resolver = getattr(self, "_resolve_nu_from_time", None)
                    if callable(resolver):
                        try:
                            # Make a clean UTC string, like the Point panel passes to resolve_nu_from_time
                            if isinstance(utc_ca, _dt.datetime):
                                utc_str = utc_ca.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
                            else:
                                utc_str = str(utc_ca)

                            nu_res = resolver(utc_str, is_et=False)
                            if nu_res is not None:
                                nu_f = float(nu_res) % 360.0
                                try:
                                    M_ca = float(_true2mean(nu_f, float(e))) % 360.0
                                except Exception:
                                    # If ν→M fails for some reason, just approximate M ≈ ν
                                    M_ca = nu_f
                                _dbg("eid=%s: resolved ν from time: nu=%.6f → M=%.6f",
                                    eid, nu_f, M_ca)
                        except Exception as err:
                            _dbg("eid=%s: time-based ν resolve failed: %r", eid, err)

                # 3) Normalize M if we have it
                if M_ca is not None:
                    try:
                        M_ca = float(M_ca) % 360.0
                    except Exception:
                        _dbg("eid=%s: mean_anom_deg not floaty: %r", eid, M_ca)
                        M_ca = None

                # 4) Final hard default so the panel never sees None
                if M_ca is None:
                    M_ca = 0.0

                # ---- Back-propagate M_ca into encounters_by_id so other code can use it ----
                base_rec = self.encounters_by_id.get(eid)
                if isinstance(base_rec, dict):
                    # Make a shallow copy so we don't accidentally share references
                    base_rec = base_rec.copy()

                    # Preserve any ν we already had
                    if base_rec.get("true_anom_deg") is None and base_rec.get("nu_deg") is not None:
                        base_rec["true_anom_deg"] = base_rec["nu_deg"]

                    base_rec["mean_anom_deg"] = M_ca
                    self.encounters_by_id[eid] = base_rec


                # ---- Period hours ----
                P_h = rec.get("period_hours", P_H_DEFAULT)
                try:
                    P_h = float(P_h)
                    if not (P_h > 0):
                        raise ValueError("non-positive period")
                except Exception:
                    _dbg("eid=%s: invalid period_hours=%r; using default %.6g", eid, P_h, P_H_DEFAULT)
                    P_h = float(P_H_DEFAULT)

                meta_by_id[eid] = {
                    "utc_ca": utc_ca,
                    "M_ca_deg": M_ca,
                    "period_hours": P_h,
                }

                _dbg(
                    "eid=%s -> utc_ca=%s | M_ca_deg=%.6f | period_hours=%.6f",
                    eid,
                    utc_ca.isoformat() if isinstance(utc_ca, _dt.datetime) else str(utc_ca),
                    float(M_ca),
                    float(P_h),
                )

            # Register (GUI thread)
            if hasattr(sp, "register_encounters"):
                wx.CallAfter(sp.register_encounters, meta_by_id)
            else:
                _dbg("Scalar panel missing register_encounters; cannot set M_ca/period metadata.")

            # Then wire the list/selection API
            if hasattr(sp, "set_encounters_from_ids"):
                wx.CallAfter(
                    sp.set_encounters_from_ids,
                    ids,
                    lambda eid: self.encounters_by_id.get(eid),
                    lambda eid, rec: str(eid),
                    select_id,
                )
            else:
                encs = []
                for eid in ids:
                    rec = (self.encounters_by_id.get(eid) or {}).copy()
                    encs.append({
                        "id": eid,
                        "name": str(eid),
                        "nu_deg": rec.get("nu_deg") or rec.get("true_anom_deg"),
                        "M_deg": rec.get("mean_anom_deg"),
                        "utc_iso": rec.get("utc_iso") or rec.get("utc"),
                        "lat_deg": rec.get("lat_deg") or rec.get("lat"),
                        "lon_deg": rec.get("lon_deg") or rec.get("lon"),
                    })
                wx.CallAfter(sp.set_encounters, encs, select_id)

            # Optional confirmation ping (async)
            def _after_register_ping():
                try:
                    cur = getattr(sp, "_selected_encounter_id", None)
                    if cur:
                        m = getattr(sp, "_enc_meta_by_id", {}).get(cur, {})
                        _dbg("ScalarPanel accepted id=%s: M_ca_deg=%r, P_h=%r, utc_ca=%r",
                            cur, m.get("M_ca_deg"), m.get("period_hours"), m.get("utc_ca"))
                except Exception:
                    traceback.print_exc()

            wx.CallAfter(_after_register_ping)

        except Exception:
            traceback.print_exc()


    # --- resolvers used by point panel ---
    def _resolve_encounter(self, enc_id: str):
        return getattr(self, "encounters_by_id", {}).get(enc_id)

    '''
    def _resolve_true_anomaly(self, utc_iso: str) -> float | None:
        """
        Resolve Europa true anomaly ν (deg) at a UTC time using JPL Horizons.
        Stores diagnostics in self._last_horizons_diag.
        """
        self._last_horizons_diag = {
            "ok": False,
            "error": None,
            "stage": None,
            "utc_iso": utc_iso,
            "jd_tdb": None,
            "target_id": None,
            "center_id": None,
            "colnames": None,
            "value": None,
        }

        # Defaults that tend to work for Europa around Jupiter barycenter.
        # Override by setting self.horizons_target_id / self.horizons_center_id if you want.
        target_id = getattr(self, "horizons_target_id", "502")      # Europa (major body id)
        center_id = getattr(self, "horizons_center_id", "500@5")    # Jupiter system barycenter

        self._last_horizons_diag["target_id"] = target_id
        self._last_horizons_diag["center_id"] = center_id

        try:
            from astropy.time import Time
            from astroquery.jplhorizons import Horizons
        except Exception as e:
            self._last_horizons_diag["error"] = f"Missing deps (astropy/astroquery): {e}"
            self._last_horizons_diag["stage"] = "import"
            return None

        # Parse UTC -> JD(TDB)
        try:
            t_utc = self._utc_string_to_astropy_time(str(utc_iso))
            jd_tdb = float(t_utc.tdb.jd)
            self._last_horizons_diag["jd_tdb"] = jd_tdb
        except Exception as e:
            self._last_horizons_diag["error"] = f"Bad UTC '{utc_iso}': {e}"
            self._last_horizons_diag["stage"] = "time_parse"
            return None

        # Horizons query helper
        def _extract_nu_from_table(tab) -> float | None:
            try:
                colnames = list(getattr(tab, "colnames", []) or [])
                self._last_horizons_diag["colnames"] = colnames
            except Exception:
                pass

            # Try common true anomaly column names
            for key in ("TA", "true_anom", "true_anomaly", "TruAnom", "tru_anom"):
                if hasattr(tab, "colnames") and key in tab.colnames:
                    try:
                        return float(tab[key][0])
                    except Exception:
                        pass
            return None

        # 1) Try elements() (sometimes works depending on target/center)
        try:
            obj = Horizons(id=target_id, id_type=None, location=center_id, epochs=[jd_tdb])
            tab = obj.elements()
            self._last_horizons_diag["stage"] = "elements"
            nu = _extract_nu_from_table(tab)
            if nu is not None:
                nu = float(nu) % 360.0
                self._last_horizons_diag.update(ok=True, value=nu)
                return nu
        except Exception as e:
            self._last_horizons_diag["error"] = f"Horizons elements() failed: {e}"
            self._last_horizons_diag["stage"] = "elements_error"

        # 2) Try ephemerides() (often more reliable for major bodies)
        try:
            obj = Horizons(id=target_id, id_type=None, location=center_id, epochs=[jd_tdb])
            tab = obj.ephemerides()
            self._last_horizons_diag["stage"] = "ephemerides"
            nu = _extract_nu_from_table(tab)
            if nu is not None:
                nu = float(nu) % 360.0
                self._last_horizons_diag.update(ok=True, value=nu)
                return nu

            # If we got here, Horizons returned successfully but didn't include TA.
            self._last_horizons_diag["error"] = (
                "Horizons returned ephemerides but no true anomaly column was found. "
                f"Columns: {getattr(tab, 'colnames', None)}"
            )
            self._last_horizons_diag["stage"] = "ephemerides_no_TA"
            return None

        except Exception as e:
            self._last_horizons_diag["error"] = f"Horizons ephemerides() failed: {e}"
            self._last_horizons_diag["stage"] = "ephemerides_error"
            return None
    '''


    def get_last_horizons_diag(self):
        return getattr(self, "_last_horizons_diag", None)
    '''
    def _resolve_true_anomaly_from_et(self, et_sec: float) -> Optional[float]:
        if not hasattr(self, "_et_list") or not self._et_list:
            return None
        import bisect
        et = float(et_sec)
        i = bisect.bisect_left(self._et_list, et)
        cands = []
        if i > 0: cands.append(self._et_list[i-1])
        if i < len(self._et_list): cands.append(self._et_list[i])
        if not cands: return None
        best = min(cands, key=lambda t: abs(t - et))
        return self._nu_by_et.get(best)

    def _resolve_nu_from_time(self, time_val, *, is_et: bool, allow_fallback: bool = True):
        """
        Resolve ν (true anomaly, deg). Prefer Horizons; fallback to nearest encounter ν.

        IMPORTANT: Forces center to Jupiter planet (@599) for consistency with encounter geometry.
        """
        from astropy.time import Time

        self._last_horizons_diag = {"ok": False, "error": None, "fallback": None}

        target_id = getattr(self, "horizons_target_id", _HORIZONS_TARGET)

        # Force Jupiter *planet* center for consistency (avoid @5 barycenter ambiguity)
        center_id = "@599"

        # --- parse/convert time to jd_tdb and et_sec ---
        try:
            if is_et:
                et_sec = float(time_val)
                jd_tdb = self._etsec_to_jd_tdb(et_sec)
            else:
                t_utc = self._utc_string_to_astropy_time(str(time_val))
                jd_tdb = float(t_utc.tdb.jd)
                et_sec = float((t_utc.tdb - Time("J2000", scale="tdb")).to_value("s"))
        except Exception as e:
            self._last_horizons_diag.update(error=f"time parse/convert failed: {e}")
            return None

        self._last_horizons_diag.update(
            is_et=bool(is_et),
            et_sec=et_sec,
            jd_tdb=jd_tdb,
            target_id=str(target_id),
            center_id=str(center_id),
        )

        # --- 1) Horizons elements() TA ---
        try:
            from astroquery.jplhorizons import Horizons
            obj = Horizons(id=target_id, id_type="id", location=center_id, epochs=[jd_tdb])
            tab = obj.elements()

            nu = None
            for key in ("TA", "true_anom", "true_anomaly"):
                if key in tab.colnames:
                    nu = float(tab[key][0])
                    break

            if nu is not None:
                nu = float(nu) % 360.0
                self._last_horizons_diag.update(ok=True, nu_deg=nu, src="horizons_elements_TA")
                return nu

            self._last_horizons_diag.update(
                error=f"No true anomaly column in Horizons return. colnames={list(tab.colnames)}",
                src="horizons_elements",
            )

        except Exception as e:
            self._last_horizons_diag.update(error=f"Horizons request failed: {e}", src="horizons_elements")

        # --- 2) Fallback: nearest encounter ν (by ET) ---
        if allow_fallback:
            nu_fb = self._resolve_true_anomaly_from_et(et_sec)
            if nu_fb is not None:
                nu_fb = float(nu_fb) % 360.0
                self._last_horizons_diag.update(ok=True, fallback="encounter_nearest", nu_deg=nu_fb, src="fallback_encounter")
                return nu_fb

        return None
    '''

    # --- canonical time conversions ---
    def _utc_to_jd_tdb(self, utc_iso: str) -> float:
        t = self._utc_string_to_astropy_time(utc_iso)
        return float(t.tdb.jd)

    def _etsec_to_jd_tdb(self, et_sec: float) -> float:
        return 2451545.0 + float(et_sec) / 86400.0

    def _utc_string_to_astropy_time(self, utc_iso: str):
        from astropy.time import Time
        s = str(utc_iso).strip().replace(" ", "T")
        if s.endswith("Z"):
            s = s[:-1] + "+00:00"
        dt = _dt.datetime.fromisoformat(s).astimezone(timezone.utc)
        return Time(dt, scale="utc")

    # ====================== Point compute / Map ===================

    def _on_point_compute(self, *, lat, lon, t, sat, Diurnal, nu=None, **_):
        """
        Compute point stresses using SatStress conventions.
        Returns Pa-level values in a dict.
        """
        import math, numpy as np

        if sat is None or Diurnal is None:
            return None

        diurn = Diurnal(sat)
        holder = getattr(diurn, "stresses", diurn)
        for k, v in (
            ("diurnal", True),
            ("tidal",   False),
            ("eccentricity", True),
            ("ecc",     True),
            ("nsr",     False),
            ("polar_wander", False),
            ("pw",      False),
            ("obliquity", False),
        ):
            if hasattr(holder, k):
                try: setattr(holder, k, bool(v))
                except Exception: pass

        lat_deg = float(lat)
        lon_deg = float(lon)
        theta = math.radians(90.0 - lat_deg)
        phi   = math.radians((lon_deg % 360.0))

        # mean motion
        omega = None
        for key in ("omega", "n", "mean_motion"):
            if hasattr(diurn, key):
                try:
                    val = float(getattr(diurn, key))
                    if np.isfinite(val) and val != 0.0:
                        omega = val; break
                except Exception:
                    pass
        if omega is None:
            return None

        # time
        if nu is not None:
            t_sec = (math.radians(float(nu)) % (2.0*math.pi)) / omega
        else:
            if t is None:
                return None
            t_sec = float(t)

        StressCalc = _import_stresscalc()
        calc = StressCalc([diurn])

        Ttt, Tpt, Tpp = calc.tensor(theta, phi, t_sec)
        import numpy as _np
        Ttt = float(_np.asarray(Ttt).reshape(()))
        Tpt = float(_np.asarray(Tpt).reshape(()))
        Tpp = float(_np.asarray(Tpp).reshape(()))

        tens_mag, tens_az, comp_mag, comp_az = calc.principal_components(theta, phi, t_sec)
        s1 = float(_np.asarray(tens_mag).reshape(()))
        s3 = float(_np.asarray(comp_mag).reshape(()))
        az_cw_from_n = float(_np.asarray(tens_az).reshape(()))

        alpha_deg = (90.0 - math.degrees(az_cw_from_n)) % 180.0

        try:
            self._last_point_inputs = {
                "lat": float(lat_deg),
                "lon": float(lon_deg),
                "nu": None if nu is None else float(nu),
                "utc": t if (t is not None and not isinstance(t, (float, int))) else None,
                "et":  None if (t is None or isinstance(t, str)) else float(t),
            }
            if getattr(self, "update_btn", None):
                self.update_btn.Enable(True)
        except Exception:
            pass

        return {
            "Ttt": Ttt, "Tpt": Tpt, "Tpp": Tpp,
            "sigma1": s1, "sigma3": s3,
            "alpha_deg": alpha_deg,
        }

    def on_open_map(self, _evt):
        """Open the ScalarPlotPanel popup and wire it to current controls."""
        if not getattr(self, "satellite", None):
            wx.MessageBox("Load/apply a satellite first.", "StressViz map",
                        wx.OK | wx.ICON_INFORMATION)
            return

        # Import Diurnal
        try:
            Diurnal = _import_diurnal()
        except Exception as e:
            wx.MessageBox(f"Couldn't import SatStress Diurnal:\n{e}",
                        "StressViz map", wx.OK | wx.ICON_ERROR)
            return


        # ---------- read range inputs ----------
        def _f(ctrl, default):
            try: return float(ctrl.GetValue())
            except Exception: return float(default)

        def _i(ctrl, default):
            try: return int(ctrl.GetValue())
            except Exception: return int(default)

        lat_min = _f(self.txt_lat_min, -90.0)
        lat_max = _f(self.txt_lat_max,  90.0)
        if lat_min > lat_max: lat_min, lat_max = lat_max, lat_min
        lat_min = max(-90.0, min( 90.0, lat_min))
        lat_max = max(-90.0, min( 90.0, lat_max))

        lon_min = _f(self.txt_lon_min, -180.0)
        lon_max = _f(self.txt_lon_max,  180.0)
        if lon_min > lon_max: lon_min, lon_max = lon_max, lon_min
        lon_min = max(-180.0, min( 180.0, lon_min))
        lon_max = max(-180.0, min( 180.0, lon_max))

        nu_min  = _f(self.txt_nu_min, 0.0)
        nu_max  = _f(self.txt_nu_max, 360.0)
        if nu_min > nu_max: nu_min, nu_max = nu_max, nu_min
        n_nu    = max(2, _i(self.sp_nu_n, 10))

        nlat_vec = max(2, _i(self.sp_lat_n, 10))
        nlon_vec = max(2, _i(self.sp_lon_n, 10))

        N_LAT = 91
        N_LON = 181
        lats  = np.linspace(lat_min, lat_max, N_LAT)
        lons  = np.linspace(lon_min, lon_max, N_LON)
        nus   = np.linspace(nu_min,  nu_max,  n_nu)

        # ---------- build / show panel ----------
        frame, panel = get_or_create_scalar_popup(self)  # parent = control panel
        self._push_encounters_to_scalar_panel(
            select_id=getattr(self, "get_selected_encounter_id", lambda: None)()
        )
        # Hand getters so Save-series can name "<SystemID>_<EncounterID>"
        try: panel._get_system_id = self.sat_panel.get_system_id
        except Exception: panel._get_system_id = lambda: "System"
        try: panel._get_encounter_id = self.point_panel.get_selected_encounter_id
        except Exception: panel._get_encounter_id = lambda: None

        try:
            panel.set_fixed_stress_range_kpa(-100.0, 100.0)
        except Exception:
            pass

        if hasattr(panel, "set_vector_density"):
            panel.set_vector_density(nlat_vec, nlon_vec)
        else:
            panel.vec_n_lat = nlat_vec
            panel.vec_n_lon = nlon_vec

        # ---------- stress evaluator ----------
        try:
            sat = self.satellite
            diurn = Diurnal(sat)

            omega = None
            for key in ("omega", "n", "mean_motion"):
                if hasattr(diurn, key):
                    try:
                        val = float(getattr(diurn, key))
                        if np.isfinite(val) and val != 0.0:
                            omega = val; break
                    except Exception:
                        pass
            if omega is None:
                raise RuntimeError("Diurnal has no valid mean motion (omega/n/mean_motion).")

            thetas = np.radians(90.0 - lats)
            phis   = np.mod(np.radians(lons), 2.0 * np.pi)

            StressCalc = _import_stresscalc()

            def eval_fn(nu_deg: float):
                import numpy as _np
                t_sec = (_np.radians(float(nu_deg)) % (2.0 * np.pi)) / omega

                holder = getattr(diurn, "stresses", diurn)
                def _set(k, v):
                    if hasattr(holder, k):
                        try: setattr(holder, k, bool(v))
                        except Exception: pass
                _set("eccentricity", True); _set("ecc", True)
                _set("diurnal", True); _set("tidal", True)
                _set("nsr", False); _set("non_synchronous_rotation", False)
                _set("polar_wander", False); _set("pw", False)
                _set("obliquity", False)

                calc = StressCalc([diurn])

                Ny = thetas.size; Nx = phis.size
                Ttt = _np.empty((Ny, Nx), float)
                Tpt = _np.empty((Ny, Nx), float)
                Tpp = _np.empty((Ny, Nx), float)
                for i in range(Ny):
                    th = float(thetas[i])
                    for j in range(Nx):
                        ph = float(phis[j])
                        a, b, c = calc.tensor(th, ph, t_sec)
                        Ttt[i, j] = float(_np.asarray(a).reshape(()))
                        Tpt[i, j] = float(_np.asarray(b).reshape(()))
                        Tpp[i, j] = float(_np.asarray(c).reshape(()))
                return Ttt, Tpt, Tpp

            # initial ν
            nu0 = None
            try:
                if hasattr(self, "_last_point_inputs"):
                    nu0 = self._last_point_inputs.get("nu")
                if nu0 is None and hasattr(self.point_panel, "txt_nu"):
                    s = self.point_panel.txt_nu.GetValue().strip()
                    nu0 = float(s) if s else None
            except Exception:
                nu0 = None
            if nu0 is not None:
                nu0 = float(np.clip(nu0, nu_min, nu_max))

            panel.set_axes(lats, lons)
            panel.bind_orbit_series(
                nus_deg=list(map(float, nus)),
                eval_fn=eval_fn,
                initial_nu_deg=(float(nu0) if nu0 is not None else None),
                title="StressViz Map",
            )

        except Exception as e:
            wx.MessageBox(f"Failed to initialize map:\n{e}",
                        "StressViz map", wx.OK | wx.ICON_ERROR)
            try:
                frame.Destroy()
            except Exception:
                pass
































