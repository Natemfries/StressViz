# stressviz/point_panel.py
from __future__ import annotations

import wx
import numpy as np
from astropy.time import Time
from typing import Optional, Callable, Iterable, Union
from .ui_style import instr_label, style_staticbox, apply_font
from .utils import true_to_mean_anomaly_deg  
from .horizons_phase import europa_phase_from_horizons

__all__ = ["PointStressPanel"]

Float = Optional[float]
_KPA = 1e-3  # Pa -> kPa


def _try_float(s: str) -> Float:
    try:
        return float(str(s).strip())
    except Exception:
        return None
    
# Robust importer for SatStress StressCalc (matches the map path)
def _import_stresscalc():
    try:
        from SatStress.satstress.stresscalc import StressCalc  # pip/namespace layout
        return StressCalc
    except Exception:
        pass
    try:
        from satstress.stresscalc import StressCalc  # flat package layout
        return StressCalc
    except Exception:
        pass

    # Last-ditch: import module then grab attribute if present
    import importlib
    for modpath in ("SatStress.satstress", "satstress.stresscalc", "SatStress.satstress.stresscalc"):
        try:
            mod = importlib.import_module(modpath)
            if hasattr(mod, "StressCalc"):
                return getattr(mod, "StressCalc")
        except Exception:
            pass

    raise ImportError("Could not import SatStress StressCalc.")



class PointStressPanel(wx.Panel):
    """
    Point-wise stress input + results panel (SatStress-style).

    This class builds the UI content but does NOT pack it into its own root sizer.
    Instead, use:
        - mount_inputs_into(parent_panel)
        - mount_results_into(parent_panel)
    to place the two groups where you want (e.g., side-by-side in AnalysisControlPanel).

    Inputs:
      Encounter ID  -> fills lat/lon/ν and ET(TDB) seconds via resolve_encounter()
      UTC (ISO) and/or ET/TDB seconds (if both filled, ET wins)
      Angle: colat θ (deg) OR latitude (deg) [checkbox selects interpretation]
      Lon (deg E), ν (deg)

    Actions:
      - "Resolve ν": resolves true anomaly from whichever time is filled (prefers ET)
      - "Calculate Stress": calls on_compute(lat, lon, t, sat, Diurnal, nu)
    """

    def __init__(
        self,
        parent,
        get_satellite: Optional[Callable] = None,
        get_diurnal_cls: Optional[Callable] = None,
        resolve_encounter: Optional[Callable[[str], Optional[dict]]] = None,
        resolve_true_anomaly: Optional[Callable[[str], Optional[float]]] = None,          # UTC string -> ν
        resolve_true_anomaly_from_et: Optional[Callable[[float], Optional[float]]] = None, # ET seconds -> ν
        resolve_nu_from_time: Optional[Callable[[Union[str, float], bool], Optional[float]]] = None,  # (value, is_et)
        on_compute: Optional[Callable] = None,
        **kwargs,
    ):
        super().__init__(parent, id=wx.ID_ANY, **kwargs)

        # Hooks
        self.get_satellite = get_satellite or (lambda: None)
        self.get_diurnal_cls = get_diurnal_cls or (lambda: None)
        self.get_eccentricity = (kwargs.pop("get_eccentricity", None) or (lambda: 0.0))
        self.resolve_encounter = resolve_encounter or (lambda _enc: None)
        self.resolve_true_anomaly = resolve_true_anomaly or (lambda _utc: None)
        self.resolve_true_anomaly_from_et = resolve_true_anomaly_from_et or (lambda _et: None)
        self.resolve_nu_from_time = resolve_nu_from_time  # may be None (we’ll bridge below)
        self.on_compute = on_compute or (lambda **_k: None)

        # Where we’re currently mounted (if anywhere)
        self._inputs_host: Optional[wx.Window] = None
        self._results_host: Optional[wx.Window] = None

        # Build inner content (not mounted yet)
        self._build_inputs_inner()
        self._build_results_inner()
    from .utils import true_to_mean_anomaly_deg

    # ---------------- Build the two groups (inner content only) ----------------
    def _build_inputs_inner(self):
        """Create the inner panel for 'Point Inputs' (not mounted yet)."""
        self._in_inner = wx.Panel(self)
        col = wx.BoxSizer(wx.VERTICAL)

        # ---- Encounter row (dropdown + button) ----
        row1 = wx.BoxSizer(wx.HORIZONTAL)
        row1.Add(wx.StaticText(self._in_inner, label="Encounter ID:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.choice_enc = wx.Choice(self._in_inner, choices=["(Load encounters)"])
        self.choice_enc.Enable(False)
        self.choice_enc.Bind(wx.EVT_CHOICE, self._on_choice_enc)
        row1.Add(self.choice_enc, 0, wx.RIGHT, 8)
        self.btn_fill = wx.Button(self._in_inner, label="Load from encounters")
        self.btn_fill.Bind(wx.EVT_BUTTON, self._on_load_selected_encounter)
        row1.Add(self.btn_fill, 0)
        col.Add(row1, 0, wx.ALL, 4)

        # ---- Time inputs (always visible; ET wins if both filled) ----
        help_row = wx.BoxSizer(wx.HORIZONTAL)
        tip = wx.StaticText(self._in_inner, label="Optional: Manually enter UTC (ISO), then resolve true anomaly.")
        tip.SetForegroundColour(wx.Colour(90, 90, 90))
        help_row.Add(tip, 0)
        col.Add(help_row, 0, wx.LEFT | wx.RIGHT | wx.TOP, 4)

        row_utc = wx.BoxSizer(wx.HORIZONTAL)
        row_utc.Add(wx.StaticText(self._in_inner, label="UTC (ISO):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.txt_utc = wx.TextCtrl(self._in_inner, value="", size=(220, -1))
        self.txt_utc.SetMinSize((self.FromDIP(220), -1))
        self.txt_utc.SetHint("e.g. 2025-08-29T17:00:00Z")
        row_utc.Add(self.txt_utc, 0, wx.RIGHT, 8)
        col.Add(row_utc, 0, wx.ALL, 4)

        '''
        row_et = wx.BoxSizer(wx.HORIZONTAL)
        row_et.Add(wx.StaticText(self._in_inner, label="ET/TDB (s past J2000):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.txt_et = wx.TextCtrl(self._in_inner, value="", size=(220, -1))
        self.txt_et.SetMinSize((self.FromDIP(220), -1))
        self.txt_et.SetHint("e.g. 123456789.0")
        row_et.Add(self.txt_et, 0)
        col.Add(row_et, 0, wx.ALL, 4)
        '''

        def _mark_utc(_evt):
            self._manual_time_dirty = True
            self._manual_time_src = "utc"
            _evt.Skip()


        self.txt_utc.Bind(wx.EVT_TEXT, _mark_utc)
        #self.txt_et.Bind(wx.EVT_TEXT, _mark_et)


        # ---- Resolve ν row (button + fixed-width status) ----
        row2c = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_nu = wx.Button(self._in_inner, label="Resolve ν")
        self.btn_nu.Bind(wx.EVT_BUTTON, self.autofill_phase_from_time)
        row2c.Add(self.btn_nu, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)

        self.lbl_nu_status = wx.StaticText(self._in_inner, label="", style=getattr(wx, "ST_ELLIPSIZE_END", 0))
        self.lbl_nu_status.SetMinSize((self.FromDIP(180), -1))
        self.lbl_nu_status.SetForegroundColour(wx.Colour(90, 90, 90))
        row2c.Add(self.lbl_nu_status, 0, wx.ALIGN_CENTER_VERTICAL)
        col.Add(row2c, 0, wx.ALL, 4)

        # ---- Lat/Lon (angle as colat θ by default) ----
        theta = "\N{GREEK SMALL LETTER THETA}"
        deg   = "\N{DEGREE SIGN}"

        row3 = wx.BoxSizer(wx.HORIZONTAL)

        self.lbl_angle = wx.StaticText(self._in_inner, label=f"Lat ({deg}):")
        row3.Add(self.lbl_angle, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.txt_angle = wx.TextCtrl(self._in_inner, value="", size=(100, -1))
        row3.Add(self.txt_angle, 0, wx.RIGHT, 12)

        row3.Add(wx.StaticText(self._in_inner, label=f"Lon ({deg}E):"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.txt_lon = wx.TextCtrl(self._in_inner, value="", size=(100, -1))
        row3.Add(self.txt_lon, 0)
        col.Add(row3, 0, wx.ALL, 4)

        # ---- ν (true) and M (mean) on one row ----
        row4 = wx.FlexGridSizer(0, 4, 0, 6)   # 4 cols: label ν | ν | label M | M
        row4.AddGrowableCol(1, 1)              # let ν field stretch a bit
        row4.AddGrowableCol(3, 1)              # let M field stretch a bit

        lbl_nu = wx.StaticText(self._in_inner, label=f"True Anomaly ({deg}):")
        row4.Add(lbl_nu, 0, wx.ALIGN_CENTER_VERTICAL)

        self.txt_nu = wx.TextCtrl(self._in_inner, value="", size=(100, -1))
        row4.Add(self.txt_nu, 1, wx.EXPAND)

        lbl_M = wx.StaticText(self._in_inner, label=f"Mean Anomaly (M) ({deg}):")
        row4.Add(lbl_M, 0, wx.ALIGN_CENTER_VERTICAL)

        self.txt_M = wx.TextCtrl(self._in_inner, value="", size=(100, -1), style=wx.TE_READONLY)
        row4.Add(self.txt_M, 1, wx.EXPAND)

        # update M live when ν changes
        self.txt_nu.Bind(wx.EVT_TEXT, self._on_nu_changed)

        col.Add(row4, 0, wx.ALL | wx.EXPAND, 4)


        # ---- Action (Compute) ----
        row5 = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_compute = wx.Button(self._in_inner, label="Calculate Stress")
        self.btn_compute.Bind(wx.EVT_BUTTON, self._on_compute_clicked)
        row5.Add(self.btn_compute, 0)
        col.Add(row5, 0, wx.ALL, 6)

        self._in_inner.SetSizer(col)

    def _build_results_inner(self):
        """Create the inner panel for 'Results' (not mounted yet)."""
        self._out_inner = wx.Panel(self)
        grid = wx.FlexGridSizer(0, 6, 6, 8)
        grid.AddGrowableCol(1, 1); grid.AddGrowableCol(3, 1); grid.AddGrowableCol(5, 1)

        def R(lbl): return wx.StaticText(self._out_inner, label=lbl)

        grid.Add(R("Stt [kPa]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_ttt = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_ttt, 1, wx.EXPAND)

        grid.Add(R("Spt [kPa]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_tpt = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_tpt, 1, wx.EXPAND)

        grid.Add(R("Spp [kPa]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_tpp = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_tpp, 1, wx.EXPAND)

        grid.Add(R("σ1 [kPa]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_s1 = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_s1, 1, wx.EXPAND)

        grid.Add(R("σ3 [kPa]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_s3 = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_s3, 1, wx.EXPAND)

        grid.Add(R("α [°]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.out_alpha = wx.TextCtrl(self._out_inner, style=wx.TE_READONLY)
        grid.Add(self.out_alpha, 1, wx.EXPAND)

        self._out_inner.SetSizer(grid)

    # ---------------- Mount helpers (safe re-parent/re-sizer) ----------------
    def _detach_from_old_host(self, which: str):
        """Detach inner from previous host sizer if it was mounted."""
        if which == "inputs" and self._inputs_host:
            old = self._inputs_host
            old.SetSizer(None)
            self._inputs_host = None
        if which == "results" and self._results_host:
            old = self._results_host
            old.SetSizer(None)
            self._results_host = None

    def mount_inputs_into(self, parent_panel: wx.Window):
        """Place the 'Point Inputs' group into parent_panel inside a StaticBoxSizer."""
        self._detach_from_old_host("inputs")

        box = wx.StaticBox(parent_panel, label="Select Encounter or Input Location/Time:")
        style_staticbox(box)
        sbs = wx.StaticBoxSizer(box, wx.VERTICAL)

        # Reparent inner content so the StaticBox owns it
        if self._in_inner.GetParent() is not box:
            self._in_inner.Reparent(box)

        sbs.Add(self._in_inner, 0, wx.EXPAND | wx.ALL, self.FromDIP(6))
        parent_panel.SetSizer(sbs)
        parent_panel.Layout()
        self._inputs_host = parent_panel

    def mount_results_into(self, parent_panel: wx.Window):
        """Place the 'Results' group into parent_panel inside a StaticBoxSizer."""
        self._detach_from_old_host("results")

        box = wx.StaticBox(parent_panel, label="Computed Stresses at Specified Location:")
        style_staticbox(box)
        sbs = wx.StaticBoxSizer(box, wx.VERTICAL)

        if self._out_inner.GetParent() is not box:
            self._out_inner.Reparent(box)

        sbs.Add(self._out_inner, 0, wx.EXPAND | wx.ALL, self.FromDIP(6))
        parent_panel.SetSizer(sbs)
        parent_panel.Layout()
        self._results_host = parent_panel

    # ---------- Encounter dropdown API ----------
    def set_encounter_choices(self, ids: Iterable[str], select: Optional[str] = None):
        ids = list(ids) if ids else []
        print(f"[point] set_encounter_choices self={self}")
        print(f"[point] received ids count={len(ids)} sample={ids[:10]} select={select!r}")

        self.choice_enc.Clear()
        if ids:
            self.choice_enc.AppendItems(ids)
            self.choice_enc.Enable(True)
            if select and select in ids:
                self.choice_enc.SetStringSelection(select)
            else:
                self.choice_enc.SetSelection(0)
        else:
            self.choice_enc.Append("(Load encounters)")
            self.choice_enc.SetSelection(0)
            self.choice_enc.Enable(False)

        print(f"[point] choice count after fill={self.choice_enc.GetCount()} enabled={self.choice_enc.IsEnabled()}")

    def get_selected_encounter_id(self) -> Optional[str]:
        return self.choice_enc.GetStringSelection() if self.choice_enc.IsEnabled() else None

    # ---------- Event wiring ----------
    def _on_choice_enc(self, _evt):
        # no immediate action; click "Load from encounters"
        pass

    def _on_load_selected_encounter(self, _evt):
        enc = self.get_selected_encounter_id()
        if not enc:
            wx.MessageBox("Select an encounter first.", "Point stress", wx.OK | wx.ICON_INFORMATION)
            return
        try:
            rec = self.resolve_encounter(enc)
        except Exception as e:
            wx.MessageBox(f"Lookup failed: {e}", "Point stress", wx.OK | wx.ICON_ERROR)
            return
        if not rec:
            wx.MessageBox(f"Encounter '{enc}' not found", "Point stress", wx.OK | wx.ICON_WARNING)
            return

        # Fill available fields; spreadsheet ET is authoritative
        et = rec.get("et_tdb_sec")
        lat = rec.get("lat_deg")
        lon = rec.get("lon_deg")
        nu = rec.get("true_anom_deg")

        # Show seconds
        #self.txt_et.SetValue("" if et is None else str(et))

        # Lat/Lon
        if lat is not None:
            self.txt_angle.SetValue(f"{float(lat):.2f}")
        if lon is not None:
            self.txt_lon.SetValue(f"{float(lon):.2f}")
        if nu is not None:
            self.txt_nu.SetValue(f"{float(nu):.1f}")
            M = rec.get("mean_anom_deg") or rec.get("M_deg")
            if M is None and (nu is not None):
                try:
                    e = float(self.get_eccentricity())
                except Exception:
                    e = 0.0
                M = true_to_mean_anomaly_deg(float(nu), e)
            self.txt_M.SetValue("" if M is None else f"{float(M)%360.0:.1f}")


        # Clear UTC if we populated ET
        if et is not None:
            self.txt_utc.SetValue("")
        self._manual_time_dirty = False
        self._manual_time_src = None


    # ---------- ν resolver ----------
    def autofill_phase_from_time(self, evt=None, *, prefer_horizons_M: bool = False) -> None:
        """
        Read UTC ISO from UI, query Horizons for Europa ν (about Jupiter), compute M, update UI.
        """
        utc_iso = self.txt_utc.GetValue().strip()
        if not utc_iso:
            wx.MessageBox("Enter a UTC ISO time first (e.g., 2026-02-18T12:00:00Z).",
                          "Missing time", wx.OK | wx.ICON_WARNING)
            return

        try:
            phase = europa_phase_from_horizons(
                utc_iso,
                center="JUPITER",                # uses 500@599
                prefer_M_from_horizons=prefer_horizons_M,
                debug_payload=False,
            )

            nu_deg = float(phase.nu_deg) % 360.0
            e = float(phase.e)

            if prefer_horizons_M and (phase.M_deg is not None) and np.isfinite(phase.M_deg):
                M_deg = float(phase.M_deg) % 360.0
                src = "horizons:M"
            else:
                M_deg = float(true_to_mean_anomaly_deg(nu_deg, e)) % 360.0
                src = "horizons:nu->M"

            # Update UI fields (format however you like)
            self.txt_nu.SetValue(f"{nu_deg:.1f}")
            self.txt_M.SetValue(f"{M_deg:.1f}")
            #self.e_ctrl.SetValue(f"{e:.12g}")
            if hasattr(self, "phase_src_ctrl"):
                self.phase_src_ctrl.SetValue(src)

            # Optional: stash in panel state for downstream computation
            self.current_phase = dict(
                utc_iso=phase.utc_iso,
                nu_deg=nu_deg,
                M_deg=M_deg,
                e=e,
                src=src,
            )

            # Optional: notify parent to re-render plots / recompute stresses
            evt = wx.CommandEvent(wx.EVT_TEXT.typeId, self.GetId())
            wx.PostEvent(self, evt)

        except Exception as err:
            wx.MessageBox(f"Horizons query failed:\n\n{err!r}",
                          "Horizons error", wx.OK | wx.ICON_ERROR)
      
    '''       
    def _on_resolve_nu(self, _evt):
        sval_et  = (self.txt_et.GetValue()  or "").strip()
        sval_utc = (self.txt_utc.GetValue() or "").strip()

        # Debug helper (optional, but useful once)
        # print("DEBUG utc=", repr(sval_utc), "et=", repr(sval_et), "dirty=", getattr(self, "_manual_time_dirty", None), "src=", getattr(self, "_manual_time_src", None))

        # Decide whether we're resolving from manual inputs.
        manual_dirty = bool(getattr(self, "_manual_time_dirty", False))
        manual_src   = getattr(self, "_manual_time_src", None)

        # Manual behavior:
        # - if user typed in UTC last, use UTC even if ET still has something
        # - if user typed in ET last, use ET
        # - if user typed but src is unknown, prefer UTC if present else ET
        if manual_dirty:
            if manual_src == "utc":
                sval_et = ""  # force UTC path
            elif manual_src == "et":
                sval_utc = ""  # force ET path
            else:
                if sval_utc:
                    sval_et = ""
                elif sval_et:
                    sval_utc = ""

        # If both are still filled (non-manual or ambiguous), keep your original rule: ET wins.
        if sval_et:
            try:
                et = float(sval_et)
            except ValueError:
                wx.MessageBox("ET must be a number: seconds past J2000 (TDB).", "Point stress",
                            wx.OK | wx.ICON_WARNING)
                return
            try:
                if callable(getattr(self, "resolve_nu_from_time", None)):
                    # IMPORTANT: disable fallback when resolving from manual inputs
                    nu = self.resolve_nu_from_time(et, is_et=True, allow_fallback=not manual_dirty)
                elif callable(getattr(self, "resolve_true_anomaly_from_et", None)):
                    nu = self.resolve_true_anomaly_from_et(et)
                else:
                    raise RuntimeError("No ν resolver for ET seconds is available.")
            except Exception as e:
                wx.MessageBox(f"ν resolve failed:\n{e}", "Point stress", wx.OK | wx.ICON_ERROR)
                return

        elif sval_utc:
            try:
                if callable(getattr(self, "resolve_nu_from_time", None)):
                    nu = self.resolve_nu_from_time(sval_utc, is_et=False, allow_fallback=not manual_dirty)
                elif callable(getattr(self, "resolve_true_anomaly", None)):
                    nu = self.resolve_true_anomaly(sval_utc)
                else:
                    raise RuntimeError("No ν resolver for UTC is available.")
            except Exception as e:
                wx.MessageBox(f"ν resolve failed:\n{e}", "Point stress", wx.OK | wx.ICON_ERROR)
                return

        else:
            wx.MessageBox("Enter UTC or ET/TDB first.", "Point stress", wx.OK | wx.ICON_INFORMATION)
            return

        if nu is None:
            diag = None
            if hasattr(self, "get_last_horizons_diag"):
                diag = self.get_last_horizons_diag()
            elif hasattr(self, "_last_horizons_diag"):
                diag = self._last_horizons_diag

            msg = "Could not resolve ν for the given time."
            if isinstance(diag, dict):
                msg += "\n\nDiag:\n" + "\n".join([f"{k}: {v}" for k, v in diag.items() if v is not None])

            wx.MessageBox(msg, "Point stress", wx.OK | wx.ICON_WARNING)
            print(getattr(self, "_last_horizons_diag", None))


            wx.MessageBox("Could not resolve ν for the given time.", "Point stress", wx.OK | wx.ICON_WARNING)
            self.lbl_nu_status.SetLabel("ν not resolved")
            if hasattr(self, "txt_M"):
                self.txt_M.SetValue("")
            self.Layout()
            return

        # Set ν and update M
        try:
            nu = float(nu) % 360.0
            self.txt_nu.SetValue(f"{nu:.6f}")

            if hasattr(self, "txt_M"):
                try:
                    e = float(self.get_eccentricity())
                except Exception:
                    e = 0.0
                M = _true_to_mean_anomaly_deg(nu, e)
                self.txt_M.SetValue(f"{M:.2f}")

            self.lbl_nu_status.SetLabel("ν resolved")
            self.Layout()

        except Exception:
            self.txt_nu.SetValue(str(nu))
            if hasattr(self, "txt_M"):
                try:
                    nu_f = float(nu)
                    try:
                        e = float(self.get_eccentricity())
                    except Exception:
                        e = 0.0
                    M = _true_to_mean_anomaly_deg(nu_f, e)
                    self.txt_M.SetValue(f"{M:.2f}")
                except Exception:
                    self.txt_M.SetValue("")
            self.lbl_nu_status.SetLabel("ν resolved")
            self.Layout()

        print(self.get_last_horizons_diag())
    '''



    def _update_M_from_nu(self):
        s = self.txt_nu.GetValue().strip()
        if not s:
            self.txt_M.SetValue("")
            return
        try:
            nu = float(s)
        except ValueError:
            self.txt_M.SetValue("")
            return
        try:
            e = float(self.get_eccentricity())
        except Exception:
            e = 0.0
        M = true_to_mean_anomaly_deg(nu, e)
        self.txt_M.SetValue(f"{M:.2f}")

    def _on_nu_changed(self, _evt):
        self._update_M_from_nu()



    # ---------- Compute click ----------
    def _on_compute_clicked(self, _evt):
        def _num_or_none(ctrl):
            try:
                s = ctrl.GetValue().strip()
            except AttributeError:
                s = str(ctrl).strip()
            if not s:
                return None
            try:
                return float(s)
            except Exception:
                return None

        # -------- latitude (from angle or lat) --------
        ang = _num_or_none(self.txt_angle)
        if ang is None:
            wx.MessageBox("Angle is required.", "Point stress", wx.OK | wx.ICON_WARNING)
            return
        #is_colat = bool(self.chk_is_colat.GetValue())
        #lat = 90.0 - ang if is_colat else ang
        lat = ang

        # -------- longitude --------
        lon = _num_or_none(self.txt_lon)
        if lon is None:
            wx.MessageBox("Lon is required.", "Point stress", wx.OK | wx.ICON_WARNING)
            return

        # -------- ν and time inputs --------
        nu = _num_or_none(self.txt_nu)

        t_sec: Optional[float] = None
        utc_str: Optional[str] = None

        utc_str = self.txt_utc.GetValue().strip() or None
        t_sec = None  # Point panel no longer supports ET


        if nu is None and t_sec is None:
            wx.MessageBox("Provide ν (deg)",
                          "Point stress", wx.OK | wx.ICON_WARNING)
            return

        # Upstream handles
        try:
            sat = self.get_satellite()
        except Exception:
            sat = None
        try:
            Diurnal = self.get_diurnal_cls()
        except Exception:
            Diurnal = None

        # Compute
        try:
            result = self.on_compute(lat=lat, lon=lon, t=t_sec, sat=sat, Diurnal=Diurnal, nu=nu)
        except Exception as e:
            wx.MessageBox(f"Compute failed:\n{e}", "Point stress", wx.OK | wx.ICON_ERROR)
            return

        if result is None:
            wx.MessageBox("No result returned by compute (check ν/time inputs).",
                          "Point stress", wx.OK | wx.ICON_WARNING)
            return

        # Display results
        def _scalar(val):
            try:
                import numpy as np
                if isinstance(val, np.ndarray):
                    return float(val.reshape(())) if val.ndim == 0 else float(val.ravel()[0])
            except Exception:
                pass
            return None if val is None else float(val)

        try:
            if isinstance(result, dict):
                Ttt = result.get("Ttt") or result.get("sigma_tt") or result.get("Stt")
                Tpt = result.get("Tpt") or result.get("sigma_tp") or result.get("Spt")
                Tpp = result.get("Tpp") or result.get("sigma_pp") or result.get("Spp")
                s1 = result.get("sigma1")
                s3 = result.get("sigma3")
                a = result.get("alpha_deg") or result.get("alpha") or result.get("alpha_degrees")
                if None in (Ttt, Tpt, Tpp):
                    wx.MessageBox("Compute returned a dict but was missing Ttt/Tpt/Tpp.",
                                  "Point stress", wx.OK | wx.ICON_WARNING)
                    return
                self.set_results(_scalar(Ttt), _scalar(Tpt), _scalar(Tpp),
                                 _scalar(s1) if s1 is not None else None,
                                 _scalar(s3) if s3 is not None else None,
                                 _scalar(a) if a is not None else None)
                return

            if isinstance(result, (tuple, list)) and len(result) >= 3:
                vals = [_scalar(v) if v is not None else None for v in result[:6]]
                while len(vals) < 6:
                    vals.append(None)
                self.set_results(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5])
                return

            wx.MessageBox("Compute returned an unexpected format.",
                          "Point stress", wx.OK | wx.ICON_WARNING)
        except Exception as e:
            wx.MessageBox(f"Display failed:\n{e}", "Point stress", wx.OK | wx.ICON_ERROR)

    # ---------- Results helper ----------
    def set_results(
        self,
        Ttt_pa: float,
        Tpt_pa: float,
        Tpp_pa: float,
        sigma1_pa: Optional[float] = None,
        sigma3_pa: Optional[float] = None,
        alpha_deg: Optional[float] = None,
    ):
        Ttt = float(Ttt_pa)
        Tpt = float(Tpt_pa)
        Tpp = float(Tpp_pa)

        if sigma1_pa is None or sigma3_pa is None:
            m = 0.5 * (Ttt + Tpp)
            r = float(np.hypot(0.5 * (Ttt - Tpp), Tpt))
            s1 = m + r
            s3 = m - r
        else:
            s1 = float(sigma1_pa)
            s3 = float(sigma3_pa)

        if alpha_deg is None:
            alpha_deg = 0.5 * np.degrees(np.arctan2(2.0 * Tpt, (Ttt - Tpp)))

        self.out_ttt.SetValue(f"{Ttt * _KPA: .6g}")
        self.out_tpt.SetValue(f"{Tpt * _KPA: .6g}")
        self.out_tpp.SetValue(f"{Tpp * _KPA: .6g}")
        self.out_s1.SetValue(f"{s1 * _KPA: .6g}")
        self.out_s3.SetValue(f"{s3 * _KPA: .6g}")
        self.out_alpha.SetValue(f"{float(alpha_deg): .6g}")

























