# stressviz/scalar_plot_panel.py
import wx, wx.adv, math, os, re
import datetime as _dt
import importlib
import numpy as np
import traceback
import matplotlib.patheffects as pe
from typing import Callable, List, Tuple, Optional, Dict, Any
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import transforms as mtransforms
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas, NavigationToolbar2WxAgg as Toolbar
from .orbital_panel import OrbitalPositionPanel
from .encounters_io import load_plume_observations_txt

__all__ = ["ScalarPlotPanel"]

_KPA = 1e-3  # Pa -> kPa

from .utils import resolve_stresscalc
from .utils import true_to_mean_anomaly_deg as _true_to_mean_anomaly_deg

STRESS_PLOT_ID = "__STRESS_PLOT__"
STRESS_PLOT_LABEL = "Stress Plot"


def _nu_to_time_seconds(nu_deg, diurn):
    """Map ν (deg) -> time (sec) using diurn.omega (rad/s) or close cousin."""
    if nu_deg is None:
        return 0.0
    for key in ("omega", "mean_motion", "n"):
        if hasattr(diurn, key):
            val = float(getattr(diurn, key))
            if val != 0.0:
                return float(np.deg2rad(nu_deg)) / val
    if hasattr(diurn, "deg_per_sec"):
        dps = float(getattr(diurn, "deg_per_sec")) or 1.0
        return float(nu_deg) / dps
    return 0.0


def _tensor_grid(calc, theta_2d, phi_2d, t_sec):
    """Get (Ttt,Tpt,Tpp) on a grid; broadcast if supported, else loop."""
    try:
        Ttt, Tpt, Tpp = calc.tensor(theta_2d, phi_2d, t_sec)
        return np.asarray(Ttt, float), np.asarray(Tpt, float), np.asarray(Tpp, float)
    except Exception:
        ny, nx = theta_2d.shape
        Ttt = np.empty((ny, nx), float)
        Tpt = np.empty((ny, nx), float)
        Tpp = np.empty((ny, nx), float)
        for i in range(ny):
            for j in range(nx):
                a, b, c = calc.tensor(float(theta_2d[i, j]), float(phi_2d[i, j]), float(t_sec))
                Ttt[i, j] = float(a)
                Tpt[i, j] = float(b)
                Tpp[i, j] = float(c)
        return Ttt, Tpt, Tpp
    
def _mean_to_true_anomaly_deg(M_deg: float, e: float) -> float:
    """
    Invert Kepler's equation to get true anomaly ν (deg) from mean anomaly M (deg) and eccentricity e.
    Steps:
      - Solve M = E - e sin E for eccentric anomaly E via Newton's method
      - Convert to ν with tan(ν/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    Robust for small e typical of Europa.
    """
    M = np.deg2rad(float(M_deg)) % (2*np.pi)
    e = float(e)

    # Initial guess for E
    E = M if e < 0.8 else np.pi

    # Newton iterations
    for _ in range(12):
        f  = E - e*np.sin(E) - M
        fp = 1.0 - e*np.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < 1e-12:
            break

    # E -> ν
    tan_half_nu = np.sqrt((1.0 + e) / (1.0 - e)) * np.tan(0.5 * E)
    nu = 2.0 * np.arctan(tan_half_nu)
    # normalize to [0, 2π)
    if nu < 0:
        nu += 2*np.pi * np.ceil(-nu/(2*np.pi))
    return (np.rad2deg(nu) % 360.0)



class ScalarPlotPanel(wx.Panel):
    """
    Native StressViz scalar field plotter with SatStress-like controls.

    Usage:
      panel.set_axes(lats_deg, lons_deg)
      panel.bind_orbit_series(
          nus_deg=[...],                           # ν samples, internal
          eval_fn=lambda M_deg: (Ttt_pa, Tpt_pa, Tpp_pa),   # <-- now expects MEAN anomaly (deg)
          initial_nu_deg=...,
      )

    You can lock a fixed color range (kPa) at any time:
      panel.set_fixed_stress_range_kpa(-100, 100)
      # or revert:
      panel.clear_fixed_stress_range()
    """
    def __init__(
        self, 
        parent, 
        get_system_id: Optional[Callable[[], str]] = None, 
        get_encounter_id: Optional[Callable[[], Optional[str]]] = None,
        get_eccentricity: Optional[Callable[[], float]] = None,   # NEW: supply e for ν->M
    ):
        super().__init__(parent)
        self._lats = None        # (Ny,)
        self._lons = None        # (Nx,)
        self._T = None           # last (Ttt, Tpt, Tpp) in Pa
        self._Z_custom = None    # optional custom scalar in kPa
        self._field = "tens"     # tens|comp|mean|diff|custom
        self._east_positive = True

        self._ax_pos = None  # axes rectangle (figure coords) captured after first draw

        self._get_system_id = get_system_id
        self._get_encounter_id = get_encounter_id
        self._get_eccentricity = get_eccentricity or (lambda: 0.0)  # safe default

        # Scale state (kPa)
        self._auto_range = True
        self._vmin = -100.0
        self._vmax = 100.0
        self._fixed_clim_kpa: Optional[Tuple[float, float]] = None

        # Vector sampling
        self._vec_nlat = 10
        self._vec_nlon = 10
        self._principal_artists = []
        self._vec_max_deg = 10.0   # longest half-length ~10°

        # Vector scaling: 100 kPa → this many degrees of total line length
        self._vec_ref_kpa = 100.0
        self._vec_len_for_ref_deg = 18.0   # total length at 100 kPa (so half-length = 4°)

        # Vector magnitude scaling (shared with the reference bar)
        self._ref_kpa = 100.0
        self._vec_scale_deg_per_100kpa = 20.0
        self._legend_art = []

        # Encounters UI/state
        self._encounters = []  # list of dicts: {"name": str, "nu_deg": float|None, "utc_iso": str|None}
        self._nu_resolver = None  # callable(utc_iso:str) -> float (deg), optional
        self._add_marker = None
        self._remove_marker = None
        self._goto_nu = None     # kept for back-compat; we pass M through it
        self._plume_M_by_encounter: dict[tuple[str, str], float] = {}  # (enc_id, obs_id) -> M_deg
        self._plumes_loaded = False 
        self._plume_meta_by_id: dict[str, dict] = {}
        self._juice_loaded = False
        self._juice_meta_by_id: dict[str, dict] = {}


        # Matplotlib handles
        self._im = None
        self._cbar = None
        self._nus = None          # list of ν values (deg)
        self._eval_fn = None      # Callable[[M_deg], (Ttt,Tpt,Tpp) in Pa]
        self._idx = 0

        self._init_snap_done = False          # once we snap at startup, don't re-override
        self._startup_target_M = None         # target M for the initial snap


        self._build_ui()

        # --- Return-to-Encounter gating state ---
        self._current_M_deg: float | None = None
        self._selected_encounter_id: str | None = None
        self._selected_encounter_M_deg: float | None = None

        # tolerance for “same M”
        self._M_TOL_DEG = 1e-2  # 0.01°

        self._enc_meta_by_id: dict[str, dict] = {}   # {enc_id: {'utc_ca': aware_dt, 'M_ca_deg': float, 'period_hours': float|None}}
        self._selected_encounter_ca_utc = None       # aware UTC datetime for active encounter


    # ---- small helper for eccentricity ----
    def _ecc_safe(self) -> float:
        try:
            return float(self._get_eccentricity()) if callable(self._get_eccentricity) else 0.0
        except Exception:
            return 0.0

    # ---------------- UI ----------------
    def _build_ui(self):
        outer = wx.BoxSizer(wx.VERTICAL)

        # Header (now references M)
        outer.Add(
            wx.StaticText(
                self,
                label="Tension positive (kPa). Use M slider to scrub orbit; Save series exports one PNG per step.",
            ),
            0, wx.ALL | wx.EXPAND, 6,
        )

        main = wx.BoxSizer(wx.HORIZONTAL)

        # =========================
        # LEFT: map + slider + orbit plot
        # =========================
        left_panel = wx.Panel(self)
        left = wx.BoxSizer(wx.VERTICAL)
        left_panel.SetSizer(left)

        # --- Matplotlib figure/canvas ---
        self.fig = Figure(figsize=(6, 5), constrained_layout=False)
        self.ax = self.fig.add_subplot(111)
        self.figure = self.fig
        self.canvas = FigureCanvas(left_panel, -1, self.fig)
        self.canvas.SetMinSize((self.FromDIP(440), self.FromDIP(440)))

        # Keep a comfortable axes rectangle (l, b, w, h) in fig fractions
        self._ax_rect = [0.10, 0.16, 0.72, 0.70]
        self.ax.set_position(self._ax_rect)

        # Toolbar under the canvas
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        # Persistent colorbar axis tied to main map axis
        self._cbar_size = "4%"   # width of the colorbar
        self._cbar_pad  = "2%"   # GAP between map and colorbar
        divider = make_axes_locatable(self.ax)
        self._cax = divider.append_axes("right", size=self._cbar_size, pad=self._cbar_pad)

        left.Add(self.canvas, 3, wx.EXPAND)
        left.Add(self.toolbar, 0, wx.EXPAND)

        # ----- Keep a strict 2:1 x:y data ratio while staying interactive -----
        self._target_aspect = 0.5  # aspect = (y/x); 0.5 => x twice y
        self.ax.set_aspect(self._target_aspect, adjustable="datalim", anchor="C")
        try:
            self.ax.set_box_aspect(0.5)
        except Exception:
            pass

        self._aspect_guard = False

        def _enforce_2to1(ax):
            if self._aspect_guard:
                return
            self._aspect_guard = True
            try:
                x0, x1 = ax.get_xlim()
                y0, y1 = ax.get_ylim()
                xspan = abs(x1 - x0)
                yspan = abs(y1 - y0)
                desired_yspan = xspan * self._target_aspect
                if abs(yspan - desired_yspan) > 1e-12:
                    cy = 0.5 * (y0 + y1)
                    new_y0 = cy - desired_yspan / 2.0
                    new_y1 = cy + desired_yspan / 2.0
                    if y1 < y0:
                        new_y0, new_y1 = new_y1, new_y0
                    ax.set_ylim(new_y0, new_y1)
            finally:
                self._aspect_guard = False

        self.ax.callbacks.connect("xlim_changed", _enforce_2to1)
        self.ax.callbacks.connect("ylim_changed", _enforce_2to1)
        self.ax.autoscale(enable=True, axis="both", tight=False)

        # --- Orbit controls: slider row + save below (aligned under slider) ---
        ctl = wx.FlexGridSizer(rows=2, cols=3, vgap=3, hgap=6)
        ctl.AddGrowableCol(1, 1)  # slider column stretches

        self.btn_prev = wx.Button(left_panel, label="◀")
        self.btn_next = wx.Button(left_panel, label="▶")
        self.sld_orbit = wx.Slider(left_panel, value=0, minValue=0, maxValue=0, style=wx.SL_HORIZONTAL)
        self.sld_orbit.SetMinSize((self.FromDIP(320), -1))
        self.lbl_nu = wx.StaticText(left_panel, label="Orbital Position: --.--°")
        self.lbl_nu.SetMinSize((self.FromDIP(88), -1))
        self.btn_save_series = wx.Button(left_panel, label="Save Stress Series")
        self.btn_save_snap = wx.Button(left_panel, label="Save Current Plots")
        self.btn_save_orbit = wx.Button(left_panel, label="Save Orbit Plot")

        # Wire up
        self.btn_prev.Bind(wx.EVT_BUTTON, lambda e: self._nudge(-1))
        self.btn_next.Bind(wx.EVT_BUTTON, lambda e: self._nudge(+1))
        self.sld_orbit.Bind(wx.EVT_SLIDER, self._on_slide)
        self.btn_save_series.Bind(wx.EVT_BUTTON, self._on_save_series)
        self.btn_save_orbit.Bind(wx.EVT_BUTTON, self._on_save_orbit)

        # Row 0 layout
        ctl.Add(self.btn_prev, 0, wx.ALIGN_CENTER_VERTICAL)
        ctl.Add(self.sld_orbit, 0, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        rpack = wx.BoxSizer(wx.HORIZONTAL)
        rpack.Add(self.btn_next, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        rpack.Add(self.lbl_nu, 0, wx.ALIGN_CENTER_VERTICAL)
        ctl.Add(rpack, 0, wx.ALIGN_CENTER_VERTICAL)

        # Row 1
        ctl.AddSpacer(1)

        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        btn_row.Add(self.btn_save_series, 0, wx.TOP, 4)
        btn_row.AddSpacer(8)
        btn_row.Add(self.btn_save_orbit, 0, wx.TOP, 4)

        ctl.Add(btn_row, 0, wx.ALIGN_LEFT)

        ctl.AddSpacer(1)

        left.Add(ctl, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 4)

        # --- Time nudges row (± minutes) ---
        trow = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_tminus = wx.Button(left_panel, label="−1h", style=wx.BU_EXACTFIT)
        self.btn_tplus  = wx.Button(left_panel, label="+1h", style=wx.BU_EXACTFIT)
        self.btn_tminus.Bind(wx.EVT_BUTTON, lambda e: self._nudge_time_minutes(-60))
        self.btn_tplus.Bind(wx.EVT_BUTTON, lambda e: self._nudge_time_minutes(+60))
        trow.Add(self.btn_tminus, 0, wx.RIGHT, 6)
        trow.Add(self.btn_tplus,  0)
        left.Add(trow, 0, wx.ALIGN_LEFT | wx.LEFT | wx.RIGHT | wx.BOTTOM, 4)

        # Orbital position mini-plot + side legend
        self._orbit_box = wx.StaticBox(left_panel, label="Orbital position (M)")
        self._orbit_sizer = wx.StaticBoxSizer(self._orbit_box, wx.VERTICAL)

        orbit_row = wx.BoxSizer(wx.HORIZONTAL)

        # --- orbit plot panel (left) ---
        self.orb_panel = OrbitalPositionPanel(left_panel)
        self.orb_panel._use_side_legend = True
        ORBIT_H = self.FromDIP(250)
        self.orb_panel.SetMinSize((self.FromDIP(360), ORBIT_H))
        self.orb_panel.SetMaxSize((-1, ORBIT_H))

        orbit_row.Add(self.orb_panel, 3, wx.EXPAND | wx.ALL, 4)

        # --- legend panel (right) ---
        self.orbit_legend_panel = wx.Panel(left_panel)
        self.orbit_legend_panel.SetMinSize((self.FromDIP(190), ORBIT_H))

        legend_outer = wx.StaticBoxSizer(
            wx.StaticBox(self.orbit_legend_panel, label="Legend"),
            wx.VERTICAL,
        )

        # This is the grey box that scrolls
        self.orbit_legend_scroll = wx.ScrolledWindow(
            self.orbit_legend_panel,
            style=wx.VSCROLL | wx.BORDER_SIMPLE,
        )
        self.orbit_legend_scroll.SetScrollRate(0, self.FromDIP(12))
        self.orbit_legend_scroll.SetMinSize((-1, self.FromDIP(140)))

        # Put legend content inside the scroll window
        self.orbit_legend_sizer = wx.BoxSizer(wx.VERTICAL)
        self.orbit_legend_scroll.SetSizer(self.orbit_legend_sizer)

        legend_outer.Add(self.orbit_legend_scroll, 1, wx.EXPAND | wx.ALL, 6)
        self.orbit_legend_panel.SetSizer(legend_outer)

        orbit_row.Add(self.orbit_legend_panel, 0, wx.EXPAND | wx.TOP | wx.RIGHT | wx.BOTTOM, 4)

        self._orbit_sizer.Add(orbit_row, 1, wx.EXPAND | wx.ALL, 0)
        left.Add(self._orbit_sizer, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)

        # Add left column to main
        main.Add(left_panel, 3, wx.EXPAND | wx.ALL, 4)

        # =========================
        # RIGHT: controls
        # =========================
        right_panel = wx.Panel(self)
        right = wx.BoxSizer(wx.VERTICAL)
        right_panel.SetSizer(right)
        rp = right_panel

        # Important sizing constraints for the whole right column
        right_panel.SetMinSize((self.FromDIP(320), -1))   # a sane minimum for your controls


        grid = wx.FlexGridSizer(0, 2, 6, 8)
        grid.AddGrowableCol(1, 1)

        grid.Add(wx.StaticText(rp, label="Longitude"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.rb_dir = wx.RadioBox(
            rp,
            choices=["East Positive", "West Positive"],
            majorDimension=1,
            style=wx.RA_SPECIFY_ROWS,
        )
        self.rb_dir.SetSelection(0)
        self.rb_dir.Bind(wx.EVT_RADIOBOX, self._on_dir_change)
        grid.Add(self.rb_dir, 0, wx.EXPAND)

        grid.Add(wx.StaticText(rp, label="Plot gradient"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.choice_field = wx.Choice(rp, choices=["σ1", "σ3", "(σ1+σ3)/2", "σ1-σ3"])
        self.choice_field.SetSelection(0)
        self.choice_field.Bind(wx.EVT_CHOICE, self._on_field_change)
        grid.Add(self.choice_field, 0, wx.EXPAND)

        grid.Add(wx.StaticText(rp, label="Stress range (kPa)"), 0, wx.ALIGN_CENTER_VERTICAL)

        hr = wx.BoxSizer(wx.HORIZONTAL)

        self.sp_l = wx.TextCtrl(rp, value="-100", style=wx.TE_PROCESS_ENTER)
        self.sp_u = wx.TextCtrl(rp, value="100", style=wx.TE_PROCESS_ENTER)

        self.sp_l.Bind(wx.EVT_TEXT_ENTER, self._on_bounds_change)
        self.sp_u.Bind(wx.EVT_TEXT_ENTER, self._on_bounds_change)

        hr.Add(wx.StaticText(rp, label="L:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        hr.Add(self.sp_l, 1, wx.RIGHT, 10)
        hr.Add(wx.StaticText(rp, label="U:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        hr.Add(self.sp_u, 1)

        grid.Add(hr, 0, wx.EXPAND)

        grid.Add(wx.StaticText(rp, label="Principal lines"), 0, wx.ALIGN_CENTER_VERTICAL)
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.chk_principal = wx.CheckBox(rp, label="Show principal lines"); self.chk_principal.SetValue(True)
        self.chk_principal.Bind(wx.EVT_CHECKBOX, lambda e: self._replot()); vbox.Add(self.chk_principal, 0, wx.BOTTOM, 2)
        self.chk_s1 = wx.CheckBox(rp, label="σ1"); self.chk_s1.SetValue(True)
        self.chk_s1.Bind(wx.EVT_CHECKBOX, lambda e: self._replot()); vbox.Add(self.chk_s1, 0, wx.BOTTOM, 2)
        self.chk_s3 = wx.CheckBox(rp, label="σ3"); self.chk_s3.SetValue(True)
        self.chk_s3.Bind(wx.EVT_CHECKBOX, lambda e: self._replot()); vbox.Add(self.chk_s3, 0, wx.BOTTOM, 2)
        self.chk_scalevec = wx.CheckBox(rp, label="Scale by |σ| magnitude"); self.chk_scalevec.SetValue(True)
        self.chk_scalevec.Bind(wx.EVT_CHECKBOX, lambda e: self._replot()); vbox.Add(self.chk_scalevec, 0)
        grid.Add(vbox, 0, wx.EXPAND)

        grid.Add(wx.StaticText(rp, label="Vector density"), 0, wx.ALIGN_CENTER_VERTICAL)
        dens = wx.BoxSizer(wx.HORIZONTAL)
        self.sp_vlat = wx.SpinCtrl(rp, min=2, max=181, initial=10)
        self.sp_vlon = wx.SpinCtrl(rp, min=2, max=361, initial=10)
        self.vec_n_lat = 10; self.vec_n_lon = 10

        def _on_vec_density(_e):
            self.vec_n_lat = int(self.sp_vlat.GetValue())
            self.vec_n_lon = int(self.sp_vlon.GetValue())
            self._replot()

        self.sp_vlat.Bind(wx.EVT_SPINCTRL, _on_vec_density)
        self.sp_vlon.Bind(wx.EVT_SPINCTRL, _on_vec_density)
        dens.Add(wx.StaticText(rp, label="Lat:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        dens.Add(self.sp_vlat, 1, wx.RIGHT, 8)
        dens.Add(wx.StaticText(rp, label="Lon:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        dens.Add(self.sp_vlon, 1)
        grid.Add(dens, 0, wx.EXPAND)

        grid.Add(wx.StaticText(rp, label="Cursor (deg / kPa)"), 0, wx.ALIGN_CENTER_VERTICAL)
        rd = wx.FlexGridSizer(0, 2, 4, 6)
        rd.AddGrowableCol(1, 1)
        rd.Add(wx.StaticText(rp, label="Lat:"))
        self.txt_lat = wx.TextCtrl(rp, style=wx.TE_READONLY)
        rd.Add(self.txt_lat, 1, wx.EXPAND)
        rd.Add(wx.StaticText(rp, label="Lon:"))
        self.txt_lon = wx.TextCtrl(rp, style=wx.TE_READONLY)
        rd.Add(self.txt_lon, 1, wx.EXPAND)
        rd.Add(wx.StaticText(rp, label="Stress [kPa]:"))
        self.txt_val = wx.TextCtrl(rp, style=wx.TE_READONLY)
        rd.Add(self.txt_val, 1, wx.EXPAND)
        grid.Add(rd, 0, wx.EXPAND)

        right.Add(grid, 0, wx.EXPAND | wx.ALL, 6)

        # ----- Encounters controls -----
        enc_box = wx.StaticBox(rp, label="Encounters")
        enc_outer = wx.StaticBoxSizer(enc_box, wx.VERTICAL)

        # row: add + plot
        enc_row1 = wx.BoxSizer(wx.HORIZONTAL)

        self.btn_add_enc = wx.Button(rp, label="Select encounters")
        self.btn_nearby = wx.Button(rp, label="Show nearby events (+/-10°)")
        self.btn_move_stress = wx.Button(rp, label="Move Stress Plot")
        self.btn_plot_enc = wx.Button(rp, label="Plot on Orbit")
        #self.btn_clear_enc = wx.Button(rp, label="Clear")  # optional but nice

        enc_row1.Add(self.btn_add_enc, 0, wx.RIGHT, 6)
        #enc_row1.AddStretchSpacer(1)
        #enc_row.Add(self.btn_clear_enc, 0, wx.RIGHT, 6)
        enc_row1.Add(self.btn_nearby, 0)

        enc_outer.Add(enc_row1, 0, wx.EXPAND | wx.ALL, 4)

        enc_row2 = wx.BoxSizer(wx.HORIZONTAL)

        enc_row2.Add(self.btn_move_stress, 0)
        #enc_row2.AddStretchSpacer(1)
        enc_row2.Add(self.btn_plot_enc, 0)

        enc_outer.Add(enc_row2, 0, wx.EXPAND | wx.ALL, 4)

        # rows: chosen encounters to plot (checked = enabled) + per-encounter color picker
        if not hasattr(self, "_selected_encounter_ids"):
            self._selected_encounter_ids = []     # list[str] in user-chosen order
        if not hasattr(self, "_selected_enabled_by_id"):
            self._selected_enabled_by_id = {}     # dict[str, bool]
        if not hasattr(self, "_selected_color_by_id"):
            self._selected_color_by_id = {}       # dict[str, "#RRGGBB"]

        self.enc_rows = self._build_selected_encounter_rows(rp)
        self.enc_rows.SetMinSize((-1, self.FromDIP(170)))
        enc_outer.Add(self.enc_rows, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)
        self._ensure_multi_enc_state()
        self._sync_rows_from_model()

        # wire (no checklist bind anymore)
        self.btn_add_enc.Bind(wx.EVT_BUTTON, self._on_choose_multiple_encounters)
        self.btn_plot_enc.Bind(wx.EVT_BUTTON, self._on_plot_selected_encounters)
        self.btn_nearby.Bind(wx.EVT_BUTTON, self._on_show_nearby_events)
        self.btn_move_stress.Bind(wx.EVT_BUTTON, self._on_move_stress_plot_to_selected_encounter)

        right.Add(enc_outer, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 6)


        # Add right column to main
        main.Add(right_panel, 2, wx.EXPAND | wx.ALL, 4)

        # finalize
        outer.Add(main, 1, wx.EXPAND, 0)
        self.SetSizer(outer)

        # mouse readout + axis labels
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.ax.set_xlabel("Longitude [°E]")
        self.ax.set_ylabel("Latitude [°]")


    

    # ---------------- public API ----------------
    @staticmethod
    def _slug(s: str) -> str:
        s = re.sub(r"\s+", "_", str(s).strip())
        return re.sub(r"[^A-Za-z0-9._-]", "", s) or "unnamed"

    def _proposed_folder_name(self) -> str:
        sys_id = ""
        if callable(getattr(self, "_get_system_id", None)):
            try:
                sys_id = (self._get_system_id() or "").strip()
            except Exception:
                pass

        enc_id = ""
        if callable(getattr(self, "_get_encounter_id", None)):
            try:
                enc_id = (self._get_encounter_id() or "").strip()
            except Exception:
                pass

        sys_id = self._slug(sys_id or "System")
        enc_id = self._slug(enc_id) if enc_id else ""
        return f"{sys_id}_{enc_id}" if enc_id else sys_id

    def set_axes(self, lats_deg: np.ndarray, lons_deg: np.ndarray):
        self._lats = np.asarray(lats_deg, float)
        self._lons = np.asarray(lons_deg, float)

    def bind_orbit_series(
        self,
        nus_deg: List[float],
        eval_fn: Callable[[float], Tuple[np.ndarray, np.ndarray, np.ndarray]],
        initial_nu_deg: Optional[float] = None,
        title: Optional[str] = None,
        initial_M_deg: Optional[float] = None,   
    ):

        """nus_deg: list of ν values (deg). eval_fn must accept MEAN anomaly M (deg)."""
        self._Z_custom = None
        self._field = "tens" if self._field not in ("tens", "comp", "mean", "diff", "custom") else self._field
        self.choice_field.SetSelection(["tens", "comp", "mean", "diff", "custom"].index(self._field))

        self._nus = list(map(float, nus_deg))
        self._eval_fn = eval_fn

        # If we've already snapped to a startup M, DON'T override the index
        if getattr(self, "_init_snap_done", False) and (self._startup_target_M is not None):
            # Keep whatever _idx the snap chose; only adjust slider range/value to match
            self.sld_orbit.SetRange(0, max(0, len(self._nus) - 1))
            try:
                self.sld_orbit.SetValue(int(self._idx))
            except Exception:
                pass
        else:
            if initial_M_deg is not None:
                self._idx = self._closest_index_for_M(float(initial_M_deg) % 360.0)
            elif initial_nu_deg is not None:
                self._idx = int(np.argmin(np.abs(np.array(self._nus) - float(initial_nu_deg))))
            else:
                self._idx = 0
            self.sld_orbit.SetRange(0, max(0, len(self._nus) - 1))
            self.sld_orbit.SetValue(self._idx)

        if title:
            self.ax.set_title(title)
        self._evaluate_and_plot()


        try:
            idx0 = int(self._idx)
            nu0  = float(self._nus[idx0]) % 360.0
            M0   = float(_true_to_mean_anomaly_deg(nu0, self._ecc_safe())) % 360.0
            self._update_orbit_viz(M0)       # mini-orbit marker = M
            self.lbl_nu.SetLabel(f"M: {M0:.2f}°")  # readout = M
            self.canvas.draw_idle()
        except Exception:
            pass

    
    def _ang_diff_deg(a: float, b: float) -> float:
        """Smallest absolute difference on a circle (deg)."""
        d = (a - b + 180.0) % 360.0 - 180.0
        return abs(d)

    def _update_return_button_state(self):
        """Enable/disable 'Return to Encounter' based on current M vs selected encounter M."""
        btn = getattr(self, "btn_return_to_enc", None)
        if not btn:
            return
        curM = self._current_M_deg
        encM = self._selected_encounter_M_deg
        if (curM is not None) and (encM is not None):
            same = self._ang_diff_deg(curM, encM) <= self._M_TOL_DEG
            btn.Enable(not same)
        else:
            btn.Enable(False)

    def _ensure_plumes_loaded_into_dropdown(self):
        """
        Load plume observations into self._plume_meta_by_id.

        Stores the FULL row dict (including precomputed phase columns like
        true_anom_deg / mean_anom_deg / eccentricity / phase_src), so the scalar
        panel can read M without querying Horizons.
        """
        if not hasattr(self, "_plume_meta_by_id"):
            self._plume_meta_by_id = {}

        # If already loaded and non-empty, don't reload
        if getattr(self, "_plumes_loaded", False) and self._plume_meta_by_id:
            return

        try:
            from . import encounters_io
            default_path = os.path.join(encounters_io.DATA_DIR, "plume_observations.txt")
            print(f"[ScalarPanel] plume default path would be: {default_path}")

            from .encounters_io import load_plume_observations_txt
            rows = load_plume_observations_txt()
            print(f"[ScalarPanel] load_plume_observations_txt returned {len(rows)} rows")

            if rows:
                r0 = rows[0]
                print(f"[ScalarPanel] first plume row keys: {list(r0.keys())}")
                print(f"[ScalarPanel] first plume row mean_anom_deg={r0.get('mean_anom_deg')!r} true_anom_deg={r0.get('true_anom_deg')!r}")

        except Exception as e:
            print(f"[ScalarPanel] plume load failed (will retry later): {repr(e)}")
            self._plumes_loaded = False
            self._plume_meta_by_id = {}
            return

        def _clean(x):
            if x is None:
                return None
            s = str(x).strip()
            if not s or s.lower() in ("nan", "none", "null"):
                return None
            return s

        def _to_float(x):
            s = _clean(x)
            if s is None:
                return None
            try:
                return float(s)
            except Exception:
                return None

        plume_meta: dict[str, dict] = {}
        for r in rows or []:
            if not isinstance(r, dict):
                continue

            enc_id = _clean(r.get("encounter_id")) or _clean(r.get("id"))
            if not enc_id:
                continue
            enc_id = str(enc_id)

            # Start with a shallow copy of the entire row so we keep all columns
            meta = dict(r)

            # Normalize the fields we rely on
            meta["kind"] = "plume"
            meta["encounter_tag"] = "plume"
            meta["encounter_id"] = enc_id
            meta["utc_iso"] = _clean(meta.get("utc_iso") or meta.get("utc"))
            meta["observer"] = _clean(meta.get("observer")) or "UNKNOWN"

            # Normalize phase values (key part!)
            nu = _to_float(meta.get("true_anom_deg") or meta.get("nu_deg"))
            M  = _to_float(meta.get("mean_anom_deg") or meta.get("M_deg"))
            e  = _to_float(meta.get("eccentricity"))

            if nu is not None:
                nu = nu % 360.0
                meta["true_anom_deg"] = nu
                meta["nu_deg"] = nu
            else:
                meta["true_anom_deg"] = None
                meta["nu_deg"] = None

            if M is not None:
                M = M % 360.0
                meta["mean_anom_deg"] = M
                meta["M_deg"] = M
            else:
                meta["mean_anom_deg"] = None
                meta["M_deg"] = None

            meta["eccentricity"] = e
            meta["exposure_s"] = _to_float(meta.get("exposure_s"))
            meta["phase_src"] = _clean(meta.get("phase_src"))

            plume_meta[enc_id] = meta

        self._plume_meta_by_id = plume_meta
        self._plumes_loaded = True
        print(f"[ScalarPanel] plume_meta_by_id size={len(self._plume_meta_by_id)}")


    def _ensure_juice_loaded_into_dropdown(self):
        """
        Load JUICE Europa flybys into self._juice_meta_by_id.

        Mirrors _ensure_plumes_loaded_into_dropdown:
        - Stores FULL row dict
        - Normalizes: encounter_id, utc_iso, observer, true_anom_deg/nu_deg, mean_anom_deg/M_deg, eccentricity, phase_src
        - Converts true anomaly -> mean anomaly using the SAME converter you already use for Clipper (best-effort resolver)
        """
        import os

        if not hasattr(self, "_juice_meta_by_id"):
            self._juice_meta_by_id = {}

        # If already loaded and non-empty, don't reload
        if getattr(self, "_juice_loaded", False) and self._juice_meta_by_id:
            return

        try:
            from . import encounters_io
            default_path = os.path.join(encounters_io.DATA_DIR, "JUICE_Europa_flybys.txt")
            print(f"[ScalarPanel] JUICE default path would be: {default_path}")

            from .encounters_io import load_juice_europa_flybys_txt
            rows = load_juice_europa_flybys_txt()
            print(f"[ScalarPanel] load_juice_europa_flybys_txt returned {len(rows)} rows")

            if rows:
                r0 = rows[0]
                print(f"[ScalarPanel] first JUICE row keys: {list(r0.keys())}")
                print(f"[ScalarPanel] first JUICE row true_anom_deg={r0.get('true_anom_deg')!r} encounter={r0.get('encounter')!r}")

        except Exception as e:
            print(f"[ScalarPanel] JUICE load failed (will retry later): {repr(e)}")
            self._juice_loaded = False
            self._juice_meta_by_id = {}
            return

        def _clean(x):
            if x is None:
                return None
            s = str(x).strip()
            if not s or s.lower() in ("nan", "none", "null"):
                return None
            return s

        def _to_float(x):
            s = _clean(x)
            if s is None:
                return None
            try:
                return float(s)
            except Exception:
                return None

        # ----------------- resolve the same ν->M converter you already use for Clipper -----------------
        # We try several common patterns WITHOUT breaking if they don't exist.
        _clipper_true_to_mean = None

        # (A) method on self
        for name in (
            "true_anom_deg_to_mean_anom_deg",
            "true_to_mean_anom_deg",
            "nu_to_M_deg",
            "nu_deg_to_M_deg",
            "_true_to_mean_anom_deg",
        ):
            fn = getattr(self, name, None)
            if callable(fn):
                _clipper_true_to_mean = fn
                break

        # (B) method on an orbital/phase helper object
        if _clipper_true_to_mean is None:
            for obj_name in ("_orbital_model", "orbital_model", "_phase_model", "phase_model", "_phase_helper", "phase_helper"):
                obj = getattr(self, obj_name, None)
                if obj is None:
                    continue
                for name in (
                    "true_anom_deg_to_mean_anom_deg",
                    "true_to_mean_anom_deg",
                    "nu_to_M_deg",
                    "nu_deg_to_M_deg",
                    "_true_to_mean_anom_deg",
                ):
                    fn = getattr(obj, name, None)
                    if callable(fn):
                        _clipper_true_to_mean = fn
                        break
                if _clipper_true_to_mean is not None:
                    break

        # (C) function in encounters_io (common place for clipper utilities)
        if _clipper_true_to_mean is None:
            try:
                from .encounters_io import true_anom_deg_to_mean_anom_deg as _clipper_true_to_mean  # type: ignore
            except Exception:
                pass

        # (D) fallback inline converter (only used if we couldn't find your clipper one)
        if _clipper_true_to_mean is None:
            import math

            def _fallback_true_to_mean(nu_deg: float, e: float) -> float:
                nu = math.radians(nu_deg)
                t = math.tan(nu / 2.0)
                fac = math.sqrt((1.0 - e) / (1.0 + e))
                E = 2.0 * math.atan2(fac * t, 1.0)
                M = E - e * math.sin(E)
                return math.degrees(M) % 360.0

            _clipper_true_to_mean = _fallback_true_to_mean
            print("[ScalarPanel] WARNING: could not find Clipper ν→M converter; using fallback.")

        # ----------------- choose eccentricity consistent with your app -----------------
        def _get_default_e():
            # Prefer whatever your panel already has
            for attr in ("eccentricity", "e"):
                try:
                    v = getattr(self, attr, None)
                    if v is not None:
                        return float(v)
                except Exception:
                    pass

            # If you store sat/orbit params, try those
            sat = getattr(self, "sat_params", None)
            if sat is not None:
                for attr in ("eccentricity", "e"):
                    try:
                        v = getattr(sat, attr, None)
                        if v is not None:
                            return float(v)
                    except Exception:
                        pass

            # Final fallback (only if nothing else exists)
            return 0.0094

        default_e = _get_default_e()

        juice_meta: dict[str, dict] = {}
        for r in rows or []:
            if not isinstance(r, dict):
                continue

            # Accept both normalized and raw header keys
            observer = _clean(r.get("observer") or r.get("Observer")) or "UNKNOWN"
            enc = _clean(r.get("encounter") or r.get("Encounter"))
            utc = _clean(r.get("utc_iso") or r.get("Start(ISO)") or r.get("start_iso") or r.get("start"))

            # True anomaly from file
            nu = _to_float(r.get("true_anom_deg") or r.get("true anom deg") or r.get("nu_deg"))
            e = _to_float(r.get("eccentricity"))
            if e is None:
                e = default_e

            if not enc:
                continue

            # Naming convention requested: JUICE-E1
            enc_id = f"JUICE-{enc}"

            # Keep ALL original columns
            meta = dict(r)

            # Normalize fields we rely on (match plume style)
            meta["kind"] = "juice"
            meta["encounter_tag"] = "juice"
            meta["encounter_id"] = enc_id
            meta["utc_iso"] = utc
            meta["observer"] = observer

            # Normalize phase values
            if nu is not None:
                nu = nu % 360.0
                meta["true_anom_deg"] = nu
                meta["nu_deg"] = nu
            else:
                meta["true_anom_deg"] = None
                meta["nu_deg"] = None

            # Compute M using the same converter style as Clipper
            M = None
            if (nu is not None) and (e is not None):
                try:
                    M = _clipper_true_to_mean(nu, e)
                except Exception:
                    M = None

            if M is not None:
                try:
                    M = float(M) % 360.0
                except Exception:
                    M = None

            if M is not None:
                meta["mean_anom_deg"] = M
                meta["M_deg"] = M
            else:
                meta["mean_anom_deg"] = None
                meta["M_deg"] = None

            meta["eccentricity"] = e
            meta["phase_src"] = _clean(meta.get("phase_src")) or "file_true_to_mean"

            juice_meta[enc_id] = meta

        self._juice_meta_by_id = juice_meta
        self._juice_loaded = True
        print(f"[ScalarPanel] juice_meta_by_id size={len(self._juice_meta_by_id)}")



    def _on_enc_choice(self, evt):
        enc_id = self.cmb_enc.GetStringSelection()
        self._selected_encounter_id = enc_id

        meta = self._enc_meta_by_id.get(enc_id, {})

        # If this is a plume observation, do NOT use encounter CA metadata.
        if enc_id.startswith("PLUME_"):
            # Defer: resolve orbital state only for this plume selection.
            self._selected_encounter_M_deg = None
            self._selected_encounter_ca_utc = meta.get("utc_iso") or meta.get("utc_ca")

            try:
                self._ensure_plume_state_for_selected(enc_id)  # new helper (below)
            except Exception:
                pass

        else:
            # Normal Clipper encounter behavior
            M_ca = meta.get("M_ca_deg")
            if M_ca is None:
                M_ca = self._lookup_encounter_M(enc_id)
            self._selected_encounter_M_deg = None if M_ca is None else float(M_ca) % 360.0

            self._selected_encounter_ca_utc = meta.get("utc_ca")

        try:
            self.clear_secondary_marker()
        except Exception:
            pass

        self._update_return_button_state()


    def _lookup_encounter_M(self, enc_id: str) -> float | None:
        """
        Return the mean anomaly M (deg) for enc_id, or None if not available.
        Update this to match your data source.
        """
        try:
            # 1) If you keep a DataFrame around:
            df = getattr(self, "_encounters_df", None)
            if df is not None:
                # Try M directly
                if "M_deg" in df.columns:
                    row = df.loc[df["enc_label"] == enc_id]
                    if not row.empty:
                        return float(row.iloc[0]["M_deg"]) % 360.0
                # Or convert from true anomaly if available
                if "true_anom_deg" in df.columns:
                    from .utils import true_to_mean_anomaly_deg as _true_to_mean_anomaly_deg
                    row = df.loc[df["enc_label"] == enc_id]
                    if not row.empty:
                        M = _true_to_mean_anomaly_deg(float(row.iloc[0]["true_anom_deg"]))
                        return float(M) % 360.0

            # 2) If another panel knows (e.g., a point/encounter panel API):
            pp = getattr(self, "point_panel", None)
            if pp and hasattr(pp, "get_selected_encounter_M"):
                return float(pp.get_selected_encounter_M(enc_id)) % 360.0

            # 3) No idea
            return None
        except Exception:
            return None
        
    def set_current_M(self, M_deg: float):
        """
        Move to the exact requested mean anomaly M (deg), without snapping
        the stress computation to the nearest slider sample.
        """
        try:
            target = float(M_deg) % 360.0
            self._M_override = target

            # Keep slider roughly in sync visually, if desired
            try:
                self._idx = int(self._closest_index_for_M(target))
                self.sld_orbit.SetValue(self._idx)
            except Exception:
                pass

            self._evaluate_and_plot()
        except Exception as e:
            print("ScalarPlotPanel.set_current_M error:", e)
        

    def _on_orbit_M_changed(self, M_deg: float | None):
        """Register the currently displayed mean anomaly; call after any M update."""
        self._current_M_deg = None if M_deg is None else float(M_deg) % 360.0
        self._update_return_button_state()

    def _nudge_time_minutes(self, dmin: int):
        """
        Shift the displayed time relative to closest approach by Δt (minutes).
        Converts Δt → ΔM using mean motion, then snaps slider to nearest sample.
        """

        def _dbg(msg, *a):
            print("[StressViz/_nudge_time_minutes] " + (msg % a if a else msg))

        try:
            enc_id = getattr(self, "_selected_encounter_id", None)
            meta   = getattr(self, "_enc_meta_by_id", {}).get(enc_id, {}) if enc_id else {}

            M_ca = meta.get("M_ca_deg")
            P_h  = meta.get("period_hours")

            # Fallbacks
            if P_h is None or not (isinstance(P_h, (int, float)) and P_h > 0):
                P_h = (
                    getattr(self, "period_hours", None) or
                    getattr(self, "default_period_hours", None) or
                    85.228  # Europa default
                )

            if M_ca is None:
                # Last-known selection value if present
                M_ca = getattr(self, "_selected_encounter_M_deg", None)

            if M_ca is None:
                # As a last resort, try to derive from currently shown ν (if the widget exists)
                try:
                    e = float(self._get_eccentricity()) if hasattr(self, "_get_eccentricity") else 0.0
                except Exception:
                    e = 0.0
                try:
                    from .utils import true_to_mean_anomaly_deg as _true2mean
                    nu_now = None
                    if hasattr(self, "txt_nu") and self.txt_nu:
                        try:
                            nu_now = float(self.txt_nu.GetValue())
                        except Exception:
                            nu_now = None
                    if nu_now is not None:
                        M_ca = float(_true2mean(nu_now, e)) % 360.0
                except Exception:
                    M_ca = None

            # If we still don't have what we need, warn & bail
            if M_ca is None or P_h is None or not (P_h > 0):
                import wx
                wx.LogWarning(
                    f"Time nudge unavailable: missing M_ca ({M_ca}) or period_hours ({P_h}) for this encounter."
                )
                _dbg("enc=%r M_ca=%r P_h=%r -> cannot nudge (dmin=%r)", enc_id, M_ca, P_h, dmin)
                return

            # Mean motion (deg/min) and new mean anomaly
            n_deg_per_min = 360.0 / (float(P_h) * 60.0)
            M_new = (float(M_ca) + n_deg_per_min * float(dmin)) % 360.0

            # Persist the updated M for the current encounter
            try:
                if enc_id:
                    md = getattr(self, "_enc_meta_by_id", {}).get(enc_id)
                    if isinstance(md, dict):
                        md["M_ca_deg"] = float(M_new)
                self._selected_encounter_M_deg = float(M_new)
            except Exception:
                pass

            _dbg("enc=%s M_ca=%.6f P_h=%.6f dmin=%+d → M_new=%.6f (n=%.6f deg/min)",
                enc_id, float(M_ca), float(P_h), int(dmin), float(M_new), float(n_deg_per_min))

            # Snap slider to nearest sample and trigger the usual pipeline
            k = self._closest_index_for_M(M_new)
            self.sld_orbit.SetValue(int(k))
            self._on_slide(None)

        except Exception:
            traceback.print_exc()



    def _wrapdiff_deg(self, a: float, b: float) -> float:
        """Smallest signed difference a-b in degrees, wrapped to [-180,180]."""
        d = (float(a) - float(b)) % 360.0
        if d > 180.0:
            d -= 360.0
        return d

    def _enc_meta(self) -> dict:
        enc_id = getattr(self, "_selected_encounter_id", None)
        if not enc_id:
            return {}
        enc_id = str(enc_id)

        if enc_id.startswith("PLUME_"):
            return (getattr(self, "_plume_meta_by_id", {}) or {}).get(enc_id, {}) or {}

        return (getattr(self, "_enc_meta_by_id", {}) or {}).get(enc_id, {}) or {}

    def _enc_display_label(self, e: dict) -> str:
        """What the user sees in the dropdown."""
        name = str(e.get("name") or e.get("id") or "")
        eid  = str(e.get("id") or "")
        # If this is a plume ID, hide the PLUME_ prefix in the label
        if eid.startswith("PLUME_"):
            # show everything after PLUME_
            return name.replace("PLUME_", "", 1) if name.startswith("PLUME_") else name
        return name

    def _enc_sort_key(self, e: dict):
        """Sort by spacecraft first; plumes should sort with their spacecraft."""
        eid = str(e.get("id") or "")
        label = self._enc_display_label(e)

        # If your IDs look like: PLUME_<Spacecraft>_<time> ...
        # this will group plumes under the same <Spacecraft> prefix.
        spacecraft = label.split("_", 1)[0].split("-", 1)[0].strip()

        # Put non-plumes before plumes if tied
        is_plume = 1 if eid.startswith("PLUME_") else 0
        return (spacecraft, is_plume, label)
    
    # --- Multi-encounter UI + handlers ------------------------------------------

    def _ensure_multi_enc_state(self) -> None:
        if not hasattr(self, "_selected_encounter_ids"):
            self._selected_encounter_ids = []
        if not hasattr(self, "_selected_enabled_by_id"):
            self._selected_enabled_by_id = {}
        if not hasattr(self, "_enc_display_by_id"):
            self._enc_display_by_id = {}

        # --- Always ensure Stress Plot row exists ---
        if STRESS_PLOT_ID not in self._selected_encounter_ids:
            self._selected_encounter_ids.insert(0, STRESS_PLOT_ID)
        if STRESS_PLOT_ID not in self._selected_enabled_by_id:
            self._selected_enabled_by_id[STRESS_PLOT_ID] = True  # default = show ν
        self._enc_display_by_id[STRESS_PLOT_ID] = STRESS_PLOT_LABEL

    def _enc_id_of(self, enc: dict) -> Optional[str]:
        eid = enc.get("id") or enc.get("encounter_id")
        eid = str(eid).strip() if eid is not None else ""
        return eid or None

    def _enc_display_label(self, enc: dict) -> str:
        name = str(enc.get("name") or "").strip()
        if name:
            return name
        eid = self._enc_id_of(enc)
        return eid or "Encounter"

    def _enc_by_id_local(self, enc_id: str) -> Optional[dict]:
        encs = getattr(self, "_encounters", []) or []
        for e in encs:
            if str(e.get("id") or e.get("encounter_id") or "") == str(enc_id):
                return e
        return None

    def _find_chk_index_by_enc_id(self, enc_id: str) -> Optional[int]:
        enc_id = str(enc_id)
        for i in range(self.chk_selected_enc.GetCount()):
            if str(self.chk_selected_enc.GetClientData(i) or "") == enc_id:
                return i
        return None

    def _sync_checklist_from_model(self) -> None:
        """Rebuild the CheckListBox to match the internal subset + enabled state."""
        self._ensure_multi_enc_state()

        self.chk_selected_enc.Freeze()
        try:
            self.chk_selected_enc.Clear()
            for enc_id in self._selected_encounter_ids:
                if enc_id == STRESS_PLOT_ID:
                    label = STRESS_PLOT_LABEL
                else:
                    enc = self._enc_by_id_local(enc_id) or {"id": enc_id}
                    label = self._enc_display_by_id.get(enc_id) or self._enc_display_label(enc)
                    self._enc_display_by_id[enc_id] = label

                self.chk_selected_enc.Append(label)
                idx = self.chk_selected_enc.GetCount() - 1
                self.chk_selected_enc.SetClientData(idx, enc_id)
                self.chk_selected_enc.Check(idx, bool(self._selected_enabled_by_id.get(enc_id, True)))
        finally:
            self.chk_selected_enc.Thaw()
    
    def _ensure_stress_plot_row(self):
        self._ensure_multi_enc_state()

        if STRESS_PLOT_ID not in self._selected_encounter_ids:
            self._selected_encounter_ids.insert(0, STRESS_PLOT_ID)

        # Default checked = ν marker shown 
        if STRESS_PLOT_ID not in self._selected_enabled_by_id:
            self._selected_enabled_by_id[STRESS_PLOT_ID] = True

        # Label cache
        self._enc_display_by_id[STRESS_PLOT_ID] = STRESS_PLOT_LABEL

    def _build_selected_encounter_rows(self, parent):
        sc = wx.ScrolledWindow(parent, style=wx.VSCROLL)
        sc.SetScrollRate(0, 12)

        self._enc_row_widgets = {}  # enc_id -> (chk, picker, label)

        s = wx.BoxSizer(wx.VERTICAL)
        sc.SetSizer(s)

        return sc
    
    def _sync_rows_from_model(self):
        """
        Rebuild the encounter rows from:
        self._selected_encounter_ids
        self._selected_enabled_by_id
        self._selected_color_by_id
        """

        sc = getattr(self, "enc_rows", None)
        if sc is None:
            return

        s = sc.GetSizer()
        if s is None:
            s = wx.BoxSizer(wx.VERTICAL)
            sc.SetSizer(s)

        # Clear existing widgets
        try:
            s.Clear(delete_windows=True)
        except Exception:
            for child in sc.GetChildren():
                try:
                    child.Destroy()
                except Exception:
                    pass
            s = wx.BoxSizer(wx.VERTICAL)
            sc.SetSizer(s)

        self._enc_row_widgets = {}

        for enc_id in list(self._selected_encounter_ids):
            eid = str(enc_id)
            row = wx.BoxSizer(wx.HORIZONTAL)

            chk = wx.CheckBox(sc)
            chk.SetValue(bool(self._selected_enabled_by_id.get(eid, True)))

            label = self._enc_display_by_id.get(eid)
            if not label:
                enc = self._enc_by_id_local(eid) or {"id": eid}
                label = self._enc_display_label(enc)
                self._enc_display_by_id[eid] = label

            lab = wx.StaticText(sc, label=label)

            row.Add(chk, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
            row.Add(lab, 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)

            # --- Special case: Stress Plot row has no color picker ---
            if eid == STRESS_PLOT_ID:
                gray_lbl = wx.StaticText(sc)
                row.Add(gray_lbl, 0, wx.ALIGN_CENTER_VERTICAL)
                self._enc_row_widgets[eid] = (chk, None, lab)

                def _on_check(evt, _eid=eid):
                    self._selected_enabled_by_id[_eid] = bool(evt.GetEventObject().GetValue())
                    evt.Skip()

                chk.Bind(wx.EVT_CHECKBOX, _on_check)

                s.Add(row, 0, wx.EXPAND | wx.ALL, 2)
                continue

            # --- Normal encounter rows ---
            picker = wx.ColourPickerCtrl(sc)
            picker.SetMinSize((40, -1))

            col = self._get_or_create_encounter_color(eid)
            try:
                picker.SetColour(col)
            except Exception:
                pass

            row.Add(wx.StaticText(sc, label="Color:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
            row.Add(picker, 0, wx.ALIGN_CENTER_VERTICAL)

            s.Add(row, 0, wx.EXPAND | wx.ALL, 2)
            self._enc_row_widgets[eid] = (chk, picker, lab)

            def _on_check(evt, _eid=eid):
                self._selected_enabled_by_id[_eid] = bool(evt.GetEventObject().GetValue())
                evt.Skip()

            def _on_color(evt, _eid=eid):
                col = evt.GetEventObject().GetColour()
                self._selected_color_by_id[_eid] = col.GetAsString(wx.C2S_HTML_SYNTAX)
                evt.Skip()

            chk.Bind(wx.EVT_CHECKBOX, _on_check)
            picker.Bind(wx.EVT_COLOURPICKER_CHANGED, _on_color)

        # Keep Stress Plot out of the color dict entirely
        self._selected_color_by_id.pop(STRESS_PLOT_ID, None)

        sc.Layout()
        s.FitInside(sc)
        sc.Refresh()

    def _get_or_create_encounter_color(self, enc_id: str) -> str:
        eid = str(enc_id)

        col = self._selected_color_by_id.get(eid)
        if col:
            return col

        import random
        r = random.randint(64, 230)
        g = random.randint(64, 230)
        b = random.randint(64, 230)
        col = f"#{r:02X}{g:02X}{b:02X}"

        self._selected_color_by_id[eid] = col
        return col
    
    def _marker_for_encounter(self, enc: dict) -> str:
        """
        Return a matplotlib marker symbol based on the observing platform/source.
        """
        text_parts = [
            str(enc.get("observer") or ""),
            str(enc.get("platform") or ""),
            str(enc.get("mission") or ""),
            str(enc.get("spacecraft") or ""),
            str(enc.get("instrument") or ""),
            str(enc.get("name") or ""),
            str(enc.get("id") or ""),
            str(enc.get("encounter_id") or ""),
        ]
        txt = " ".join(text_parts).lower()

        if "clipper" in txt or "europa clipper" in txt:
            return "o"
        if "juice" in txt:
            return "^"
        if "hst" in txt or "hubble" in txt:
            return "s"
        if "keck" in txt:
            return "*"

        # fallback
        return "o"

    def _on_choose_multiple_encounters(self, _evt=None) -> None:
        """Add encounters into the subset (does not remove existing ones)."""
        self._ensure_multi_enc_state()

        encs = getattr(self, "_encounters", []) or []
        if not encs:
            wx.MessageBox("No encounters available.", "Encounters",
                        wx.OK | wx.ICON_INFORMATION)
            return

        ids: list[str] = []
        labels: list[str] = []
        for e in encs:
            eid = self._enc_id_of(e)
            if not eid:
                continue
            ids.append(eid)
            labels.append(self._enc_display_label(e))

        if not ids:
            wx.MessageBox("No encounters available.", "Encounters",
                        wx.OK | wx.ICON_INFORMATION)
            return

        prev_set = set(self._selected_encounter_ids)
        preselect = [i for i, eid in enumerate(ids) if eid in prev_set]

        dlg = wx.MultiChoiceDialog(
            self,
            message="Select encounters to add to the list (checkboxes below control what plots).",
            caption="Add Encounters",
            choices=labels,
        )
        try:
            if preselect:
                dlg.SetSelections(preselect)
            if dlg.ShowModal() != wx.ID_OK:
                return
            sel_idx = list(dlg.GetSelections() or [])
        finally:
            dlg.Destroy()

        picked_ids = [ids[i] for i in sel_idx if 0 <= i < len(ids)]
        if not picked_ids:
            return

        for eid in picked_ids:
            if eid in prev_set:
                continue
            prev_set.add(eid)
            self._selected_encounter_ids.append(eid)
            self._selected_enabled_by_id[eid] = True

            enc = self._enc_by_id_local(eid) or {"id": eid}
            self._enc_display_by_id[eid] = self._enc_display_label(enc)

        self._sync_rows_from_model()

    def _on_toggle_selected_encounter(self, evt) -> None:
        self._ensure_multi_enc_state()

        idx = evt.GetInt()
        enc_id = self.chk_selected_enc.GetClientData(idx)
        if not enc_id:
            return
        enc_id = str(enc_id)

        checked = bool(self.chk_selected_enc.IsChecked(idx))
        self._selected_enabled_by_id[enc_id] = checked

        if enc_id == STRESS_PLOT_ID:
            op = getattr(self, "orb_panel", None)
            if op and hasattr(op, "set_show_nu_marker"):
                op.set_show_nu_marker(checked)

    def _on_plot_selected_encounters(self, _evt=None) -> None:
        """Plot all CHECKED encounters from the subset (refreshes every button click)."""
        self._ensure_multi_enc_state()

        if not self._selected_encounter_ids:
            wx.MessageBox(
                "No encounters in the list. Use “Select encounters” first.",
                "Plot on Orbit",
                wx.OK | wx.ICON_INFORMATION,
            )
            return

        op = getattr(self, "orb_panel", None)
        if not op:
            return

        # Always apply Stress Plot toggle to ν marker on plot clicks
        show_nu = bool(self._selected_enabled_by_id.get(STRESS_PLOT_ID, True))
        if hasattr(op, "set_show_nu_marker"):
            op.set_show_nu_marker(show_nu)

        # Only real encounters should be plotted as secondary markers
        enabled_ids = [
            str(eid) for eid in self._selected_encounter_ids
            if str(eid) != STRESS_PLOT_ID and self._selected_enabled_by_id.get(str(eid), True)
        ]

        if not enabled_ids:
            wx.MessageBox(
                "No encounters are checked to plot.",
                "Plot on Orbit",
                wx.OK | wx.ICON_INFORMATION,
            )
            return

        # Clear first so every click reflects current checkbox state
        if hasattr(op, "clear_secondary_markers"):
            op.clear_secondary_markers()
        elif hasattr(op, "clear_secondary"):
            op.clear_secondary()

        supports_multi = hasattr(op, "add_secondary_M")
        supports_single = hasattr(op, "set_secondary_M")

        any_ok = False
        failed: list[str] = []

        def _get_or_create_color(eid: str) -> str:
            col = self._selected_color_by_id.get(eid)
            if col:
                return col

            palette = [
                "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            ]
            try:
                idx = max(0, self._selected_encounter_ids.index(eid))
            except Exception:
                idx = 0
            col = palette[idx % len(palette)]
            self._selected_color_by_id[eid] = col
            return col

        for eid in enabled_ids:
            enc = self._enc_by_id_local(eid) or {"id": eid}
            label = self._enc_display_by_id.get(eid) or self._enc_display_label(enc)
            self._enc_display_by_id[eid] = label

            M = self._M_for_encounter(enc)
            if M is None or not np.isfinite(M):
                failed.append(label)
                continue

            color = _get_or_create_color(eid)
            marker = self._marker_for_encounter(enc)

            if supports_multi:
                try:
                    op.add_secondary_M(float(M), label=label, color=color, marker=marker)
                    any_ok = True
                except TypeError:
                    try:
                        op.add_secondary_M(float(M), label=label, color=color)
                        any_ok = True
                    except Exception:
                        failed.append(label)
                except Exception:
                    failed.append(label)
                continue

            if supports_single:
                try:
                    op.set_secondary_M(float(M), label=label, color=color, marker=marker)
                    any_ok = True
                    break
                except TypeError:
                    try:
                        op.set_secondary_M(float(M), label=label, color=color)
                        any_ok = True
                        break
                    except Exception:
                        failed.append(label)
                except Exception:
                    failed.append(label)
                continue

            failed.append(label)

        try:
            self._sync_rows_from_model()
        except Exception:
            pass

        try:
            self._update_orbit_side_legend()
        except Exception:
            pass

        if not any_ok:
            wx.MessageBox(
                "No valid checked encounters to plot.\n\n"
                + ("Mean anomaly (M) failed for:\n" + "\n".join(failed[:15]) if failed else ""),
                "Nothing to plot",
                wx.OK | wx.ICON_INFORMATION,
            )
        elif failed:
            wx.MessageBox(
                "Some encounters could not resolve mean anomaly (M):\n\n"
                + "\n".join(failed[:15]) + ("\n…" if len(failed) > 15 else ""),
                "Partial Plot",
                wx.OK | wx.ICON_WARNING,
            )

    def _on_clear_selected_encounters(self, _evt=None) -> None:
        self._ensure_multi_enc_state()

        self._selected_encounter_ids.clear()
        self._selected_enabled_by_id.clear()
        self._enc_display_by_id.clear()

        self.chk_selected_enc.Clear()

        op = getattr(self, "orb_panel", None)
        if op and hasattr(op, "clear_secondary_markers"):
            try:
                op.clear_secondary_markers()
            except Exception:
                pass
        elif op and hasattr(op, "clear_secondary"):
            try:
                op.clear_secondary()
            except Exception:
                pass


    def _M_from_delta_minutes(self, dmin: float):
        import math
        meta = self._enc_meta()
        # fallbacks:
        M_ca = meta.get("M_ca_deg")
        if M_ca is None:
            M_ca = self._selected_encounter_M_deg   # use the Return-to-Encounter M if set
        P_h = meta.get("period_hours")
        if P_h is None:
            P_h = getattr(self, "period_hours", None)  # allow a global period on the panel
        if P_h is None:
            P_h = 85.228  # Europa default fallback

        if M_ca is None or P_h is None:
            return None

        n = 2*math.pi / (float(P_h)*3600.0)  # rad/s
        dM_deg = math.degrees(n * (float(dmin)*60.0))
    def _minutes_from_M(self, M_deg: float):
        import math
        meta = self._enc_meta()
        M_ca = meta.get("M_ca_deg")
        if M_ca is None:
            M_ca = self._selected_encounter_M_deg
        P_h = meta.get("period_hours")
        if P_h is None:
            P_h = getattr(self, "period_hours", None)
        if P_h is None:
            P_h = 85.228

        if M_ca is None or P_h is None:
            return None

        # smallest signed ΔM in degrees
        d = (float(M_deg) - float(M_ca)) % 360.0
        if d > 180.0:
            d -= 360.0

        n = 2*math.pi / (float(P_h)*3600.0)
        dt_sec = math.radians(d) / n
        return dt_sec / 60.0



    def _closest_index_for_M(self, M_target_deg: float) -> int:
        """Return index k whose mapped M(nu_k) is closest to M_target."""
        if not self._nus:
            return 0
        e = self._ecc_safe()
        Ms = [float(_true_to_mean_anomaly_deg(nu, e)) % 360.0 for nu in self._nus]
        diffs = [abs(self._wrapdiff_deg(M, M_target_deg)) for M in Ms]
        return int(np.argmin(diffs))
    
    def _in_periodic_window(self, M, M0, win, period=360.0):
        M0 = float(M0) % period
        M = float(M) % period
        lo = (M0 - win) % period
        hi = (M0 + win) % period

        if lo <= hi:
            return lo <= M <= hi
        else:
            return (M >= lo) or (M <= hi)
        
    def _get_target_encounter_from_checked_subset(self):
        """
        Returns: (enc_id, enc_dict) for the first CHECKED real encounter
        in the current selected subset, or (None, None) if none.
        """
        self._ensure_multi_enc_state()

        for enc_id in self._selected_encounter_ids:
            enc_id = str(enc_id)

            # skip the synthetic Stress Plot row
            if enc_id == STRESS_PLOT_ID:
                continue

            if not self._selected_enabled_by_id.get(enc_id, True):
                continue

            enc = self._enc_by_id_local(enc_id)
            if enc is None and hasattr(self, "_enc_by_id"):
                try:
                    enc = self._enc_by_id(enc_id)
                except Exception:
                    enc = None

            return enc_id, enc

        return None, None
    
    def _get_scalar_panel_ref(self):
        sp = getattr(self, "scalar_panel_ref", None)
        if sp is not None:
            return sp

        try:
            host = wx.GetTopLevelParent(self)
        except Exception:
            host = None

        if host is not None:
            sp = getattr(host, "scalar_panel_ref", None)
            if sp is not None:
                return sp

        return None
    
    def _on_show_nearby_events(self, _evt=None):
        self._ensure_multi_enc_state()

        encs = getattr(self, "_encounters", []) or []
        if not encs:
            wx.MessageBox(
                "No encounters available.",
                "Nearby events",
                wx.OK | wx.ICON_INFORMATION,
                self,
            )
            return

        # Target = first checked item in the selected-encounters rectangle
        target_id, target = self._get_target_encounter_from_checked_subset()
        if not target_id or not target:
            wx.MessageBox(
                "Check exactly one reference encounter in the Selected Encounters list first.",
                "Nearby events",
                wx.OK | wx.ICON_INFORMATION,
                self,
            )
            return

        M0 = self._M_for_encounter(target)
        if M0 is None or not np.isfinite(M0):
            wx.MessageBox(
                "Could not compute mean anomaly (M) for the reference encounter.",
                "Nearby events",
                wx.OK | wx.ICON_WARNING,
                self,
            )
            return

        M0 = float(M0) % 360.0


        # ------------------------------------------------------------
        # Find nearby encounters around the reference M
        # ------------------------------------------------------------
        win = 10.0  # degrees, half-width

        nearby_ids = []
        for e in encs:
            eid = self._enc_id_of(e)
            if not eid:
                continue

            M = self._M_for_encounter(e)
            if M is None or not np.isfinite(M):
                continue

            if self._in_periodic_window(M, M0, win):
                nearby_ids.append(str(eid))

        if not nearby_ids:
            wx.MessageBox(
                f"No encounters within ±{win:.0f}° of the reference encounter.",
                "Nearby events",
                wx.OK | wx.ICON_INFORMATION,
                self,
            )
            return

        # Add to subset + enable them
        prev = set(self._selected_encounter_ids)

        for eid in nearby_ids:
            if eid not in prev:
                self._selected_encounter_ids.append(eid)
                prev.add(eid)

            self._selected_enabled_by_id[eid] = True

            enc = self._enc_by_id_local(eid) or {"id": eid}
            self._enc_display_by_id[eid] = self._enc_display_label(enc)

        # Rebuild checklist so newly added items show up checked
        self._sync_rows_from_model()

        # Plot immediately
        self._on_plot_selected_encounters(None)

        wx.MessageBox(
            f"Added {len(nearby_ids)} encounter(s) within ±{win:.0f}° of "
            f"{self._enc_display_by_id.get(target_id, target_id)}.",
            "Nearby events",
            wx.OK | wx.ICON_INFORMATION,
            self,
        )

    def _on_move_stress_plot_to_selected_encounter(self, _evt=None):
        self._ensure_multi_enc_state()

        target_id, target = self._get_target_encounter_from_checked_subset()
        if not target_id or not target:
            wx.MessageBox(
                "Check exactly one reference encounter in the Selected Encounters list first.",
                "Move Stress Plot",
                wx.OK | wx.ICON_INFORMATION,
                self,
            )
            return

        M0 = self._M_for_encounter(target)
        if M0 is None or not np.isfinite(M0):
            wx.MessageBox(
                "Could not compute mean anomaly (M) for the selected encounter.",
                "Move Stress Plot",
                wx.OK | wx.ICON_WARNING,
                self,
            )
            return

        M0 = float(M0) % 360.0

        sp = self._get_scalar_panel_ref()
        if sp is None:
            try:
                host = wx.GetTopLevelParent(self)
            except Exception:
                host = None
            if host is not None:
                sp = getattr(host, "scalar_panel_ref", None)

        if sp is None:
            wx.MessageBox(
                "Stress plot window is not open.",
                "Move Stress Plot",
                wx.OK | wx.ICON_INFORMATION,
                self,
            )
            return

        try:
            sp.set_current_M(M0)
        except Exception as e:
            wx.MessageBox(
                f"Could not move stress plot.\n\n{e}",
                "Move Stress Plot",
                wx.OK | wx.ICON_WARNING,
                self,
            )
            return

        try:
            self._on_orbit_M_changed(M0)
        except Exception:
            pass

        try:
            op = getattr(self, "orb_panel", None)
            if op is not None and hasattr(op, "set_show_nu_marker"):
                op.set_show_nu_marker(True)
        except Exception:
            pass

        try:
            self._update_orbit_side_legend()
        except Exception:
            pass

    def wire_satstress(
        self,
        get_satellite: Callable[[], object],
        Diurnal_cls: type,
        lats_deg: np.ndarray,
        lons_deg: np.ndarray,
        nus: Optional[np.ndarray] = None,
        title: Optional[str] = None,
    ):
        """
        One-liner: plug SatStress into this panel.
        NOTE: We pass MEAN anomaly (M) to calc; we still index with ν samples.
        """
        self.set_axes(lats_deg, lons_deg)
        if nus is None:
            nus = np.linspace(0.0, 360.0, 13)
        nus = np.asarray(nus, float)

        StressCalc = resolve_stresscalc()

        def _eval_fn(M_deg: float):
            sat = get_satellite()
            if sat is None:
                raise RuntimeError("No satellite set.")
            diurn = Diurnal_cls(sat)

            holder = getattr(diurn, "stresses", diurn)
            for k, v in (("diurnal", True), ("tidal", True), ("nsr", False),
                        ("polar_wander", False), ("obliquity", False)):
                if hasattr(holder, k):
                    try: setattr(holder, k, bool(v))
                    except Exception: pass

            StressCalc = resolve_stresscalc()
            calc = StressCalc([diurn])

            lat = np.asarray(self._lats, float)
            lon = np.asarray(self._lons, float)
            TH, PH = np.meshgrid(np.radians(90.0 - lat), np.radians(lon))

            # Map mean anomaly to time (sec) via diurnal mean motion
            t_sec = _nu_to_time_seconds(M_deg, diurn)
            return _tensor_grid(calc, TH, PH, t_sec)  # Pa

        # NEW: make sure ν→M uses a real e right away
        def _ecc_provider():
            try:
                sat = get_satellite()
                if sat is None:
                    return 0.0
                d = Diurnal_cls(sat)
                return float(getattr(d, "e", 0.0))
            except Exception:
                return 0.0

        self._get_eccentricity = _ecc_provider  # <-- this powers _ecc_safe()

        self.set_axes(lats_deg, lons_deg)
        if nus is None:
            nus = np.linspace(0.0, 360.0, 13)
        nus = np.asarray(nus, float)

        self.bind_orbit_series(
            list(nus), _eval_fn, initial_nu_deg=float(nus[0]), title=title
        )

    # Back-compat helpers
    def show_components(self, lats_deg, lons_deg, Ttt, Tpt, Tpp, title: Optional[str] = None):
        self.set_axes(lats_deg, lons_deg)
        self._Z_custom = None
        self._T = (np.asarray(Ttt, float), np.asarray(Tpt, float), np.asarray(Tpp, float))
        self._field = "tens"
        self.choice_field.SetSelection(0)
        if title:
            self.ax.set_title(title)
        self._replot()

    def plot_scalar_map(self, lats_deg, lons_deg, Z_kPa, title: Optional[str] = None):
        self.set_axes(lats_deg, lons_deg)
        self._T = None
        self._Z_custom = np.asarray(Z_kPa, float)
        self._field = "custom"
        self.choice_field.SetSelection(4)
        if title:
            self.ax.set_title(title)
        self._replot()

    def _plot_vectors(self):
        """
        Draw σ1 and σ3 line glyphs at sampled grid points.
        """
        if getattr(self, "chk_principal", None) and not self.chk_principal.GetValue():
            return
        if self._T is None or self._lats is None or self._lons is None:
            return

        Ttt, Tpt, Tpp = self._T
        Ny, Nx = Ttt.shape

        m = 0.5 * (Ttt + Tpp)
        r = np.hypot(0.5 * (Ttt - Tpp), Tpt)
        s1 = m + r  # Pa
        s3 = m - r  # Pa
        ang = 0.5 * np.arctan2(2.0 * Tpt, (Ttt - Tpp))
        ang_perp = ang + np.pi / 2.0

        nlat = int(getattr(self, "vec_n_lat", getattr(self, "_vec_nlat", 10)))
        nlon = int(getattr(self, "vec_n_lon", getattr(self, "_vec_nlon", 10)))
        i_idx = np.linspace(0, Ny - 1, max(2, min(nlat, Ny)), dtype=int)
        j_idx = np.linspace(0, Nx - 1, max(2, min(nlon, Nx)), dtype=int)
        ii, jj = np.meshgrid(i_idx, j_idx, indexing="ij")

        s1_kpa = (s1[ii, jj]) * _KPA
        s3_kpa = (s3[ii, jj]) * _KPA
        ang1   = ang[ii, jj]
        ang3   = ang_perp[ii, jj]

        lats = self._lats; lons = self._lons
        if lats[0] > lats[-1]: lats = lats[::-1]
        if lons[0] > lons[-1]: lons = lons[::-1]
        X, Y = np.meshgrid(lons, lats, indexing="xy")
        Xs = X[ii, jj].astype(float)
        Ys = Y[ii, jj].astype(float)

        width_deg = abs(float(lons[-1]) - float(lons[0])) or 360.0
        base_len  = 0.10 * width_deg
        ref_kpa   = 100.0

        hl1 = 0.5 * base_len * np.clip(np.abs(s1_kpa) / ref_kpa, 0.0, 10.0)
        hl3 = 0.5 * base_len * np.clip(np.abs(s3_kpa) / ref_kpa, 0.0, 10.0)

        flip = (self.rb_dir.GetSelection() == 1) if hasattr(self, "rb_dir") else (not self._east_positive)
        if flip:
            xmin, xmax = float(lons[0]), float(lons[-1])
            Xs = xmax - (Xs - xmin)
            sign_x = -1.0
        else:
            sign_x = 1.0

        dx1 = sign_x * hl1 * np.cos(ang1);  dy1 = hl1 * np.sin(ang1)
        dx3 = sign_x * hl3 * np.cos(ang3);  dy3 = hl3 * np.sin(ang3)

        def _segments(dx, dy):
            x0 = (Xs - dx).ravel(); y0 = (Ys - dy).ravel()
            x1 = (Xs + dx).ravel(); y1 = (Ys + dy).ravel()
            return np.stack([np.stack([x0, y0], 1),
                            np.stack([x1, y1], 1)], 1)

        seg1 = _segments(dx1, dy1)
        seg3 = _segments(dx3, dy3)

        colors1 = np.where((s1_kpa.ravel() >= 0.0), 'red', 'blue').tolist()
        colors3 = np.where((s3_kpa.ravel() >= 0.0), 'red', 'blue').tolist()

        xlim = self.ax.get_xlim(); ylim = self.ax.get_ylim()

        for art in getattr(self, "_principal_artists", []):
            try: art.remove()
            except Exception: pass
        self._principal_artists = []

        lw     = 1.0
        under  = 1.0
        alphaU = 0.35
        z_fg   = 10
        z_bg   = 9

        arts = []
        show_s1 = (getattr(self, "chk_s1", None) is None) or self.chk_s1.GetValue()
        show_s3 = (getattr(self, "chk_s3", None) is None) or self.chk_s3.GetValue()

        if show_s1 and seg1.size:
            bg1 = LineCollection(seg1, colors="#000000", linewidths=lw + under,
                                 capstyle='round', joinstyle='round', antialiased=True,
                                 zorder=z_bg, alpha=alphaU)
            fg1 = LineCollection(seg1, colors=colors1, linewidths=lw,
                                 capstyle='round', joinstyle='round', antialiased=True,
                                 zorder=z_fg)
            arts.append(self.ax.add_collection(bg1))
            arts.append(self.ax.add_collection(fg1))

        if show_s3 and seg3.size:
            bg3 = LineCollection(seg3, colors="#000000", linewidths=lw + under,
                                 capstyle='round', joinstyle='round', antialiased=True,
                                 zorder=z_bg, alpha=alphaU)
            fg3 = LineCollection(seg3, colors=colors3, linewidths=lw,
                                 capstyle='round', joinstyle='round', antialiased=True,
                                 zorder=z_fg)
            arts.append(self.ax.add_collection(bg3))
            arts.append(self.ax.add_collection(fg3))

        self._principal_artists = arts
        self.ax.set_xlim(xlim); self.ax.set_ylim(ylim)

    def _draw_scale_bar(self):
        """Prominent 100 kPa reference bar, left-bottom inside the axes."""
        if not hasattr(self, "_ref_artists"): self._ref_artists = []
        for a in self._ref_artists:
            try: a.remove()
            except Exception: pass
        self._ref_artists = []

        if not (getattr(self, "chk_principal", None) and self.chk_principal.GetValue()):
            self.canvas.draw_idle(); return

        xmin, xmax = self.ax.get_xlim(); ymin, ymax = self.ax.get_ylim()
        dx = xmax - xmin; dy = ymax - ymin

        ref_k = float(getattr(self, "_ref_kpa", 100.0))
        Lh = self._len_from_kpa(ref_k)
        if isinstance(Lh, np.ndarray):
            Lh = float(Lh)

        x_center = xmin + 0.22 * dx
        y_line   = ymin + 0.06 * dy
        x0 = x_center - Lh; x1 = x_center + Lh

        self._ref_artists.append(self.ax.plot([x0, x1], [y_line, y_line], color="black", lw=2.4, zorder=4)[0])
        cap = 0.012 * dy
        self._ref_artists.append(self.ax.plot([x0, x0], [y_line - cap, y_line + cap], color="black", lw=1.8, zorder=4)[0])
        self._ref_artists.append(self.ax.plot([x1, x1], [y_line - cap, y_line + cap], color="black", lw=1.8, zorder=4)[0])
        self._ref_artists.append(self.ax.text(x_center, y_line + 0.03 * dy, f"{int(ref_k)} kPa",
                                            ha="center", va="bottom", fontsize=9, color="black"))

        self.ax.set_xlim(xmin, xmax); self.ax.set_ylim(ymin, ymax)

    # ---------------- internals ----------------
    def _nudge(self, step: int):
        if self._nus is None or len(self._nus) == 0:
            return
        try:
            self._M_override = None
            self._idx = int(np.clip(self._idx + step, 0, len(self._nus) - 1))
            self.sld_orbit.SetValue(self._idx)
            self._evaluate_and_plot()
        except Exception as e:
            print("ScalarPlotPanel._nudge error:", e)

    def _on_slide(self, _evt):
        try:
            self._M_override = None
            self._idx = int(self.sld_orbit.GetValue())
            self._evaluate_and_plot()
        except Exception as e:
            print("ScalarPlotPanel._on_slide error:", e)

    def _evaluate_and_plot(self):
        """
        Compute (Ttt,Tpt,Tpp) for the current M and trigger a redraw.
        """
        if self._nus is None or self._eval_fn is None or self._lats is None or self._lons is None:
            return

        try:
            # Use exact override M if present; otherwise derive M from current slider index
            if getattr(self, "_M_override", None) is not None:
                M = float(self._M_override) % 360.0
            else:
                nu = float(self._nus[self._idx]) % 360.0
                e = self._ecc_safe()
                M = float(_true_to_mean_anomaly_deg(nu, e)) % 360.0

            # 1) label as M with 2 decimals
            try:
                self.lbl_nu.SetLabel(f"M: {M:.2f}°")
            except Exception:
                pass

            # 2) mini-orbit sync using M value
            try:
                self._update_orbit_viz(M)
            except Exception:
                pass

            # 3) compute stress components (Pa) for this M
            Ttt, Tpt, Tpp = self._eval_fn(M)
            self._T = (np.asarray(Ttt, float), np.asarray(Tpt, float), np.asarray(Tpp, float))
            self._Z_custom = None

            # 4) replot
            try:
                self._replot()
            except Exception as replot_err:
                print("ScalarPlotPanel._replot error:", replot_err)
                try:
                    self.canvas.draw_idle()
                except Exception:
                    pass

        except Exception as e:
            print("ScalarPlotPanel._evaluate_and_plot error:", e)

    def _update_orbit_side_legend(self):
        s = getattr(self, "orbit_legend_sizer", None)
        panel = getattr(self, "orbit_legend_panel", None)
        scroll = getattr(self, "orbit_legend_scroll", None)
        op = getattr(self, "orb_panel", None)

        if s is None or panel is None or scroll is None or op is None:
            return

        try:
            s.Clear(delete_windows=True)
        except Exception:
            for child in scroll.GetChildren():
                try:
                    child.Destroy()
                except Exception:
                    pass

        secondary = getattr(op, "_secondary_markers", []) or []
        seen = set()

        symbol_map = {
            "o": "●",
            "^": "▲",
            "s": "■",
            "*": "★",
            "+": "+",
            "x": "×",
        }

        for item in secondary:
            try:
                lbl = item[1] if len(item) > 1 else None
                color = item[2] if len(item) > 2 else "#000000"
                marker = item[3] if len(item) > 3 else "o"
            except Exception:
                continue

            label = (str(lbl).strip() if lbl is not None else "") or "Encounter"
            if label in seen:
                continue
            seen.add(label)

            row = wx.BoxSizer(wx.HORIZONTAL)

            sym = wx.StaticText(scroll, label=symbol_map.get(marker, "●"))
            try:
                sym.SetForegroundColour(wx.Colour(color))
            except Exception:
                pass

            txt = wx.StaticText(scroll, label=label)

            row.Add(sym, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.Add(txt, 1, wx.ALIGN_CENTER_VERTICAL)

            s.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 3)

        if (getattr(op, "_nu_deg", None) is not None) and bool(getattr(op, "_show_nu_marker", True)):
            row = wx.BoxSizer(wx.HORIZONTAL)

            sym = wx.StaticText(scroll, label="+")
            sym.SetForegroundColour(wx.Colour("#000000"))

            txt = wx.StaticText(scroll, label=f"M = {op._nu_deg:.2f}°")

            row.Add(sym, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.Add(txt, 1, wx.ALIGN_CENTER_VERTICAL)

            s.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 3)

        scroll.Layout()
        scroll.FitInside()
        panel.Layout()

    def _on_dir_change(self, _):
        self._east_positive = (self.rb_dir.GetSelection() == 0)
        self._replot()

    def _on_field_change(self, _):
        self._field = ["tens", "comp", "mean", "diff", "custom"][self.choice_field.GetSelection()]
        self._replot()

    def _on_auto_range(self, _):
        self._auto_range = self.chk_auto.IsChecked()
        self.sp_l.Enable(not self._auto_range)
        self.sp_u.Enable(not self._auto_range)
        self._replot()

    def _on_bounds_change(self, _evt):
        try:
            vmin = float(self.sp_l.GetValue())
            vmax = float(self.sp_u.GetValue())
        except Exception:
            return

        if vmin >= vmax:
            return

        self._vmin = vmin
        self._vmax = vmax
        self._replot()

    def _update_orbit_viz(self, M_deg: float):
        op = getattr(self, "orb_panel", None)
        if not op:
            return
        # Prefer set_M; else reuse set_nu to place the marker at angle M
        setter = getattr(op, "set_M", None) or getattr(op, "set_nu", None)
        if callable(setter):
            try:
                setter(float(M_deg))
            except Exception:
                pass
        try:
            self._update_orbit_side_legend()
        except Exception:
            pass


    def _apply_axes_layout(self):
        self.ax.set_position(self._ax_rect)
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        div = make_axes_locatable(self.ax)
        self._cax.set_axes_locator(div.new_locator(self.ax, "right",
                                                size=self._cbar_size,
                                                pad=self._cbar_pad))
        self.canvas.draw_idle()

    def _relayout_overlay(self, _evt=None):
        try:
            w, h = self._plate.GetClientSize()
            pad_lr = self.FromDIP(16)
            pad_bottom = self.FromDIP(10)
            ov_h = self.FromDIP(36)
            self._overlay.SetSize(w - 2*pad_lr, ov_h)
            self._overlay.SetPosition((pad_lr, h - ov_h - pad_bottom))
            self._overlay.Layout()
        except Exception:
            pass

    def _derive_field_kPa(self):
        if self._field == "custom":
            return None if self._Z_custom is None else np.asarray(self._Z_custom, float)
        if self._T is None:
            return None
        Ttt, Tpt, Tpp = self._T
        m = 0.5 * (Ttt + Tpp)
        r = np.hypot(0.5 * (Ttt - Tpp), Tpt)
        s1 = m + r
        s3 = m - r
        if self._field == "tens":
            Z = s1
        elif self._field == "comp":
            Z = s3
        elif self._field == "mean":
            Z = 0.5 * (s1 + s3)
        else:
            Z = (s1 - s3)
        return Z * _KPA  # kPa

    def _extent_and_data(self, Z):
        lats = self._lats
        lons = self._lons
        Zd = np.array(Z, copy=True)
        if lats.size > 1 and lats[1] < lats[0]:
            lats = lats[::-1]
            Zd = Zd[::-1, :]
        if lons.size > 1 and lons[1] < lons[0]:
            lons = lons[::-1]
            Zd = Zd[:, ::-1]
        if not self._east_positive:
            Zd = Zd[:, ::-1]
            extent = [float(lons[-1]), float(lons[0]), float(lats[0]), float(lats[-1])]
        else:
            extent = [float(lons[0]), float(lons[-1]), float(lats[0]), float(lats[-1])]
        return extent, Zd
    
    # --- vector magnitude → drawn half-length (in degrees of lon) ---
    def _len_from_kpa(self, kpa: float) -> float:
        try:
            ref = float(getattr(self, "_vec_ref_kpa", 100.0))
            vmax_half = float(getattr(self, "_vec_max_deg", 10.0))
            if ref <= 0:
                ref = 100.0
            L = vmax_half * min(abs(float(kpa)) / ref, 1.0)
            return float(L)
        except Exception:
            return float(getattr(self, "_vec_max_deg", 10.0))

    def _ensure_scale_ax(self):
        pos = self.ax.get_position()
        fig_w_in, fig_h_in = self.fig.get_size_inches()
        size_pct  = 9.0
        pad_inch  = 0.90
        h_frac   = (size_pct / 100.0) * pos.height
        pad_frac = (pad_inch / fig_h_in)
        y0 = pos.y0 - pad_frac - h_frac
        rect = [pos.x0, max(0.0, y0), pos.width, h_frac]
        if getattr(self, "_scale_ax", None) is None or self._scale_ax.figure is not self.fig:
            self._scale_ax = self.fig.add_axes(rect)
        else:
            self._scale_ax.set_position(rect)
        self._scale_ax.set_axis_off()
        return self._scale_ax

    def _clear_legend(self):
        for art in getattr(self, "_legend_art", []):
            try: art.remove()
            except Exception: pass
        self._legend_art = []

    def _draw_vector_scale_bar(self):
        self._ensure_scale_ax()
        ax = self._scale_ax
        ax.clear()
        ax.set_axis_off()

        ref_kpa  = float(getattr(self, "_vec_ref_kpa", 10.0))
        half_deg = float(self._len_from_kpa(ref_kpa))
        total_deg = 2.0 * half_deg

        try:
            x0, x1 = self.ax.get_xlim()
            map_width = abs(float(x1 - x0)) or 360.0
        except Exception:
            map_width = 360.0

        frac = max(0.18, min(0.42, total_deg / map_width))

        y            = 0.30
        label_y      = y - 0.12
        left_label_x = 0.02
        left_margin  = 0.14
        x0n          = left_margin
        x1n          = min(0.92, x0n + frac)
        right_label_x = min(0.98, x1n + 0.02)

        tick_h = 0.12
        lw     = 1.6

        ax.text(left_label_x, label_y, "Scale", ha="left", va="center", fontsize=9, color="black")
        ax.plot([x0n, x1n], [y, y], color="black", linewidth=lw, solid_capstyle="butt", clip_on=False)
        ax.plot([x0n, x0n], [y - tick_h/2, y + tick_h/2], color="black", linewidth=lw, clip_on=False)
        ax.plot([x1n, x1n], [y - tick_h/2, y + tick_h/2], color="black", linewidth=lw, clip_on=False)
        ax.text(right_label_x, y, f"{int(ref_kpa)} kPa", ha="left", va="center", fontsize=9, color="black")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    def _ensure_cbar_axes(self, pad_frac=0.015, cbar_w_frac=0.03):
        if hasattr(self, "_ax_rect"):
            self.ax.set_position(self._ax_rect)
        pos = self.ax.get_position()
        x0, y0, w, h = pos.x0, pos.y0, pos.width, pos.height
        cax_rect = [x0 + w + pad_frac, y0, cbar_w_frac, h]
        if getattr(self, "_cax", None) is None or self._cax.figure is not self.fig or self._cax not in self.fig.axes:
            self._cax = self.fig.add_axes(cax_rect)
        else:
            self._cax.set_position(cax_rect)

    def _replot(self):
        if self._field == "custom" and self._Z_custom is not None:
            Z_kPa = np.asarray(self._Z_custom, float)
        else:
            Z_kPa = self._derive_field_kPa()
        if Z_kPa is None:
            return

        if getattr(self, "_fixed_clim_kpa", None) is not None:
            vmin, vmax = map(float, self._fixed_clim_kpa)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax) if vmin < 0 < vmax else Normalize(vmin=vmin, vmax=vmax)
        else:
            vmin = float(self.sp_l.GetValue()); vmax = float(self.sp_u.GetValue())
            if vmin >= vmax:
                vmax = vmin + 1.0
                norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax) if vmin < 0 < vmax else Normalize(vmin=vmin, vmax=vmax)

        extent, Zdisp = self._extent_and_data(Z_kPa)

        first = (self._im is None)
        if first:
            self.ax.clear()
            self.ax.set_xlabel("Longitude [°E]" if self._east_positive else "Longitude [°W]", labelpad=-1.8)
            self.ax.set_ylabel("Latitude [°]")
            if hasattr(self, "_ax_rect"):
                self.ax.set_position(self._ax_rect)
            self._im = self.ax.imshow(
                Zdisp, extent=extent, origin="lower", aspect="auto",
                cmap="gist_rainbow_r", norm=norm, zorder=0
            )
        else:
            self._im.set_data(Zdisp)
            self._im.set_extent(extent)
            self._im.set_norm(norm)

        self._ensure_cbar_axes()

        if getattr(self, "_cbar", None) is None or self._cbar.ax is None:
            self._cbar = self.fig.colorbar(self._im, cax=self._cax, label="kPa")
        else:
            try:
                self._cbar.update_normal(self._im)
                self._cbar.ax.set_ylabel("kPa")
            except Exception:
                try: self._cbar.remove()
                except Exception: pass
                self._cbar = self.fig.colorbar(self._im, cax=self._cax, label="kPa")

        for coll in getattr(self, "_principal_artists", []):
            try: coll.remove()
            except Exception: pass
        self._principal_artists = []
        self._plot_vectors()

        try:
            self._draw_vector_scale_bar()
        except Exception as e:
            print("scale bar draw error:", e)

        self.canvas.draw_idle()

    # mouse readout
    def _on_motion(self, evt):
        if evt.inaxes != self.ax:
            return
        try:
            lon = float(evt.xdata)
            lat = float(evt.ydata)
        except Exception:
            return
        self.txt_lat.SetValue(f"{lat:7.2f}")
        self.txt_lon.SetValue(f"{lon:7.2f}")

        Z_kPa = None
        if self._field == "custom" and self._Z_custom is not None:
            Z_kPa = np.asarray(self._Z_custom, float)
        else:
            Z_kPa = self._derive_field_kPa()
        if Z_kPa is None:
            self.txt_val.SetValue("")
            return

        lats = self._lats
        lons = self._lons
        Z = Z_kPa.copy()
        if lats[0] > lats[-1]:
            lats = lats[::-1]
            Z = Z[::-1, :]
        if lons[0] > lons[-1]:
            lons = lons[::-1]
            Z = Z[:, ::-1]
        if not (lats[0] <= lat <= lats[-1]) or not (lons[0] <= lon <= lons[-1]):
            self.txt_val.SetValue("")
            return
        i = np.interp(lat, lats, np.arange(lats.size))
        j = np.interp(lon, lons, np.arange(lons.size))
        i0, j0 = int(np.floor(i)), int(np.floor(j))
        i1, j1 = min(i0 + 1, lats.size - 1), min(j0 + 1, lons.size - 1)
        di, dj = i - i0, j - j0
        v = ((1 - di) * (1 - dj) * Z[i0, j0] +
             (1 - di) * dj * Z[i0, j1] +
             di * (1 - dj) * Z[i1, j0] +
             di * dj * Z[i1, j1])
        self.txt_val.SetValue(f"{v:9.2f}")

    # save series
    def _on_save_series(self, _evt):
        if self._nus is None or self._eval_fn is None:
            return

        with wx.DirDialog(self, "Choose output folder") as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            base_outdir = dlg.GetPath()

        default_name = self._proposed_folder_name()
        with wx.TextEntryDialog(self,
                                "Subfolder name (edit if needed):",
                                "Save series",
                                value=default_name) as ed:
            if ed.ShowModal() != wx.ID_OK:
                return
            subname = ed.GetValue().strip() or default_name

        outdir = os.path.join(base_outdir, self._slug(subname))
        os.makedirs(outdir, exist_ok=True)

        old_idx = self._idx
        try:
            for k, nu in enumerate(self._nus):
                self._idx = k
                self._evaluate_and_plot()
                fname = os.path.join(
                    outdir,
                    f"orbit_{int(nu):03d}.{int(round((nu % 1.0) * 100)):02d}.png"  # filenames unchanged
                )
                self.fig.savefig(fname, dpi=150, bbox_inches="tight")
        finally:
            self._idx = old_idx
            self._evaluate_and_plot()

        wx.MessageBox(f"Saved series to:\n{outdir}",
                    "Export complete", wx.OK | wx.ICON_INFORMATION, self)
    
    def _on_save_orbit(self, _evt=None):
        """Save the current orbital plot as displayed."""
        
        op = getattr(self, "orb_panel", None)
        if not op:
            wx.MessageBox("Orbital plot panel not found.", "Save orbit",
                        wx.OK | wx.ICON_WARNING, self)
            return

        with wx.FileDialog(
            self,
            message="Save orbital plot",
            wildcard="PNG (*.png)|*.png|PDF (*.pdf)|*.pdf|SVG (*.svg)|*.svg|JPEG (*.jpg;*.jpeg)|*.jpg;*.jpeg",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
            defaultFile="orbit_plot.png",
        ) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            path = dlg.GetPath()

        fig = getattr(op, "figure", None) or getattr(op, "fig", None)
        if fig is None:
            wx.MessageBox("Could not find the orbital figure to save.", "Save orbit",
                        wx.OK | wx.ICON_WARNING, self)
            return

        # --- Temporarily switch to floating legend for saved figure ---
        old_side_flag = getattr(op, "_use_side_legend", False)

        try:
            op._use_side_legend = False
            op._update_legend()

            # Make sure the orbit canvas is rendered with the mpl legend present
            try:
                if hasattr(op, "canvas"):
                    op.canvas.draw()
            except Exception:
                pass

            ext = os.path.splitext(path)[1].lower()
            save_kwargs = {"bbox_inches": "tight"}
            if ext in (".png", ".jpg", ".jpeg"):
                save_kwargs["dpi"] = 200

            # Include legend if it lives outside axes
            extra = []
            leg = getattr(op, "_legend_obj", None)
            if leg is not None:
                extra.append(leg)
            if extra:
                save_kwargs["bbox_extra_artists"] = extra

            fig.savefig(path, **save_kwargs)

        except Exception as e:
            wx.MessageBox(f"Could not save:\n{e}", "Save orbit",
                        wx.OK | wx.ICON_ERROR, self)
            return

        finally:
            # Restore live UI mode
            op._use_side_legend = old_side_flag
            try:
                op._update_legend()
            except Exception:
                pass
            try:
                if hasattr(self, "_update_orbit_side_legend"):
                    self._update_orbit_side_legend()
            except Exception:
                pass
            try:
                if hasattr(op, "canvas"):
                    op.canvas.draw_idle()
            except Exception:
                pass

        wx.MessageBox(f"Saved orbital plot to:\n{path}", "Save orbit",
                    wx.OK | wx.ICON_INFORMATION, self)

    def _enc_by_id(self, eid) -> Optional[dict]:
        if eid is None:
            return None
        for rec in getattr(self, "_encounters", []):
            if str(rec.get("id")) == str(eid):
                return rec
        return None

    def _nu_for_encounter(self, enc_or_id) -> Optional[float]:
        """Kept for back-compat; returns ν if available (and stores it if resolved)."""
        enc = enc_or_id if isinstance(enc_or_id, dict) else self._enc_by_id(enc_or_id)
        if not enc:
            return None

        if enc.get("nu_deg") is not None:
            try:
                return float(enc["nu_deg"]) % 360.0
            except Exception:
                return None

        utc = enc.get("utc_iso") or enc.get("utc")
        resolver = getattr(self, "_nu_resolver", None)
        if utc and callable(resolver):
            try:
                nu = resolver(str(utc))
                if nu is None:
                    return None
                nu = float(nu) % 360.0

                #persist back onto the encounter so future calls don't need Horizons
                enc["nu_deg"] = nu
                enc["true_anom_deg"] = nu

                return nu
            except Exception:
                return None

        return None


    def _M_for_encounter(self, enc_or_id) -> Optional[float]:
        """
        Returns MEAN anomaly M (deg).
        Priority:
        1) Stored M in encounter dict
        2) Resolve ν (true anomaly) via stored ν or self._nu_resolver(utc_iso)
        3) Convert ν -> M using current eccentricity
        Persists resolved ν and computed M back onto the encounter dict.
        """
        enc = enc_or_id if isinstance(enc_or_id, dict) else self._enc_by_id(enc_or_id)
        if not enc:
            return None

        #plume cache (only if you actually have stable keys)
        # If you key only by encounter_id, replace the tuple logic accordingly.
        try:
            if enc.get("encounter_tag") == "plume":
                enc_id = enc.get("encounter_id") or enc.get("id")
                cache = getattr(self, "_plume_M_by_encounter", None)
                if isinstance(cache, dict) and enc_id:
                    cached = cache.get(str(enc_id))
                    if cached is not None:
                        M = float(cached) % 360.0
                        enc["mean_anom_deg"] = M
                        enc["M_deg"] = M
                        return M
        except Exception:
            pass


        #Prefer stored M
        for k in ("mean_anom_deg", "M_deg"):
            v = enc.get(k)
            if v is not None:
                try:
                    M = float(v) % 360.0
                    if np.isfinite(M):
                        # normalize storage
                        enc["mean_anom_deg"] = M
                        enc["M_deg"] = M
                        return M
                except Exception:
                    pass

        #Resolve ν using existing helper 
        nu = self._nu_for_encounter(enc)
        if nu is None:
            return None

        try:
            nu = float(nu) % 360.0
        except Exception:
            return None

        # Persist ν for future calls
        enc["nu_deg"] = nu
        enc["true_anom_deg"] = nu

        #Convert ν -> M
        try:
            M = float(_true_to_mean_anomaly_deg(nu, self._ecc_safe())) % 360.0
        except Exception:
            return None

        if not np.isfinite(M):
            return None

        # Persist M
        enc["mean_anom_deg"] = M
        enc["M_deg"] = M
        enc["phase_src"] = enc.get("phase_src") or ("horizons:nu->M" if callable(getattr(self, "_nu_resolver", None)) else "nu->M")

        try:
            if enc.get("encounter_tag") == "plume":
                enc_id = enc.get("encounter_id") or enc.get("id")
                cache = getattr(self, "_plume_M_by_encounter", None)
                if isinstance(cache, dict) and enc_id:
                    cache[str(enc_id)] = float(M)
        except Exception:
            pass


        return M


    def _on_plot_encounter(self, _evt=None):
        enc = self._selected_encounter()
        if not enc:
            return

        M = self._M_for_encounter(enc)
        if M is None or not np.isfinite(M):
            wx.MessageBox(
                "Could not resolve mean anomaly (M) for this encounter.",
                "No M",
                wx.OK | wx.ICON_WARNING
            )
            return

        # Prefer:
        # - for plume rows: the dropdown "name" (e.g., "HST 2014-01-22T...")
        # - for spacecraft/Clipper: the concise id (e.g., "Clipper-E01")
        tag = str(enc.get("encounter_tag") or enc.get("kind") or "").lower()
        if tag == "plume":
            label = str(enc.get("name") or enc.get("id") or enc.get("encounter_id") or "Plume").strip()
        else:
            label = str(enc.get("id") or enc.get("name") or "Encounter").strip()

        op = getattr(self, "orb_panel", None)
        if op and hasattr(op, "set_secondary_M"):
            op.set_secondary_M(float(M), label=label)




    def _on_goto_encounter(self, _evt=None):
        eid = None
        try:
            if callable(self.get_encounter_id):
                eid = self.get_encounter_id()
        except Exception:
            pass
        enc = self._enc_by_id(eid) or self._selected_encounter()
        if not enc:
            return

        M = self._M_for_encounter(enc)
        if M is None or not np.isfinite(M):
            return

        goto = getattr(self, "_goto_nu", None)
        if callable(goto):
            try:
                goto(float(M))
                return
            except Exception:
                pass

        if hasattr(self, "_draw_orbit_position"):
            try:
                self._draw_orbit_position(float(M))
            except Exception:
                pass

    def _refresh_encounter_choices(self, select_id: str | None = None):
        """
        Populate cmb_enc from controller-provided encounters PLUS locally-loaded plumes + JUICE.

        Key behavior:
        - Uses self._encounters as the base (Clipper) list.
        - Appends JUICE encounters (loaded on demand).
        - Appends plume encounters (loaded on demand).
        - De-dupes by encounter "id".
        - Rewrites self._encounters to the combined list so _selected_encounter() stays consistent.
        - Preserves selection if possible.
        """
        encs_base = list(getattr(self, "_encounters", []) or [])

        # --- JUICE ---
        juice_encs: list[dict] = []
        try:
            juice_encs = list(self._get_juice_encounters_for_list() or [])
        except Exception as e:
            print(f"[ScalarPanel] JUICE list build failed: {e}")
            juice_encs = []

        # --- Plumes ---
        plume_encs: list[dict] = []
        try:
            plume_encs = list(self._get_plume_encounters_for_list() or [])
        except Exception as e:
            print(f"[ScalarPanel] plume list build failed: {e}")
            plume_encs = []

        # Debug (keep for now)
        print(f"[ScalarPanel] refresh: base={len(encs_base)} juice={len(juice_encs)} plumes={len(plume_encs)}")

        # Combine (choose ordering you want)
        encs_all = encs_base + juice_encs + plume_encs

        # De-dupe by stable id
        seen = set()
        dedup: list[dict] = []
        for e in encs_all:
            if not isinstance(e, dict):
                continue
            eid = str(e.get("id") or e.get("name") or "").strip()
            if not eid or eid in seen:
                continue
            seen.add(eid)
            # ensure we have both id/name
            if "id" not in e:
                e["id"] = eid
            if "name" not in e:
                e["name"] = eid
            dedup.append(e)

        encs = dedup

        # Store back so selection lookup uses the same list that populates the dropdown
        self._encounters = encs

        # Sort combined list before showing
        try:
            encs.sort(key=self._enc_sort_key)
        except Exception:
            pass

        self._encounters = encs

        labels = [self._enc_display_label(e) for e in encs]

        cmb = getattr(self, "cmb_enc", None)
        if not cmb:
            return

        # Preserve current selection label if any
        try:
            cur_label = cmb.GetStringSelection() if hasattr(cmb, "GetStringSelection") else None
        except Exception:
            cur_label = None

        cmb.Freeze()
        try:
            # Populate items
            if hasattr(cmb, "SetItems"):
                cmb.SetItems(labels)
            else:
                cmb.Clear()
                if labels:
                    if hasattr(cmb, "AppendItems"):
                        cmb.AppendItems(labels)
                    else:
                        for lbl in labels:
                            cmb.Append(lbl)

            # Choose selection index
            if not labels:
                if hasattr(cmb, "SetValue"):
                    cmb.SetValue("")
                return

            target_idx = None

            # 1) explicit select_id by id match
            if select_id is not None:
                sid = str(select_id)
                for i, e in enumerate(encs):
                    if str(e.get("id")) == sid:
                        target_idx = i
                        break

            # 2) keep current label selection if still present
            if target_idx is None and cur_label and cur_label in labels:
                target_idx = labels.index(cur_label)

            # 3) fallback: keep _selected_encounter_id if it matches
            if target_idx is None:
                sel_id = getattr(self, "_selected_encounter_id", None)
                if sel_id is not None:
                    sid = str(sel_id)
                    for i, e in enumerate(encs):
                        if str(e.get("id")) == sid:
                            target_idx = i
                            break

            # 4) default to first
            if target_idx is None:
                target_idx = 0

            cnt = cmb.GetCount() if hasattr(cmb, "GetCount") else len(labels)
            if 0 <= target_idx < cnt:
                try:
                    if hasattr(cmb, "SetStringSelection"):
                        if not cmb.SetStringSelection(labels[target_idx]):
                            cmb.SetSelection(target_idx)
                    else:
                        cmb.SetSelection(target_idx)
                except Exception:
                    pass

            # Keep internal selected id in sync
            try:
                self._selected_encounter_id = str(encs[int(target_idx)].get("id") or labels[int(target_idx)])
            except Exception:
                pass

        finally:
            cmb.Thaw()


    def _selected_encounter(self) -> Optional[dict]:
        """Return the encounter currently selected in self.cmb_enc, or None."""
        cmb = getattr(self, "cmb_enc", None)
        encs = getattr(self, "_encounters", []) or []
        if cmb is None or not encs:
            return None

        # Prefer index — we kept list order in _refresh_encounter_choices
        try:
            idx = cmb.GetSelection()
            if idx is not None and idx != -1:
                try:
                    return encs[int(idx)]
                except Exception:
                    pass
        except Exception:
            pass

        # Fallback by label match
        try:
            label = cmb.GetStringSelection()
            if label:
                for e in encs:
                    if str(e.get("id") or e.get("name") or "") == str(label):
                        return e
        except Exception:
            pass

        return None

    def _get_plume_encounters_for_list(self) -> list[dict]:
        self._ensure_plumes_loaded_into_dropdown()
        plume_meta = getattr(self, "_plume_meta_by_id", {}) or {}

        out: list[dict] = []
        for enc_id, meta in plume_meta.items():
            obs = meta.get("observer") or "UNKNOWN"
            utc = meta.get("utc_iso") or meta.get("utc") or ""

            disp_id = enc_id.replace("PLUME_", "", 1) if str(enc_id).startswith("PLUME_") else str(enc_id)

            det = str(meta.get("detection") or "").strip().upper()
            if det == "Y":
                det_label = "Y"
            elif det == "N":
                det_label = "N"
            else:
                det_label = "U"

            label = f"{obs} {utc} [{det_label}]".strip() or disp_id

            M = meta.get("mean_anom_deg")
            if M in ("", None):
                M = meta.get("M_deg")

            nu = meta.get("true_anom_deg")
            if nu in ("", None):
                nu = meta.get("nu_deg")

            ecc = meta.get("eccentricity")
            if ecc in ("", None):
                ecc = meta.get("e")

            out.append({
                "id": enc_id,
                "name": label,
                "encounter_id": enc_id,
                "encounter_tag": "plume",

                "utc_iso": utc,
                "observer": obs,
                "exposure_s": meta.get("exposure_s"),

                "mean_anom_deg": M,
                "M_deg": M,
                "true_anom_deg": nu,
                "nu_deg": nu,
                "eccentricity": ecc,
                "phase_src": meta.get("phase_src") or "file",
                "detection": det,
            })

        return out

    def _get_juice_encounters_for_list(self) -> list[dict]:
        self._ensure_juice_loaded_into_dropdown()
        juice_meta = getattr(self, "_juice_meta_by_id", {}) or {}

        out: list[dict] = []
        for enc_id, meta in juice_meta.items():
            # enc_id is JUICE-E1, JUICE-E2, ...

            obs = meta.get("observer") or "UNKNOWN"
            utc = meta.get("utc_iso") or meta.get("utc") or ""

            # For JUICE you said: "same naming convention as clipper encounters (i.e JUICE-E1)"
            # So default dropdown label should just be the encounter id.
            # If you later decide you want a richer label, swap to: f"{enc_id} {utc}".strip()
            label = str(enc_id)

            # --- pull phase fields from meta (as strings -> floats later in _M_for_encounter) ---
            M = meta.get("mean_anom_deg")
            if M in ("", None):
                M = meta.get("M_deg")

            nu = meta.get("true_anom_deg")
            if nu in ("", None):
                nu = meta.get("nu_deg")

            ecc = meta.get("eccentricity")
            if ecc in ("", None):
                ecc = meta.get("e")

            out.append({
                # identifiers
                "id": enc_id,                  # what selection returns
                "name": label,                 # what dropdown shows
                "encounter_id": enc_id,         # stable key for caches/lookup
                "encounter_tag": "juice",

                # timing/meta
                "utc_iso": utc,
                "observer": obs,

                # phase/ephemeris (propagate through to _M_for_encounter)
                "mean_anom_deg": M,
                "M_deg": M,
                "true_anom_deg": nu,
                "nu_deg": nu,
                "eccentricity": ecc,
                "phase_src": meta.get("phase_src") or "file_true_to_mean",
            })

        return out

    def set_encounters(
        self,
        encounters: List[Dict[str, Any]],
        select_id: Optional[str] = None,
    ) -> None:
        self._encounters = list(encounters or [])
        self._refresh_encounter_choices(select_id=select_id)

    def set_encounters_from_ids(
        self,
        ids: List[str],
        lookup: Callable[[str], Dict[str, Any]],
        label_builder: Optional[Callable[[str, Dict[str, Any]], str]] = None,
        select_id: Optional[str] = None,
    ) -> None:
        encs: list[dict] = []
        for eid in ids or []:
            rec = (lookup(eid) or {}).copy()
            if not rec:
                continue
            name = label_builder(eid, rec) if label_builder else str(eid)
            encs.append({
                "id": eid,
                "name": name,
                "nu_deg": rec.get("nu_deg") or rec.get("true_anom_deg"),
                "utc_iso": rec.get("utc_iso") or rec.get("utc"),
                "lat_deg": rec.get("lat_deg") or rec.get("lat"),
                "lon_deg": rec.get("lon_deg") or rec.get("lon"),
                # if your control panel stored mean_anom_deg, Point panel can display it; we resolve M here anyway
            })
        self.set_encounters(encs, select_id=select_id)


    # (mini) orbit panel helpers (unchanged draw, label made neutral)
    def _init_orbit_axes(self):
        ax = self.ax_orbit
        ax.clear()
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)

        theta = np.linspace(0, 2*np.pi, 361)
        ax.plot(np.cos(theta), np.sin(theta), linewidth=1.0, zorder=1)
        ax.plot([0.0, 1.0], [0.0, 0.0], linewidth=1.0, alpha=0.3, zorder=1)
        ax.scatter([1.0], [0.0], s=18, zorder=2)
        ax.text(1.02, -0.06, "0°", fontsize=8, zorder=2)  # neutral (no ν/M wording)

        (self._nu_line,) = ax.plot([0,0],[0,0], linewidth=2.0, zorder=4)
        self._nu_point   = ax.scatter([0],[0], s=28, zorder=5)
        self._nu_text    = ax.text(0, 0, "", fontsize=8, zorder=6)
        for a in (self._nu_line, self._nu_point, self._nu_text):
            a.set_visible(False)

        # Secondary (encounter) marker in gray — independent of active (blue) marker
        self._enc_sel_scatter = ax.scatter([], [], s=34, zorder=4,
                                        color="#7f7f7f", alpha=0.95, edgecolors="none")
        (self._enc_sel_line,) = ax.plot([], [], linewidth=1.6, zorder=4,
                                        color="#7f7f7f", alpha=0.8)


        self._enc_data   = {}
        self._enc_texts  = {}
        self._enc_scatter = ax.scatter([], [], s=24, zorder=3)
        self._enc_lines   = LineCollection([], linewidths=1.5, zorder=3)
        ax.add_collection(self._enc_lines)

        self.fig_orbit.tight_layout()
        self.canvas_orbit.draw()

    def _draw_orbit_position(self, deg_val):
        """Draw active position (treat input as M)."""
        if not hasattr(self, "_nu_line"): self._init_orbit_axes()
        if deg_val is None:
            for a in (self._nu_line, self._nu_point, self._nu_text): a.set_visible(False)
            self.canvas_orbit.draw_idle(); return
        th = np.deg2rad(float(deg_val) % 360.0); x, y = np.cos(th), np.sin(th)
        self._nu_line.set_data([0.0, x], [0.0, y])
        self._nu_point.set_offsets([[x, y]])
        self._nu_text.set_position((x*1.05, y*1.05)); self._nu_text.set_text(f"{deg_val:.0f}°")
        for a in (self._nu_line, self._nu_point, self._nu_text): a.set_visible(True)
        self.canvas_orbit.draw_idle()

    def set_orbit_encounter_marker(self, key, nu_deg, label=None, on=True):
        """Kept for API parity; you may pass M here as well."""
        if not hasattr(self, "_enc_data"): self._init_orbit_axes()
        if on:
            th = np.deg2rad(float(nu_deg) % 360.0)
            self._enc_data[key] = {"x": np.cos(th), "y": np.sin(th), "label": label or ""}
        else:
            self._enc_data.pop(key, None)
        self._redraw_orbit_encounters()

    def _set_secondary_marker(self, deg_val, label=None):
        op = getattr(self, "orb_panel", None)
        if op and hasattr(op, "set_secondary_M"):
            op.set_secondary_M(None if deg_val is None else float(deg_val), label=label)



    def clear_secondary_marker(self):
        op = getattr(self, "orb_panel", None)
        if not op:
            return
        clearer = getattr(op, "clear_secondary", None) or getattr(op, "set_secondary_M", None)
        if callable(clearer):
            clearer(None)



    def _redraw_orbit_encounters(self):
        ax = self.ax_orbit
        xs, ys, segs = [], [], []
        for key, d in self._enc_data.items():
            x, y = d["x"], d["y"]
            xs.append(x); ys.append(y); segs.append([(0,0), (x,y)])
            txt = self._enc_texts.get(key)
            if txt is None:
                self._enc_texts[key] = ax.text(x*1.05, y*1.05, d["label"], fontsize=8, zorder=6)
            else:
                txt.set_position((x*1.05, y*1.05)); txt.set_text(d["label"])

        for key in list(self._enc_texts.keys()):
            if key not in self._enc_data:
                try: self._enc_texts[key].remove()
                except Exception: pass
                del self._enc_texts[key]

        if xs:
            self._enc_scatter.set_offsets(np.c_[xs, ys])
            self._enc_lines.set_segments(segs)
        else:
            self._enc_scatter.set_offsets(np.empty((0,2)))
            self._enc_lines.set_segments([])

        self.canvas_orbit.draw_idle()

    # ---------------- fixed-range API (kPa) ----------------
    def set_fixed_stress_range_kpa(self, vmin=-100.0, vmax=100.0):
        self._fixed_clim_kpa = (float(vmin), float(vmax))
        if hasattr(self, "chk_auto"):
            self.chk_auto.SetValue(False)
        self._auto_range = False
        if hasattr(self, "sp_l"): self.sp_l.SetValue(int(vmin))
        if hasattr(self, "sp_u"): self.sp_u.SetValue(int(vmax))
        self._apply_fixed_clim()
        try:
            self._replot()
        except Exception:
            pass

    def _apply_fixed_clim(self):
        if self._im is None:
            return
        if self._fixed_clim_kpa is None:
            try:
                self._im.autoscale()
            except Exception:
                pass
        else:
            vmin_kpa, vmax_kpa = self._fixed_clim_kpa
            self._im.set_clim(float(vmin_kpa), float(vmax_kpa))
        try:
            self.canvas.draw_idle()
        except Exception:
            pass

    def clear_fixed_stress_range(self):
        self._fixed_clim_kpa = None
        self._replot()

    def set_vector_density(self, nlat: int, nlon: int):
        self.vec_n_lat = max(2, int(nlat))
        self.vec_n_lon = max(2, int(nlon))
        if hasattr(self, "sp_vlat"): self.sp_vlat.SetValue(self.vec_n_lat)
        if hasattr(self, "sp_vlon"): self.sp_vlon.SetValue(self.vec_n_lon)
        self._replot()


def get_or_create_scalar_popup(
    parent,
    get_system_id: Optional[Callable[[], str]] = None,
    get_encounter_id: Optional[Callable[[], Optional[str]]] = None,
) -> Tuple[wx.Frame, "ScalarPlotPanel"]:
    """
    Build the ScalarPlotPanel popup safely.
    """

    def _maybe(chain: Tuple[str, ...]):
        obj = parent
        for name in chain:
            obj = getattr(obj, name, None)
            if obj is None:
                return None
        return obj if callable(obj) else None

    if get_system_id is None:
        get_system_id = _maybe(("get_system_id",)) or _maybe(("global_params_panel", "get_system_id"))
    if get_encounter_id is None:
        get_encounter_id = _maybe(("get_selected_encounter_id",)) or _maybe(("point_panel", "get_selected_encounter_id"))
    get_eccentricity = _maybe(("get_eccentricity",)) or _maybe(("global_params_panel", "get_eccentricity"))

    frame = wx.Frame(parent, title="StressViz Map", size=(1100, 780))
    panel = ScalarPlotPanel(
        frame,
        get_system_id=get_system_id,
        get_encounter_id=get_encounter_id,
        get_eccentricity=(get_eccentricity or (lambda: 0.0)),
    )

    sizer = wx.BoxSizer(wx.VERTICAL)
    sizer.Add(panel, 1, wx.EXPAND)
    frame.SetSizer(sizer)
    frame.SetClientSize((panel.FromDIP(1250), panel.FromDIP(860)))  # width, height
    frame.SetMinSize((panel.FromDIP(1100), panel.FromDIP(780)))
    frame.Layout()

    try:
        parent.scalar_panel_ref = panel
    except Exception:
        pass

    try:
        host = wx.GetTopLevelParent(parent)
        if host is not None:
            host.scalar_panel_ref = panel
    except Exception:
        pass

    try:
        frame.scalar_panel_ref = panel
    except Exception:
        pass

    frame.Centre()
    frame.Show()

    def _deferred_init():
        # Hook orbit "goto" function
        try:
            orbit = getattr(parent, "orbit_panel", None)
            if orbit is not None:
                goto = getattr(orbit, "highlight_active", None) or getattr(orbit, "_draw_orbit_position", None)
                if callable(goto):
                    panel._goto_nu = lambda deg: goto(float(deg))
            if not callable(getattr(panel, "_goto_nu", None)) and hasattr(panel, "_draw_orbit_position"):
                panel._goto_nu = lambda deg: panel._draw_orbit_position(float(deg))
        except Exception:
            pass

        # True anomaly resolver (Horizons) from parent/control panel (optional)
        try:
            resolver = None

            # 1) explicit hook if parent already exposes it
            resolver = getattr(parent, "_resolve_true_anomaly", None)
            if callable(resolver):
                panel._nu_resolver = resolver
            else:
                # 2) common locations: control_panel or analysis_control_panel
                cp = getattr(parent, "control_panel", None) or getattr(parent, "analysis_control_panel", None)
                if cp is not None:
                    # your ACP method is currently named "_"
                    cand = getattr(cp, "_", None)
                    if callable(cand):
                        panel._nu_resolver = cand
        except Exception:
            pass


        # Populate encounters and sync selection
        try:
            enc_map: Dict[str, Dict[str, Any]] = getattr(parent, "encounters_by_id", {}) or {}
            ids: List[str] = sorted(enc_map.keys(), key=lambda x: str(x))

            current_sel = None
            get_sel = _maybe(("get_selected_encounter_id",)) or _maybe(("point_panel", "get_selected_encounter_id"))
            if callable(get_sel):
                try:
                    current_sel = get_sel()
                except Exception:
                    current_sel = None

            if hasattr(panel, "set_encounters_from_ids"):
                try:
                    panel.set_encounters_from_ids(
                        ids,
                        lookup=lambda eid: enc_map.get(eid),
                        label_builder=lambda eid, _rec: str(eid),
                        select_id=(current_sel or (ids[0] if ids else None)),
                    )
                except Exception:
                    pass
            else:
                encs: List[Dict[str, Any]] = []
                for eid in ids:
                    rec = enc_map.get(eid, {}) or {}
                    encs.append({
                        "id": eid,
                        "name": str(eid),
                        "nu_deg": rec.get("nu_deg") or rec.get("true_anom_deg"),
                        "utc_iso": rec.get("utc_iso") or rec.get("utc"),
                        "lat_deg": rec.get("lat_deg") or rec.get("lat"),
                        "lon_deg": rec.get("lon_deg") or rec.get("lon"),
                    })
                try:
                    panel._encounters = encs
                    if hasattr(panel, "_refresh_encounter_choices"):
                        panel._refresh_encounter_choices()
                    if current_sel and hasattr(panel, "cmb_enc"):
                        try:
                            idx = next(i for i, e in enumerate(encs) if e.get("id") == current_sel)
                            panel.cmb_enc.SetSelection(idx)
                        except StopIteration:
                            pass
                except Exception:
                    pass
        except Exception:
            pass

        # --- Snap initial position to the selected encounter's exact M by injecting ν0 ---
        try:
            # Find selected encounter id
            sel_id = None
            get_sel = _maybe(("get_selected_encounter_id",)) or _maybe(("point_panel", "get_selected_encounter_id"))
            if callable(get_sel):
                try:
                    sel_id = get_sel()
                except Exception:
                    sel_id = None
            if sel_id is None and getattr(panel, "cmb_enc", None):
                try:
                    sel_txt = panel.cmb_enc.GetStringSelection()
                    if sel_txt:
                        sel_id = sel_txt
                except Exception:
                    pass

            # Resolve target M
            M0 = None
            if sel_id is not None:
                try:
                    M0 = panel._M_for_encounter(sel_id)
                except Exception:
                    M0 = None

            if M0 is not None and np.isfinite(M0):
                def _inject_and_snap():
                    try:
                        # Wait until series is bound (nus + eval_fn ready)
                        if not (panel._eval_fn and panel._nus):
                            wx.CallLater(50, _inject_and_snap)
                            return

                        # Compute e and exact ν0 that maps to M0
                        e = panel._ecc_safe()
                        nu0 = _mean_to_true_anomaly_deg(M0, e)  # exact inverse

                        # If ν0 already present near some step, just snap to it
                        Ms = [float(_true_to_mean_anomaly_deg(nu, e)) % 360.0 for nu in panel._nus]
                        diffs = [abs(((m - M0) + 180.0) % 360.0 - 180.0) for m in Ms]
                        k = int(np.argmin(diffs))
                        if diffs[k] <= (360.0 / max(1, len(panel._nus))) * 0.25:
                            # Close enough—just snap
                            panel._idx = k
                            try: panel.sld_orbit.SetValue(k)
                            except Exception: pass
                            panel._evaluate_and_plot()
                            panel._startup_target_M = float(M0)   
                            panel._init_snap_done = True      
                            return

                        # Otherwise, inject ν0 as a first step so we hit M0 exactly
                        panel._nus = [float(nu0)] + list(panel._nus)
                        panel._idx = 0
                        try:
                            panel.sld_orbit.SetRange(0, max(0, len(panel._nus) - 1))
                            panel.sld_orbit.SetValue(0)
                        except Exception:
                            pass
                        panel._evaluate_and_plot()
                        panel._startup_target_M = float(M0)
                        panel._init_snap_done = True

                    except Exception:
                        pass

                _inject_and_snap()
        except Exception:
            pass


    wx.CallLater(1, _deferred_init)

    return frame, panel









