# stressviz/orbital_panel.py
import wx
import pandas as pd
import numpy as np
import matplotlib

# Use the wx backend
matplotlib.use("WXAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib import rcParams


class OrbitalPositionPanel(wx.Panel):
    """
    Polar plot showing orbital position ν (true anomaly, deg).
    Optional encounter data -> checklist + plotted points.

    Public API:
      - set_nu(nu_deg: float|None)
      - set_title(text: str)
      - set_secondary_M(deg_or_none, label: str|None = None)
      - clear_secondary()
      - set_legend_position(where: str)
          where ∈ {'outside-right','top-left','top-right','bottom-left','bottom-right'}
          default = 'outside-right' (floating in figure space; never overlaps labels)
    """

    def __init__(self, parent, datafile=None, title="Orbital position (M)"):
        super().__init__(parent)

        # -------- Optional data load --------
        self.df = None
        self.has_data = False
        if datafile:
            try:
                self.df = self.load_data(datafile)
                self.has_data = self.df is not None and not self.df.empty
            except Exception:
                self.df = None
                self.has_data = False

        self.selected_encounters = []
        if self.has_data:
            self.selected_encounters = list(sorted(self.df["encounter"].unique()))

        # -------- Figure / Axes --------
        self.figure = Figure(figsize=(3.0, 3.0))
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.canvas.SetMinSize((220, 220))
        self.ax = self.figure.add_subplot(111, polar=True)

        # Polar styling: 0° at South, CCW increasing
        self.ax.set_theta_zero_location("S")
        self.ax.set_theta_direction(1)
        self._configure_axes()

        # Legend settings (floating by default)
        self._legend_obj = None
        self._legend_where = "outside-right"

        # ν marker state
        self._nu_marker = None
        self._nu_deg = None
        self._show_nu_marker = True

        # Secondary (gray) marker state
        self._secondary_markers = [] # (deg_float, label_str)
        self._secondary_palette = rcParams["axes.prop_cycle"].by_key().get("color", []) or ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"]
        self._secondary_palette_i = 0
        self._secondary_color_by_key = {}

        # -------- Layout --------
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if self.has_data:
            left = wx.BoxSizer(wx.VERTICAL)
            left.Add(wx.StaticText(self, label="Select Encounters:"), 0, wx.ALL, 5)

            self.select_all_btn = wx.Button(self, label="Select All")
            self.deselect_all_btn = wx.Button(self, label="Deselect All")
            self.select_all_btn.Bind(wx.EVT_BUTTON, self.on_select_all)
            self.deselect_all_btn.Bind(wx.EVT_BUTTON, self.on_deselect_all)
            left.Add(self.select_all_btn, 0, wx.ALL | wx.EXPAND, 5)
            left.Add(self.deselect_all_btn, 0, wx.ALL | wx.EXPAND, 5)

            # NEW: encounter rows with checkbox + color picker
            if not hasattr(self, "_enc_enabled_by_id"):
                self._enc_enabled_by_id = {}
            if not hasattr(self, "_enc_color_by_id"):
                self._enc_color_by_id = {}

            self.encounter_rows = self._build_encounter_rows(self, self.selected_encounters)
            left.Add(self.encounter_rows, 1, wx.EXPAND | wx.ALL, 5)

            main_sizer.Add(left, 0, wx.EXPAND)

        main_sizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(main_sizer)

        # -------- Initial draw --------
        self.set_title(title)
        self._draw_unit_circle()
        if self.has_data:
            self.plot_selected()
        else:
            self._update_legend()
            self.canvas.draw_idle()

    # ================= Data I/O (optional) =================
    def load_data(self, filepath):
        """Load encounters from an Excel file with expected columns."""
        df = pd.read_excel(filepath)
        df = df.rename(
            columns={
                "Time of C/A (TCA) (ET/TDB)": "et_seconds",
                "True Anomaly of Body (deg)": "true_anomaly_deg",
                "Encounter Tag": "encounter",
            }
        )
        if "true_anomaly_deg" not in df or "encounter" not in df:
            return None
        df["true_anomaly_deg"] = df["true_anomaly_deg"] % 360.0
        df["true_anomaly_rad"] = np.radians(df["true_anomaly_deg"])
        df["r"] = 1.0  # unit circle radius
        return df

    # ================= Public API =================
    def set_title(self, text):
        try:
            self.ax.set_title(str(text), y=1.08, fontsize=10)
        except Exception:
            pass

    def set_nu(self, nu_deg):
        self._nu_deg = None if nu_deg is None else float(nu_deg) % 360.0

        self._configure_axes()
        self._draw_unit_circle()

        if self.has_data and self.selected_encounters:
            self._plot_encounters_only()

        self._draw_secondary_markers()

        if self._nu_marker is not None:
            try:
                self._nu_marker.remove()
            except Exception:
                pass
            self._nu_marker = None

        self._draw_nu_marker()

        self._update_legend()
        self.canvas.draw_idle()

    def _draw_nu_marker(self):
        """Draw the black ν marker at self._nu_deg (degrees)."""
        if not getattr(self, "_show_nu_marker", True):
            return
        if self._nu_deg is None:
            return

        try:
            th = np.radians(float(self._nu_deg))
        except Exception:
            return

        self._nu_marker = self.ax.scatter(
            [th], [1.0],
            s=110,                 # slightly larger so "+" is visible
            marker="+",
            color="#000000",
            alpha=0.95,
            linewidths=1.8,       # makes the plus thicker
            zorder=5,
            label=None,
        )

    def set_show_nu_marker(self, show: bool):
        self._show_nu_marker = bool(show)
        # use whatever redraw function you have; if you added _refresh_plot, call that
        if hasattr(self, "_refresh_plot"):
            self._refresh_plot()
        else:
            # fallback: just redraw canvas; but note this won't re-add artists if ax was cleared elsewhere
            try:
                self.canvas.draw_idle()
            except Exception:
                pass


    def _color_for_secondary_key(self, key: str) -> str:
        key = str(key or "").strip() or "Encounter"

        if key in self._secondary_color_by_key:
            return self._secondary_color_by_key[key]

        pal = self._secondary_palette or ["#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
        # avoid grey-like colors if they appear
        while True:
            c = pal[self._secondary_palette_i % len(pal)]
            self._secondary_palette_i += 1
            if str(c).lower() not in ("#7f7f7f", "0.5", "gray", "grey", "#000000"):
                break

        self._secondary_color_by_key[key] = c
        return c
    

    def set_secondary_M(self, M: float, label: str = ""):
        # keep backward compatibility
        try:
            deg = float(M)
        except Exception:
            self._secondary_markers = []
            self._refresh_plot()
            return
        self._secondary_markers = [(deg, str(label or ""))]
        self._refresh_plot()

    def set_legend_position(self, where: str):
        """
        Move legend:
          'outside-right' (default), or 'top-left'|'top-right'|'bottom-left'|'bottom-right'
        """
        self._legend_where = (where or "outside-right").lower()
        self._update_legend()
        self.canvas.draw_idle()

    def _refresh_plot(self):
        """Clear axes and redraw everything in correct order."""
        self.ax.clear()

        # Restore polar config
        self.ax.set_theta_zero_location("S")
        self.ax.set_theta_direction(1)
        self._configure_axes()

        # Base circle
        self._draw_unit_circle()

        # Encounter scatter (if using df)
        if self.has_data:
            self._plot_encounters_only()

        # ν marker (if you have one)
        if hasattr(self, "_draw_nu_marker"):
            self._draw_nu_marker()

        # Secondary markers (NEW multi-support)
        if hasattr(self, "_draw_secondary_markers"):
            self._draw_secondary_markers()

        # Legend
        if hasattr(self, "_update_legend"):
            self._update_legend()

        self.canvas.draw_idle()

    # ================= Internals: axes & drawing =================
    def _configure_axes(self):
        self.ax.set_rmin(0.0)
        self.ax.set_rmax(1.1)
        self.ax.set_ylim(0, 1.1)
        self.ax.set_rticks([])
        self.ax.set_yticklabels([])
        self.ax.yaxis.set_visible(False)
        self.ax.yaxis.grid(False)
        # Hide the polar spine to avoid a second ring
        self.ax.spines["polar"].set_visible(False)
        # Pull theta labels slightly inward so an inside-corner legend can coexist
        try:
            self.ax.set_thetagrids(range(0, 360, 45), frac=0.90)
        except Exception:
            pass

    def _draw_unit_circle(self):
        self.ax.clear()
        self.ax.set_theta_zero_location("S")
        self.ax.set_theta_direction(1)
        self._configure_axes()

        theta = np.linspace(0, 2 * np.pi, 360)
        self.ax.plot(theta, np.ones_like(theta), color="black", linewidth=1, zorder=1)

    def _plot_encounters_only(self):
        if not self.has_data or not self.selected_encounters:
            return

        colors = getattr(self, "_enc_color_by_id", {}) or {}

        for tag in self.selected_encounters:
            subset = self.df[self.df["encounter"] == tag]
            if subset.empty:
                continue

            kw = dict(alpha=0.8, zorder=2, label=tag)

            col = colors.get(str(tag))
            if col:
                kw["color"] = col  # "#RRGGBB"

            self.ax.scatter(subset["true_anomaly_rad"], subset["r"], **kw)

    

    def _draw_secondary_markers(self):
        markers = getattr(self, "_secondary_markers", []) or []
        for item in markers:
            # Backward compatible:
            # old: (deg, label, color)
            # new: (deg, label, color, marker)
            try:
                deg = item[0]
                label = item[1] if len(item) > 1 else None
                color = item[2] if len(item) > 2 else None
                marker = item[3] if len(item) > 3 else "o"
            except Exception:
                continue

            try:
                th = np.radians(float(deg))
            except Exception:
                continue

            self.ax.scatter(
                [th], [1.0],
                s=90,
                marker=marker,
                color=color,
                alpha=0.95,
                edgecolors="none",
                zorder=4,
            )

            '''
            if label:
                self.ax.text(th, 1.05, str(label), fontsize=8, va="center", ha="left", zorder=5)
            '''

    def clear_secondary_markers(self):
        self._secondary_markers = []
        self._refresh_plot()

    def add_secondary_M(self, M: float, label: str = "", color: str = None, marker: str = "o"):
        try:
            deg = float(M)
        except Exception:
            return

        lbl = str(label or "").strip() or "Encounter"

        if not color:
            color = self._color_for_secondary_key(lbl)

        self._secondary_markers.append((deg, lbl, color, marker))

        if hasattr(self, "_refresh_plot"):
            self._refresh_plot()
        else:
            self.canvas.draw_idle()


    def set_secondary_M(self, M: float, label: str = "", color: str = None, marker: str = "o"):
        try:
            deg = float(M)
        except Exception:
            self._secondary_markers = []
            if hasattr(self, "_refresh_plot"):
                self._refresh_plot()
            else:
                self.canvas.draw_idle()
            return

        lbl = str(label or "").strip() or "Encounter"

        if not color:
            color = self._color_for_secondary_key(lbl)

        self._secondary_markers = [(deg, lbl, color, marker)]

        if hasattr(self, "_refresh_plot"):
            self._refresh_plot()
        else:
            self.canvas.draw_idle()

    # ================= Legend (floating) =================

    def _update_legend(self):
        if getattr(self, "_use_side_legend", False):
            if self._legend_obj is not None:
                try:
                    self._legend_obj.remove()
                except Exception:
                    pass
                self._legend_obj = None
            return
        
        if self._legend_obj is not None:
            try:
                self._legend_obj.remove()
            except Exception:
                pass
            self._legend_obj = None

        secondary = getattr(self, "_secondary_markers", []) or []

        has_enc = bool(secondary)
        has_nu = (self._nu_deg is not None) and bool(getattr(self, "_show_nu_marker", True))

        if not (has_enc or has_nu):
            return

        handles = []

        if has_enc:
            # One legend entry per plotted encounter (de-dup by label)
            seen = set()
            for item in secondary:
                try:
                    deg = item[0]
                    lbl = item[1] if len(item) > 1 else None
                    color = item[2] if len(item) > 2 else "#000000"
                    marker = item[3] if len(item) > 3 else "o"
                except Exception:
                    continue

                label = (str(lbl).strip() if lbl is not None else "") or "Encounter"
                if label in seen:
                    continue
                seen.add(label)

                # line-type markers like "+" should use color, not facecolor
                if marker in ["+", "x", "|", "_"]:
                    handles.append(
                        Line2D(
                            [0], [0],
                            marker=marker,
                            linestyle="None",
                            markersize=9,
                            color=color,
                            label=label,
                        )
                    )
                else:
                    handles.append(
                        Line2D(
                            [0], [0],
                            marker=marker,
                            linestyle="None",
                            markersize=8,
                            markerfacecolor=color,
                            markeredgecolor="none",
                            label=label,
                        )
                    )

        if has_nu:
            handles.append(
                Line2D(
                    [0], [0],
                    marker="+",
                    linestyle="None",
                    markersize=10,
                    color="#000000",
                    label=f"M = {self._nu_deg:.2f}°",
                )
            )

        where = (self._legend_where or "outside-right").lower()

        if where == "outside-right":
            try:
                self.canvas.draw()
                renderer = self.canvas.get_renderer()
            except Exception:
                renderer = None

            if renderer is not None:
                tight = self.ax.get_tightbbox(renderer)
                tight_fig = tight.transformed(self.figure.transFigure.inverted())
                x_right = tight_fig.x1
                y_mid = tight_fig.y0 + tight_fig.height / 2.0
            else:
                bbox = self.ax.get_position()
                x_right = bbox.x1
                y_mid = bbox.y0 + bbox.height / 2.0

            gap = 0.010
            x_anchor = min(x_right + gap, 0.985)

            anchor = (x_anchor, y_mid)
            loc = "center left"
            transform = self.figure.transFigure
        else:
            corners = {
                "top-left": ((0.02, 0.98), "upper left"),
                "top-right": ((0.98, 0.98), "upper right"),
                "bottom-left": ((0.02, 0.02), "lower left"),
                "bottom-right": ((0.98, 0.02), "lower right"),
            }
            anchor, loc = corners.get(where, corners["top-left"])
            transform = self.ax.transAxes

        leg = self.ax.legend(
            handles=handles,
            loc=loc,
            bbox_to_anchor=anchor,
            bbox_transform=transform,
            frameon=True,
            framealpha=0.9,
            borderaxespad=0.0,
            handlelength=1.2,
            handletextpad=0.6,
            labelspacing=0.4,
            fontsize=9,
        )
        try:
            leg.set_in_layout(False)
        except Exception:
            pass

        fr = leg.get_frame()
        fr.set_linewidth(0.8)
        fr.set_edgecolor("0.6")
        fr.set_facecolor("white")

        try:
            leg.set_draggable(True)
        except Exception:
            pass

        self._legend_obj = leg


    # ================= Encounter UI (optional) =================
    def _build_encounter_rows(self, parent, encounter_ids):
        """
        Scrollable rows:
        [checkbox] [label] [color picker]
        """
        sc = wx.ScrolledWindow(parent, style=wx.VSCROLL)
        sc.SetScrollRate(0, 12)

        vsz = wx.BoxSizer(wx.VERTICAL)
        self._enc_row_widgets = {}

        for eid in encounter_ids:
            eid = str(eid)

            row = wx.BoxSizer(wx.HORIZONTAL)

            # --- checkbox ---
            chk = wx.CheckBox(sc)
            if eid not in self._enc_enabled_by_id:
                self._enc_enabled_by_id[eid] = False  # change to True if you want default-on
            chk.SetValue(bool(self._enc_enabled_by_id[eid]))

            # --- label ---
            lab = wx.StaticText(sc, label=eid)

            # --- color picker ---
            picker = wx.ColourPickerCtrl(sc)
            # Optional: show an initial default color without forcing plot to use it.
            # If you want "Auto" behavior, only apply color if user has set one (i.e., exists in dict).
            if eid in self._enc_color_by_id:
                try:
                    picker.SetColour(self._enc_color_by_id[eid])  # "#RRGGBB"
                except Exception:
                    pass

            row.Add(chk, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
            row.Add(lab, 1, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
            row.Add(picker, 0, wx.ALIGN_CENTER_VERTICAL)

            vsz.Add(row, 0, wx.EXPAND | wx.ALL, 3)

            self._enc_row_widgets[eid] = (chk, picker, lab)

            # --- bindings ---
            def _on_check(evt, _eid=eid):
                self._enc_enabled_by_id[_eid] = bool(evt.GetEventObject().GetValue())
                self.on_selection_change(None)

            def _on_color(evt, _eid=eid):
                col = evt.GetEventObject().GetColour()
                self._enc_color_by_id[_eid] = col.GetAsString(wx.C2S_HTML_SYNTAX)  # "#RRGGBB"
                # no need to change selected_encounters; just redraw
                self.plot_selected()

            chk.Bind(wx.EVT_CHECKBOX, _on_check)
            picker.Bind(wx.EVT_COLOURPICKER_CHANGED, _on_color)

        sc.SetSizer(vsz)
        vsz.FitInside(sc)
        return sc


    def on_select_all(self, _evt):
        rows = getattr(self, "_enc_row_widgets", {}) or {}
        for eid, (chk, _picker, _lab) in rows.items():
            chk.SetValue(True)
            self._enc_enabled_by_id[str(eid)] = True
        self.on_selection_change(None)


    def on_deselect_all(self, _evt):
        # Disable all encounters
        for eid, (chk, _picker, _lab) in getattr(self, "_enc_row_widgets", {}).items():
            chk.SetValue(False)
            self._enc_enabled_by_id[eid] = False
        self.on_selection_change(None)

    def on_selection_change(self, _evt):
        if not self.has_data:
            return

        # Gather selected encounters from row checkboxes
        rows = getattr(self, "_enc_row_widgets", {}) or {}
        selected = []
        for eid, (chk, _picker, _lab) in rows.items():
            try:
                if chk.GetValue():
                    selected.append(str(eid))
            except Exception:
                pass

        self.selected_encounters = selected
        self.plot_selected()

    def plot_selected(self):
        self._draw_unit_circle()

        if self.has_data and self.selected_encounters:
            self._plot_encounters_only()

        self._draw_secondary_marker()

        if self._nu_deg is not None:
            nu_rad = np.radians(self._nu_deg)
            self._nu_marker = self.ax.scatter(nu_rad, 1.0, marker="o", s=100, zorder=5)

        self._update_legend()
        self.canvas.draw_idle()

    def on_select_all(self, _evt):
        if not self.has_data:
            return
        self.encounter_checklist.SetCheckedItems(range(self.encounter_checklist.GetCount()))
        self.on_selection_change(None)

    def on_deselect_all(self, _evt):
        if not self.has_data:
            return
        self.encounter_checklist.SetCheckedItems([])
        self.on_selection_change(None)

