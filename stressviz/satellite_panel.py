# stressviz/satellite_panel.py
import wx
import wx.grid as wxgrid

# Import the SatStress GUI panel + helpers you showed earlier
# Adjust the module path to wherever SatelliteLayersPanel lives in your tree.
from SatStress.satstressgui import SatelliteLayersPanel, file_dialog, error_dialog

class StressVizSatellitePanel(SatelliteLayersPanel):
    """
    Drop-in satellite panel for StressViz:
    - Reuses SatStress' SatelliteLayersPanel (load/save/manual entry).
    - Adds a read-only matrix of current parameters.
    - Fires on_satellite_changed when the satellite is loaded/applied.
    """
    def __init__(self, *args, on_satellite_changed=None, **kw):
        self.on_satellite_changed_cb = on_satellite_changed
        super().__init__(*args, **kw)

        root = self.GetSizer()

        # Read-only matrix
        mat_box = wx.StaticBox(self, label=u"Parameter matrix")
        mat_sizer = wx.StaticBoxSizer(mat_box, wx.VERTICAL)
        self.matrix_grid = wxgrid.Grid(self)
        self.matrix_grid.CreateGrid(0, 2)
        self.matrix_grid.SetColLabelValue(0, "Parameter")
        self.matrix_grid.SetColLabelValue(1, "Value")
        self.matrix_grid.EnableEditing(False)
        mat_sizer.Add(self.matrix_grid, 1, wx.EXPAND | wx.ALL, 6)

        # Apply row (optional Love-number hook)
        apply_row = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_auto_love = wx.CheckBox(self, label=u"Compute Love numbers on apply")
        self.btn_apply = wx.Button(self, label=u"Apply to model")
        apply_row.Add(self.chk_auto_love, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL, 10)
        apply_row.AddStretchSpacer(1)
        apply_row.Add(self.btn_apply, 0)

        root.AddSpacer(10)
        root.Add(mat_sizer, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        root.Add(apply_row, 0, wx.EXPAND | wx.ALL, 6)

        # Events + initial populate
        self.btn_apply.Bind(wx.EVT_BUTTON, self._on_apply_clicked)

    # Let .sat.bak show in picker too
    def load(self, evt):
        try:
            file_dialog(
                self,
                message=u"Load from satellite file",
                style=wx.FD_OPEN,
                wildcard='Satellite files (*.satellite;*.sat;*.sat.bak)|*.satellite;*.sat;*.sat.bak|All|*.*',
                action=self.load_entries
            )
        except Exception as e:
            error_dialog(self, str(e), u"Satellite Error")

    def load_entries(self, filename):
        self.sc.load_satellite(filename)
        self.update_parameters()              # SatStress method (already there)
        self._notify_satellite_changed()

    def _on_apply_clicked(self, evt):
        # If SatStress bindings already sync self.parameters -> self.sc as you type,
        # we just need to refresh the matrix.

        if self.chk_auto_love.IsChecked():
            try:
                # Hook your Love-number compute here if desired
                # e.g., Diurnal(self.sc.satellite).calcLove()
                pass
            except Exception as e:
                error_dialog(self, f"Love number compute failed:\n{e}", u"Love number error")

        self._notify_satellite_changed()

    def _refresh_matrix_from_controls(self):
        rows = []
        # Global params (built from self.sc.satellite_vars)
        for (pname, _disp) in getattr(self.sc, "satellite_vars", []):
            ctrl = self.parameters.get(pname)
            if ctrl:
                rows.append((pname, ctrl.GetValue()))
        # Layer params (self.sc.layer_vars_d x self.sc.satlayers_d)
        for (layer_id, _v) in getattr(self.sc, "satlayers_d", []):
            for (p, _disp) in getattr(self.sc, "layer_vars_d", []):
                key = f"{p}_{layer_id}"
                ctrl = self.parameters.get(key)
                if ctrl:
                    rows.append((key, ctrl.GetValue()))

        g = self.matrix_grid
        if g.GetNumberRows() > 0:
            g.DeleteRows(0, g.GetNumberRows())
        g.AppendRows(len(rows))
        for r, (k, v) in enumerate(rows):
            g.SetCellValue(r, 0, str(k))
            g.SetCellValue(r, 1, "" if v is None else str(v))
        g.AutoSizeColumns()

    def _notify_satellite_changed(self):
        if self.on_satellite_changed_cb:
            sat_obj = getattr(self.sc, "satellite", None)
            self.on_satellite_changed_cb(sat_obj)
