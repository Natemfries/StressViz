#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import wx

# ----------------- Paths -----------------
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
SATSTRESS_DIR = os.path.join(BASE_DIR, "SatStress")
if os.path.isdir(SATSTRESS_DIR) and SATSTRESS_DIR not in sys.path:
    sys.path.insert(0, SATSTRESS_DIR)

# ----------------- StressViz imports -----------------
from stressviz.control_panel import AnalysisControlPanel

# Onboarding (modeless, deferred)
from stressviz.onboarding import maybe_show_getting_started, show_getting_started

# Optional docs URL for the onboarding window (or leave as None)
DOCS_URL = None  # e.g., "https://github.com/YourOrg/StressViz#readme"

# Optional: encounters TXT loader (used by deferred autoload)
try:
    from stressviz.encounters_io import load_europa_min  # returns a DataFrame
except Exception:
    load_europa_min = None


class StressVizFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="StressViz", size=(1150, 750))

        # Central control panel
        self.control = AnalysisControlPanel(self)

        root = wx.BoxSizer(wx.VERTICAL)
        root.Add(self.control, 1, wx.EXPAND | wx.ALL, 8)
        self.SetSizer(root)
        self.Centre()

        # Menus, accelerators, status bar
        self._build_menu()
        self.CreateStatusBar()
        self.SetStatusText("Ready")

        # Defer non-essential work so the window appears instantly
        wx.CallAfter(self._post_show_deferred)

    # ----------------- UI wiring -----------------
    def _build_menu(self):
        menubar = wx.MenuBar()

        # Help menu
        help_menu = wx.Menu()
        self.ID_GETTING_STARTED = wx.NewIdRef()
        item_gs = help_menu.Append(
            self.ID_GETTING_STARTED,
            "Getting Started\tF1",
            "Open the Getting Started guide",
        )
        menubar.Append(help_menu, "&Help")
        self.SetMenuBar(menubar)

        # Handlers
        self.Bind(wx.EVT_MENU, self.on_open_getting_started, item_gs)

        # Accelerator for F1 (in case platform menus differ)
        accel_tbl = wx.AcceleratorTable([(wx.ACCEL_NORMAL, wx.WXK_F1, self.ID_GETTING_STARTED)])
        self.SetAcceleratorTable(accel_tbl)

    def on_open_getting_started(self, _evt):
        show_getting_started(self, docs_url=DOCS_URL)

    # ----------------- Deferred work after first paint -----------------
    def _post_show_deferred(self):
        """
        Runs after the frame is shown. Keep it light and non-blocking.
        """
        # 1) Auto-open Getting Started (modeless, first-run gated inside onboarding)
        maybe_show_getting_started(self, docs_url=DOCS_URL)

        # 2) Autoload encounters from TXT (optional)
        wx.CallAfter(self._autoload_encounters_txt)

    # ----------------- Data autoload (deferred) -----------------
    def _autoload_encounters_txt(self):
        """
        Optional: auto-load encounters into the Point panel dropdown from a TXT.
        Requires: stressviz.encounters_io.load_europa_min and
                  AnalysisControlPanel.load_encounters_df.
        """
        try:
            if load_europa_min is None:
                return
            if not hasattr(self.control, "load_encounters_df"):
                return

            txtfile = os.path.join(BASE_DIR, "data", "Europa_Encounters_21F31_V7_LP01_ver2.txt")
            if os.path.isfile(txtfile):
                df = load_europa_min(txtfile)
                # If this is heavy for some users, thread it in the future and CallAfter the UI update.
                self.control.load_encounters_df(df)
                try:
                    self.control.add_plume_observations()
                except Exception as e:
                    import traceback
                    print(f"[StressViz] add_plume_observations failed: {e}")
                    traceback.print_exc()
                    
                self.SetStatusText("Encounters loaded")
        except Exception as e:
            import traceback
            print(f"[StressViz] Autoload encounters failed: {e}")
            traceback.print_exc()
            self.SetStatusText("Encounter autoload failed (see console)")


class App(wx.App):
    def OnInit(self):
        # macOS nicety: set app name for menus/prefs (harmless elsewhere)
        try:
            self.SetAppName("StressViz")
        except Exception:
            pass

        frame = StressVizFrame()
        frame.Show()

        # Nothing blocking here—onboarding/autoload are deferred in the frame
        return True


if __name__ == "__main__":
    App(False).MainLoop()
