# stressviz/onboarding.py
# Lightweight, modeless "Getting Started" window with first-run logic.
# - Non-blocking (modeless) so it won't slow startup
# - Defer showing until after the main frame paints
# - Keeps a strong reference to avoid GC closing the window
# - Optional docs button
# - Safe font bolding across wx versions

import wx
import wx.html as wxhtml

# Auto-show only once per user. Set to False to never auto-open.
AUTO_SHOW_FIRST_RUN = True

# -------------- First-run persistence --------------
def _cfg():
    return wx.Config(appName="StressViz", vendorName="StressViz")

def _is_first_run() -> bool:
    if not AUTO_SHOW_FIRST_RUN:
        return False
    return not _cfg().ReadBool("first_run_seen", False)

def _mark_first_run_seen():
    if not AUTO_SHOW_FIRST_RUN:
        return
    c = _cfg()
    c.WriteBool("first_run_seen", True)
    c.Flush()

class GettingStartedPanel(wx.Panel):
    def __init__(self, parent, docs_url: str | None = None):
        super().__init__(parent)

        self._docs_url = docs_url

        vbox = wx.BoxSizer(wx.VERTICAL)

        # Header
        hdr = wx.StaticText(self, label="Getting Started with StressViz")
        f = hdr.GetFont()
        f.SetPointSize(f.GetPointSize() + 4)

        try:
            f.MakeBold()
        except AttributeError:
            f.SetWeight(wx.FONTWEIGHT_BOLD)

        hdr.SetFont(f)
        vbox.Add(hdr, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Body
        self._html_zoom = 1.0

        self.html = wxhtml.HtmlWindow(self, style=wxhtml.HW_SCROLLBAR_AUTO)
        self._render_html()
        vbox.Add(self.html, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        # Buttons
        zoom_row = wx.BoxSizer(wx.HORIZONTAL)
        zoom_row.AddStretchSpacer(1)

        btn_zoom_out = wx.Button(self, wx.ID_ANY, "Zoom Out")
        btn_zoom_reset = wx.Button(self, wx.ID_ANY, "Reset")
        btn_zoom_in = wx.Button(self, wx.ID_ANY, "Zoom In")

        zoom_row.Add(btn_zoom_out, 0, wx.ALL, 4)
        zoom_row.Add(btn_zoom_reset, 0, wx.ALL, 4)
        zoom_row.Add(btn_zoom_in, 0, wx.ALL, 4)

        self.Bind(wx.EVT_BUTTON, self._on_zoom_out, btn_zoom_out)
        self.Bind(wx.EVT_BUTTON, self._on_zoom_reset, btn_zoom_reset)
        self.Bind(wx.EVT_BUTTON, self._on_zoom_in, btn_zoom_in)

        vbox.Add(zoom_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)
        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        btn_row.AddStretchSpacer(1)

        if self._docs_url:
            btn_docs = wx.Button(self, wx.ID_ANY, "Open Documentation")
            btn_row.Add(btn_docs, 0, wx.ALL, 6)
            self.Bind(wx.EVT_BUTTON, self._on_open_docs, btn_docs)

        vbox.Add(btn_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)

        self.SetSizer(vbox)

    def _on_open_docs(self, _evt):
        url = self._docs_url
        if url:
            try:
                wx.LaunchDefaultBrowser(url)
            except Exception:
                pass

    def _render_html(self):
        base_size = int(18 * self._html_zoom)

        self.html.SetFonts(
            normal_face="",
            fixed_face="",
            sizes=[
                max(6, base_size - 4),
                max(7, base_size - 2),
                max(8, base_size),
                max(10, base_size + 2),
                max(12, base_size + 4),
                max(14, base_size + 6),
                max(16, base_size + 8),
            ],
        )

        self.html.SetPage(_HTML_BODY)


    def _on_zoom_in(self, _evt):
        self._html_zoom = min(self._html_zoom + 0.1, 2.0)
        self._render_html()


    def _on_zoom_out(self, _evt):
        self._html_zoom = max(self._html_zoom - 0.1, 0.7)
        self._render_html()


    def _on_zoom_reset(self, _evt):
        self._html_zoom = 1.0
        self._render_html()

# -------------- Modeless guide window --------------
class GettingStartedFrame(wx.Frame):
    def __init__(self, parent, docs_url: str | None = None):
        super().__init__(
            parent,
            title="Welcome to StressViz",
            size=(680, 560)
        )

        panel = GettingStartedPanel(self, docs_url=docs_url)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(sizer)

        self.Bind(wx.EVT_CLOSE, self._on_close)

        try:
            self.CentreOnParent() if parent else self.Centre()
        except Exception:
            self.Centre()

# -------------- Public API --------------
# Keep a global fallback reference if no parent is provided
_global_keepalive = set()

def show_getting_started(parent: wx.Window | None, docs_url: str | None = None):
    """
    Open the modeless Getting Started window (for Help menu / F1).
    Ensures a strong reference so it doesn't get GC'd.
    """
    # If parent exists and already has one open, focus it
    if parent is not None and getattr(parent, "_getting_started_win", None):
        try:
            win = parent._getting_started_win
            win.Raise()
            win.RequestUserAttention() if hasattr(win, "RequestUserAttention") else None
            return
        except Exception:
            parent._getting_started_win = None  # fall through to recreate

    frame = GettingStartedFrame(parent, docs_url=docs_url)

    # --- keep a strong reference so GC doesn't kill it ---
    if parent is not None:
        parent._getting_started_win = frame

        def _cleanup(_evt):
            # clear ref when the window is destroyed
            if getattr(parent, "_getting_started_win", None) is frame:
                parent._getting_started_win = None
            _evt.Skip()

        frame.Bind(wx.EVT_WINDOW_DESTROY, _cleanup)
    else:
        # No parent? Keep a module-level reference
        _global_keepalive.add(frame)

        def _cleanup_global(_evt):
            _global_keepalive.discard(frame)
            _evt.Skip()

        frame.Bind(wx.EVT_WINDOW_DESTROY, _cleanup_global)

    frame.Show()
    try:
        if parent:
            parent.Raise()
    except Exception:
        pass


def maybe_show_getting_started(parent: wx.Window | None, docs_url: str | None = None):
    """
    Defer opening until after the main frame paints.
    Shows only on first run if AUTO_SHOW_FIRST_RUN is True.
    """
    def _open():
        if _is_first_run():
            _mark_first_run_seen()
            show_getting_started(parent, docs_url=docs_url)

    # Defer so startup feels instant
    wx.CallAfter(_open)

# -------------- Static HTML content --------------
_HTML_BODY = """
<html>
  <body style="font-family:-apple-system,Segoe UI,Arial; line-height:1.35;">
    <p>StressViz initializes with the default Europa satellite parameters from /data/EuropaSample.sat.</p>
    <p>To continue with Europa defaults:</p>
    <ol>
      <li><b>Go to the Plots tab.</b></li>
      <li><b>Select Encounters</b> - Select desired subset of encounters from the dropdown menu. These will be displayed in the "Observation Events" panel.</li>
      <li><b>Plot on Orbit</b> - Plot selected encounters on the Orbit Plot</li>
      <li><b>Move Stress Plot</b> - Move the Stress Plot to match with the mean anomaly (M) of the top selected encounter. The Stress Plot marker on the Orbit Plot will also move to reflect this.</li>
      <li><b>Show Nearby Events</b> - Searches for encounters within +/-10&deg M of the top selected encounter. This plots all on the Orbit Plot and adds them to he Observation Events panel.</li>
    </ol>
    <p>To manually input an encounter not in the internal database:</p>
    <ol>
      <li><b>Go to the Controls tab.</b></li>
      <li><b>Select Encounter or Input Location/Time</b> - Option to load Clipper encounter from dropdown or manual input.</li>
      <li><b>Enter Date-Time</b> - Enter the datetime in UTC (ISO) format (yyyy-mm-ddThh:mm:ssZ).</li>
      <li><b>Resolve Orbital Phase<b> - Queries true anomaly of Europa at given date-time from JPL New Horizons and converts to mean anomaly (M) for use with SatStress code.</li>
      <li><b>Enter Point Location</b> - Enter Lat/Lon for the desired point on Europa's surface. This is used to compute stress at that specific point. Arbitrary values may be used if this is not a priority.
      <li><b>Compute Stress</b> - Compute the stress at the point of the inputted lat/lon and orbital position.</li>
      <li><b>Plot</b> - Select <i>Plot</i> to push the manually inputted encounter to the Observation Events panel.</li>
    </ol>
    <p>To use satellite parameters other than the default Europa:</p>
    <ol>
      <li><b>Go to the Controls tab.</b></li>
      <li><b>Load Satellite Parameters</b> — Select <i>Europa Preset</i> or load a <code>.sat</code> file. ***Load a file is currently not functional, but you may manually input desired values</li>
      <li><b>Compute Love Numbers</b> — or enter values manually.</li>
      <li><b>Select and load encounter</b> — or enter UTC manually and resolve true anomaly.</li>
      <li><b>Set Location</b> — Provide lat/lon if fields are empty.</li>
      <li><b>Compute Stress</b> — Compute stress at the specified location.</li>
      <li><b>Enter Grid/Orbit Ranges</b> — if you plan to sweep ν or map a region.</li>
      <li><b>Open Stress Map</b> or <b>Scalar Plot</b> — to visualize results.</li>
    </ol>
    <p><b>Tips</b>:</p>
    <ul>
      <li>Use <i>Save Stress Series</i> to export a series of stress plots across the given range of ν values.</li>
      <li>Use <i>Save Orbit Plot</i> to export the current Orbit plot.</li>
      <li>Help → Getting Started (or press <b>F1</b>) to reopen this guide.</li>
    </ul>
  </body>
</html>
"""


