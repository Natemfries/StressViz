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

# -------------- Modeless guide window --------------
class GettingStartedFrame(wx.Frame):
    def __init__(self, parent, docs_url: str | None = None):
        # Keep it light; allow resize so long content fits small screens
        super().__init__(parent,
                         title="Welcome to StressViz",
                         size=(680, 560))

        self._docs_url = docs_url

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # Header
        hdr = wx.StaticText(panel, label="Getting Started with StressViz")
        f = hdr.GetFont()
        f.SetPointSize(f.GetPointSize() + 4)
        try:
            f.MakeBold()
        except AttributeError:
            f.SetWeight(wx.FONTWEIGHT_BOLD)
        hdr.SetFont(f)
        vbox.Add(hdr, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Body (HtmlWindow)
        html = wxhtml.HtmlWindow(panel, style=wxhtml.HW_SCROLLBAR_AUTO)
        html.SetPage(_HTML_BODY)
        vbox.Add(html, 1, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)

        # Buttons
        btn_row = wx.BoxSizer(wx.HORIZONTAL)
        btn_row.AddStretchSpacer(1)

        if self._docs_url:
            btn_docs = wx.Button(panel, wx.ID_ANY, "Open Documentation")
            btn_row.Add(btn_docs, 0, wx.ALL, 6)
            self.Bind(wx.EVT_BUTTON, self._on_open_docs, btn_docs)

        btn_ok = wx.Button(panel, wx.ID_OK, "Continue")
        btn_ok.SetDefault()
        btn_row.Add(btn_ok, 0, wx.ALL, 6)
        vbox.Add(btn_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)

        panel.SetSizer(vbox)

        # Events
        self.Bind(wx.EVT_BUTTON, lambda _e: self.Close(), btn_ok)
        self.Bind(wx.EVT_CLOSE, self._on_close)

        # Placement
        try:
            self.CentreOnParent() if parent else self.Centre()
        except Exception:
            self.Centre()

    def _on_open_docs(self, _evt):
        url = self._docs_url
        if url:
            try:
                wx.LaunchDefaultBrowser(url)
            except Exception:
                pass

    def _on_close(self, evt):
        # Destroy safely; parent cleanup happens in show_getting_started
        self.Destroy()
        evt.Skip()

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
  <body style="font-family:-apple-system,Segoe UI,Arial; font-size: 12pt; line-height:1.35;">
    <p>Follow this quick path to your first plot:</p>
    <ol>
      <li><b>Load Satellite Parameters</b> — Select “Europa Preset” or load a <code>.sat</code> file.</li>
      <li><b>Compute Love Numbers</b> — or enter values manually.</li>
      <li><b>Select and load encounter</b> — or enter UTC manually and resolve true anomaly.</li>
      <li><b>Set Location</b> — Provide lat/lon if fields are empty.</li>
      <li><b>Calculate Stress</b> — Compute stress at the specified location.</li>
      <li><b>Enter Grid/Orbit Ranges</b> — if you plan to sweep ν or map a region.</li>
      <li><b>Open Stress Map</b> or <b>Scalar Plot</b> — to visualize results.</li>
    </ol>
    <p><b>Tips</b>:</p>
    <ul>
      <li>Use <i>Save Series</i> to export a sweep across ν values.</li>
      <li>Help → Getting Started (or press <b>F1</b>) to reopen this guide.</li>
    </ul>
  </body>
</html>
"""


