# stressviz/ui_style.py
import wx

DEFAULT_PT = 15

def instr_label(parent, text, pt=DEFAULT_PT):
    lbl = wx.StaticText(parent, label=text)
    f = lbl.GetFont()
    f.SetPointSize(pt)
    f.MakeBold()
    lbl.SetFont(f)
    return lbl

def style_staticbox(box_or_sizer, pt=DEFAULT_PT, bold=True):
    """
    Accepts either a wx.StaticBox or a wx.StaticBoxSizer and styles the label.
    """
    if isinstance(box_or_sizer, wx.StaticBoxSizer):
        box = box_or_sizer.GetStaticBox()
    elif isinstance(box_or_sizer, wx.StaticBox):
        box = box_or_sizer
    else:
        raise TypeError(
            f"style_staticbox expected wx.StaticBox or wx.StaticBoxSizer, got {type(box_or_sizer)}"
        )

    f = box.GetFont()
    f.SetPointSize(pt)
    if bold:
        f.MakeBold()
    box.SetFont(f)


def apply_font(widget: wx.Window, pt=DEFAULT_PT, bold=False):
    """Set font on any widget (StaticText, TextCtrl, etc.)."""
    f = widget.GetFont()
    f.SetPointSize(pt)
    if bold:
        f.MakeBold()
    widget.SetFont(f)
