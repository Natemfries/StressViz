# stressviz/satstress_interface.py
"""
Bridges StressViz to the SatStress engine and GUI.

- Robustly load a SatStress Satellite from a .sat file.
- Launch the *actual* SatStress ScalarPlotPanel (with Basemap, sliders,
  vectors, save-series, lineaments/cycloids UI, etc.) inside a wx.Frame.
"""

from __future__ import annotations

import os, sys, types
from typing import Any, Tuple

import numpy as np
import wx
import importlib


__all__ = [
    "SatFileLoadError",
    "load_satellite_from_file",
    "SSScalarFrame",
    "open_satstress_scalar_plot",
]



_CONFIG_DEFAULTS = {
    "satellite_calculation": True,
    "grid_calculation": True,
    "Diurnal": True,
    "Nonsynchronous Rotation": False,
    "Obliquity": False,
    "Polar Wander": False,
    "SCALAR_PLOT_STEP": 30,
    "PLOT_LBOUND": -100,  # kPa
    "PLOT_UBOUND":  100,  # kPa
}

def _ensure_cfg_api(mod: types.ModuleType):
    """Give a config module the load/save API and required keys."""
    # keep a private store so load/save behave like SatStress expects
    store = {k: getattr(mod, k, v) for k, v in _CONFIG_DEFAULTS.items()}

    def load(key, default=None):
        if key in store:
            return store[key]
        if default is not None:
            return default
        raise KeyError(key)

    def save(**kwargs):
        store.update(kwargs)
        for k, v in kwargs.items():
            setattr(mod, k, v)

    mod.load = getattr(mod, "load", load)
    mod.save = getattr(mod, "save", save)
    for k, v in store.items():
        setattr(mod, k, v)

def _install_config_defaults():
    """Install a SatStress-like config before the GUI imports it."""
    for name in ("SatStress.config", "satstress.config", "config"):
        if name not in sys.modules:
            mod = types.ModuleType(name.rsplit(".", 1)[-1])
            sys.modules[name] = mod
        _ensure_cfg_api(sys.modules[name])

def _patch_already_imported_config_modules():
    """After the GUI is imported, patch the actual config module it used."""
    for name, mod in list(sys.modules.items()):
        if not isinstance(mod, types.ModuleType):
            continue
        if name.endswith(".config") and (name.startswith("SatStress") or name.startswith("satstress")):
            _ensure_cfg_api(mod)
        if name == "config":  # some forks import bare 'config'
            _ensure_cfg_api(mod)


def _resolve_stresscalc_class(ss):
    """Find StressCalc across common SatStress layouts."""
    if hasattr(ss, "StressCalc"):
        return ss.StressCalc
    for modpath in ("SatStress.satstress.stresscalc", "satstress.stresscalc", "SatStress.satstress"):
        try:
            mod = importlib.import_module(modpath)
            if hasattr(mod, "StressCalc"):
                return getattr(mod, "StressCalc")
        except Exception:
            pass
    # last resort: your project helper
    try:
        from .utils import get_satstress_handles
        _pkg, StressCalc, _Diurnal = get_satstress_handles()
        if StressCalc:
            return StressCalc
    except Exception:
        pass
    raise ImportError("Found SatStress but no usable StressCalc class in known locations.")

# ======================================================================================
# Satellite loading
# ======================================================================================

class SatFileLoadError(Exception):
    """Raised when a .sat file cannot be loaded into a SatStress Satellite."""


def _import_satellite_class():
    """
    Import SatStress' Satellite class.

    Expected path (typical tree):
      SatStress/satstress/satstress.py → class Satellite
    """
    from SatStress.satstress.satstress import Satellite  # type: ignore
    return Satellite


def load_satellite_from_file(path: str) -> Any:
    """
    Best-effort loader for SatStress Satellite objects from a .sat file path.
    Handles constructor signature differences across forks.

    Parameters
    ----------
    path : str
        Path to a .sat file.

    Returns
    -------
    Any
        A SatStress Satellite instance.

    Raises
    ------
    SatFileLoadError
        If no loading strategy succeeds.
    """
    if not os.path.isfile(path):
        raise SatFileLoadError(f"File not found: {path}")

    Satellite = _import_satellite_class()

    # 1) Most likely: constructor accepts a PATH STRING
    try:
        sat = Satellite(path)  # type: ignore[call-arg]
        return sat
    except TypeError:
        pass
    except Exception:
        pass

    # 2) Some versions expect a FILE OBJECT
    try:
        with open(path, "r") as f:
            sat = Satellite(f)  # type: ignore[call-arg]
            return sat
    except TypeError:
        pass
    except Exception:
        pass

    # 3) Rare: empty init + explicit reader method
    try:
        try:
            sat = Satellite(None)  # type: ignore[call-arg]
        except Exception:
            sat = Satellite()  # type: ignore[call-arg]
        for m in ("read_sat_file", "readSatFile", "from_file", "read", "load", "open", "open_file"):
            if hasattr(sat, m):
                obj = getattr(sat, m)(path)
                return obj if obj is not None else sat
    except Exception:
        pass

    raise SatFileLoadError(
        "Could not initialize Satellite with path/file, and no reader method worked."
    )


# ======================================================================================
# SatStress ScalarPlot wrapper (imports are lazy to avoid hard deps at module import)
# ======================================================================================

def _lazy_import_satstress():
    """Import satstress (engine) and the GUI ScalarPlotPanel on demand."""
    try:
        import satstress as ss  # type: ignore
    except Exception as e:
        raise ImportError(
            "Failed to import 'satstress'. Make sure SatStress is installed/activated."
        ) from e

    try:
        # Full GUI panel with Basemap and all plot behavior
        from SatStress.satstressgui import ScalarPlotPanel as SS_ScalarPlotPanel  # type: ignore
    except Exception as e:
        raise ImportError(
            "Failed to import 'SatStress.satstressgui.ScalarPlotPanel'. "
            "Ensure the SatStress GUI package (with Basemap) is on PYTHONPATH."
        ) from e

    return ss, SS_ScalarPlotPanel

def _prime_satstress_config():
    """
    Ensure the SatStress GUI's global config has the keys it expects.
    Some forks KeyError if these are missing.
    """
    try:
        cfg = importlib.import_module("SatStress.config")
    except Exception:
        return  # best effort

    # Keys + defaults the GUI expects
    defaults = {
        "satellite_calculation": True,   # <-- the one failing now
        "grid_calculation": True,
        "Diurnal": True,                 # enables orbit slider by default
        "Nonsynchronous Rotation": False,
        "Obliquity": False,
        "Polar Wander": False,
        "SCALAR_PLOT_STEP": 30,          # tick mark increment default
        "PLOT_LBOUND": -100,             # color scale defaults (kPa)
        "PLOT_UBOUND":  100,
    }

    # config.save/load API used by SatStress
    save = getattr(cfg, "save", None)
    load = getattr(cfg, "load", None)

    for k, v in defaults.items():
        try:
            # if not set, this usually raises
            _ = load(k) if callable(load) else getattr(cfg, k)
        except Exception:
            if callable(save):
                save(**{k: v})
            else:
                setattr(cfg, k, v)


class _GridShim:
    """Minimal grid object the SatStress GUI panel expects."""
    def __init__(self, lat_min, lat_max, lat_num, lon_min, lon_max, lon_num):
        self.lat_min = float(lat_min)
        self.lat_max = float(lat_max)
        self.lat_num = int(lat_num)
        self.lon_min = float(lon_min)
        self.lon_max = float(lon_max)
        self.lon_num = int(lon_num)

        # Harmless defaults for optional fields the panel may read
        self.grid_id = "stressviz"
        self.nsr_period_min = 1.0
        self.nsr_period_max = 1.0
        self.nsr_period_num = 1
        self.orbit_min = None
        self.orbit_max = None
        self.orbit_num = None
        self.time_min = 0.0
        self.time_max = 0.0
        self.time_num = 1


class _SCShim:
    """
    Bare-minimum SatStress 'controller' the ScalarPlotPanel talks to.

    Holds a StressCalc, a grid, a Satellite, and a parameter dict the GUI reads/writes.
    """
    def __init__(
        self,
        satellite: Any,
        stresscalc: Any,
        grid: _GridShim,
        orbit_min: float = 0.0,
        orbit_max: float = 360.0,
        orbit_num: int = 10,
    ):
        self._sat = satellite
        self.calc = stresscalc
        self.grid = grid

        self.parameters: dict[str, Any] = {
            # plotting knobs
            "projection": "cyl",
            "direction": "east",
            "field": "tens",
            "to_plot_principal_vectors": True,
            "to_plot_latitude_vectors": False,
            "to_plot_longitude_vectors": False,
            "to_plot_shear_vectors": False,
            "to_plot_lineaments": False,
            "to_plot_cycloids": False,
            "to_plot_triangles": True,
            "show_cycl_names": False,
            "to_plot_pw_markers": True,

            # >>> add these SatStress “Stresses tab” defaults <<<
            "satellite_calculation": True,          # SatStress checks this
            "grid_calculation": True, 
            "Diurnal": True,                         # enables orbit slider
            "Nonsynchronous Rotation": False,        # NSR off by default
            "Obliquity": False,
            "Polar Wander": False,
            "SCALAR_PLOT_STEP": 30,                  # tick mark increment default

            # grid bounds
            "LAT_MIN": grid.lat_min,
            "LAT_MAX": grid.lat_max,
            "LAT_NUM": grid.lat_num,
            "LON_MIN": grid.lon_min,
            "LON_MAX": grid.lon_max,
            "LON_NUM": grid.lon_num,

            # orbit / nsr series
            "ORBIT_MIN": float(orbit_min),
            "ORBIT_MAX": float(orbit_max),
            "ORBIT_NUM": int(orbit_num),
            "TIME_MIN": 0.0,
            "nsr_time": 0.0,
            "TIME_NUM": 0,

            # metadata
            "SYSTEM_ID": getattr(satellite, "system_id", "system"),
        }


        # flags the panel checks
        self._changed = True
        self.calc_changed = True
        self.projection_changed = False

        # placeholders referenced by cycloid/lineament UI
        self.cycloids = {}
        self.params_for_cycloids = {}
        self.many_changed = False
        self.cycloid_changed = False
        self.cyc = None

    # --- API the panel calls ---
    def get_parameter(self, typ, key, default=None):
        v = self.parameters.get(key, default)
        try:
            return typ(v)
        except Exception:
            return default

    def set_parameter(self, key, val):
        self.parameters[key] = val
        self._changed = True

    def get_satellite(self):
        return self._sat

    def changed(self):
        was = self._changed
        self._changed = False
        return was

    # file ops stubs the panel calls (you can wire later)
    def load_netcdf(self, *_):  # noqa: D401
        """No-op: override to wire NetCDF loading."""
        pass

    def save_netcdf(self, *_):  # noqa: D401
        """No-op: override to wire NetCDF saving."""
        pass


def _make_wrapped_panel_class(SS_ScalarPlotPanel):
    """
    Build a subclass that injects sc/grid before the SatStress base __init__ runs.
    Defined dynamically so we can lazy-import the base panel.
    """
    class _WrappedSSScalarPanel(SS_ScalarPlotPanel):  # type: ignore[misc]
        def __init__(self, parent: wx.Window, scshim: _SCShim):
            self.sc = scshim
            self.grid = scshim.grid
            super().__init__(parent)

    return _WrappedSSScalarPanel


class SSScalarFrame(wx.Frame):
    """
    Standalone frame that hosts the real SatStress ScalarPlotPanel.

    Create and Show() this to get the full SatStress plot behavior in StressViz.
    """
    def __init__(
        self,
        parent: wx.Window,
        satellite: Any,
        diurnal_cls: type,
        lat_min: float = -90,
        lat_max: float = 90,
        lat_num: int = 91,
        lon_min: float = -180,
        lon_max: float = 180,
        lon_num: int = 181,
        orbit_min: float = 0,
        orbit_max: float = 360,
        orbit_num: int = 12,
        nu0_deg: float | None = None,
        title: str = "SatStress Plot",
    ):
        # 1) window
        super().__init__(parent, title=title, size=(1100, 720))

        # 2) make sure the GUI sees required config keys BEFORE importing panels
        _install_config_defaults()

        # 3) import SatStress GUI + wrap its ScalarPlotPanel
        ss, SS_ScalarPlotPanel = _lazy_import_satstress()
        WrappedPanel = _make_wrapped_panel_class(SS_ScalarPlotPanel)

        _patch_already_imported_config_modules()

        # 4) build StressCalc robustly across forks
        StressCalc = _resolve_stresscalc_class(ss)
        try:
            calc = StressCalc([diurnal_cls], satellite=satellite)
        except TypeError:
            # older signature
            calc = StressCalc([diurnal_cls])

        # 5) controller shim + GUI panel
        grid = _GridShim(lat_min, lat_max, lat_num, lon_min, lon_max, lon_num)
        sc = _SCShim(satellite, calc, grid, orbit_min, orbit_max, orbit_num)

        panel = WrappedPanel(self, sc)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(sizer)

        # 6) initial draw
        panel.update_parameters()
        panel.plot()

        # 7) set initial ν on the orbit slider if provided
        if nu0_deg is not None:
            try:
                panel.orbit_pos = float(nu0_deg) % 360.0  # SatStress stores degrees
                if hasattr(panel, "scp") and hasattr(panel.scp, "orbit_slider"):
                    panel.scp.orbit_slider.set_val(panel.orbit_pos)
                panel.plot()
            except Exception:
                pass

        self.Centre()
        self.Show()



def open_satstress_scalar_plot(
    parent: wx.Window,
    satellite: Any,
    diurnal_cls: type,
    lats_deg: np.ndarray | list[float],
    lons_deg: np.ndarray | list[float],
    orbit_series: Tuple[float, float, int] = (0.0, 360.0, 10),
    nu0_deg: float | None = None,
    title: str = "SatStress Plot (StressViz)",
) -> SSScalarFrame:
    """
    Convenience helper: derive ranges from arrays and show the real SatStress plot.

    Parameters
    ----------
    parent : wx.Window
        Parent window for the frame.
    satellite : Any
        SatStress Satellite instance.
    diurnal_cls : type
        The SatStress Diurnal class (not an instance).
    lats_deg, lons_deg : array-like
        Latitude and longitude axes in degrees.
    orbit_series : (float, float, int)
        (ORBIT_MIN, ORBIT_MAX, ORBIT_NUM) for the slider/series buttons.
    nu0_deg : float, optional
        Initial orbital position to set the slider at.
    title : str
        Frame title.

    Returns
    -------
    SSScalarFrame
        The created (and already shown) frame.
    """
    lats = np.asarray(lats_deg, float)
    lons = np.asarray(lons_deg, float)

    lat_min, lat_max = float(lats.min()), float(lats.max())
    lon_min, lon_max = float(lons.min()), float(lons.max())
    lat_num = int(lats.size)
    lon_num = int(lons.size)

    omin, omax, onum = orbit_series

    frame = SSScalarFrame(
        parent=parent,
        satellite=satellite,
        diurnal_cls=diurnal_cls,
        lat_min=lat_min,
        lat_max=lat_max,
        lat_num=lat_num,
        lon_min=lon_min,
        lon_max=lon_max,
        lon_num=lon_num,
        orbit_min=float(omin),
        orbit_max=float(omax),
        orbit_num=int(onum),
        nu0_deg=nu0_deg,
        title=title,
    )
    return frame











