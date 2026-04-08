# stressviz/stress_panel.py

import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

class StressPlotPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        self.figure = Figure(figsize=(5, 5))
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.ax = self.figure.add_subplot(111)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND | wx.ALL, 5)
        self.SetSizer(sizer)

    def plot_stress(self, stress_data, lat, lon):
        """
        Plot the stress field given lat/lon and stress values.
        """
        self.ax.clear()
        self.ax.set_title("Stress Field")
        self.ax.set_xlabel("Longitude")
        self.ax.set_ylabel("Latitude")
        self.ax.pcolormesh(lon, lat, stress_data, shading='auto')
        self.canvas.draw()

