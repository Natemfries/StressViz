#About
StressViz is a Python 3–compatible extension of the SatStressGUI program, which models tidal stresses in icy satellites. It integrates orbital position plotting with SatStress tidal stress calculations, providing a tool to analyze stress environments in support of Europa Clipper plume detections. StressViz is being developed for Europa, but could also be applied to other icy satellites that are compatible with SatStress. StressViz is being developed at the Jet Propulsion Laboratory, California Institute of Technology.


#Installation
Requirements
- Python 3.10+
- wxPython ≥ 4.2
- numpy, scipy, matplotlib
- netCDF4 (optional, but supported)
- Basemap (for map projections, legacy but supported)

#Setup
Clone the repository and create the conda environment:
git clone https://github.com/Natemfries/StressViz.git
cd StressViz
conda env create -f environment.yml
conda activate stressviz

Then run:
python main.py

#Quick Run
- "Europa Preset"
- "Compute Love Numbers"
- "Select Encounter"
- "Calculate Stress"
- "Plot" Europa 

#Known Issues
- Conda may attempt to use the default Anaconda channel and hit rate limits. Ensure `conda-forge` is used with strict priority.
- wxPython installation may fail outside of Conda environments.
- "Load from file" button to input satellite parameters doesn't work. Bypass this with the "Europa Preset" button.
- Colorbar range for stress plots doesn't update when bounds are changed.

#Encounters
Includes pre-loaded planned Europa Clipper and JUICE encounters, as well as past plume search observations. Metadata for these encounters is in data/plume_observations.txt, data/JUICE_Europa_flybys.txt, and data/Europa_Encounters_21F31_V7_LP01_ver2.txt. Positive (Y) and negative (N) detections are based on the references cited in Europa_Plume_Observations.xslx.


#Stress Plots
The stress plots show the evolution of the diurnal stress map throughout one orbit. They are generated using SatStressGUIV6.0. 


#Orbital Position Plots
Orbital plots show the mean anomaly of Europa at the time of a given encounter, which can then be matched to the corresponding stress plot. 


#Credits
StressViz builds on the foundation of SatStressGUI, a program developed at the Jet Propulsion Laboratory, California Institute of Technology, to model tidal stresses in icy satellites.
- SatStressGUI V5.0 was created through the efforts of: Zane Selvans, Mikael Beuthe, Jonathan Kay, Lee Tang, Simon Kattenhorn, C.M. Cooper, Emily S. Martin, David Dubois, Ben J. Ayton, Jessica B. Li, Andre Ismailyan, Peter Sinclair, Nhu Doan, Chad Harper, and Alex Patthoff.
- SatStressGUI V6.0 (Python 3 modernization and updates) is being developed and maintained by Nathan Friesenhahn at NASA JPL.
- StressViz is an extension of SatStressGUI V6.0, adding orbital position visualization, encounter metadata integration, and new user interface features.

