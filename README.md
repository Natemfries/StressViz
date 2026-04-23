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
Includes pre-loaded planned Europa Clipper and JUICE encounters, as well as past plume search observations. Metadata for these encounters is in data/plume_observations.txt, data/JUICE_Europa_flybys.txt, and data/Europa_Encounters_21F31_V7_LP01_ver2.txt. 

Encounters are noted as positive (Y), negative (N), or contested (*) plume findings based on the following references. Contested findings were initially reported as positive detections, and then were later reported as negative detections in a separate paper. Full references for each observation can be found in Europa_Plume_Observations.xslx.
Positive/Contested Detections:
- HST 2012-12-30T18:49:00Z - Initially reported as a positive detection according to Roth et al. (2014), but has been updated to a negative detection in Roth et al., (2026).
- HST 2014-01-26T18:05:00Z HST 2014-03-17T11:47:00Z, HST 2014-04-04T05:20:00Z, and HST 2016-02-22T00:00:00Z - All reported as positive detections in Sparks et al. (2016,2017), but were subsequently refuted by Giono et al., (2020).
- Keck 2016-04-26T05:32:00Z - Reported as a positive detection in Paganini et al., (2019).
- Galileo 1997-12-16T12:00:59Z - Reported as a positive detection in Jia et al., (2018).
- Subaru 2021-07-17T10:21:00Z - Reported as a positive detection in Kimura et al., (2024).

#Stress Plots
The stress plots show the evolution of the diurnal stress map throughout one orbit. They are generated using SatStressGUIV6.0. 


#Orbital Position Plots
Orbital plots show the mean anomaly of Europa at the time of a given encounter, which can then be matched to the corresponding stress plot. 


#Credits
StressViz builds on the foundation of SatStressGUI, a program developed at the Jet Propulsion Laboratory, California Institute of Technology, to model tidal stresses in icy satellites.
- SatStressGUI V5.0 was created through the efforts of: Zane Selvans, Mikael Beuthe, Jonathan Kay, Lee Tang, Simon Kattenhorn, C.M. Cooper, Emily S. Martin, David Dubois, Ben J. Ayton, Jessica B. Li, Andre Ismailyan, Peter Sinclair, Nhu Doan, Chad Harper, and Alex Patthoff.
- SatStressGUI V6.0 (Python 3 modernization and updates) is being developed and maintained by Nathan Friesenhahn at NASA JPL.
- StressViz is an extension of SatStressGUI V6.0, adding orbital position visualization, encounter metadata integration, and new user interface features.

#References
- G. L. Villanueva et al. ,Endogenous CO2 ice mixture on the surface of Europa and no detection of    plume activity.Science381,1305-1308(2023).DOI:10.1126/science.adg4270
- Jia, Xianzhe, et al. "Evidence of a plume on Europa from Galileo magnetic and plasma wave signatures." Nature Astronomy 2.6 (2018): 459-464.
- Kimura, Jun, et al. "A search for water vapor plumes on Europa by spatially resolved spectroscopic observation using Subaru/IRCS." Publications of the Astronomical Society of Japan 76.6 (2024): 1302-1308.
- Paganini, Lucas, et al. "A measurement of water vapour amid a largely quiescent environment on Europa." Nature Astronomy 4.3 (2020): 266-272.
- Lorenz Roth et al. ,Transient Water Vapor at Europa’s South Pole.Science343,171-174(2014).DOI:10.1126/science.1247051
- Roth, L., et al. "Europa’s Lyman-α emissions from HST/STIS observations." (2026).
- Sparks, William B., et al. "Probing for evidence of plumes on Europa with HST/STIS." The Astrophysical Journal 829.2 (2016): 121.
- Sparks, William B., et al. "Active cryovolcanism on Europa?." The Astrophysical Journal Letters 839.2 (2017): L18.
- Sparks, W. B., et al. "A search for water vapor plumes on Europa using SOFIA." The Astrophysical Journal Letters 871.1 (2019): L5.
- Villanueva, Geronimo L., et al. "Endogenous CO2 ice mixture on the surface of Europa and no detection of plume activity." Science 381.6664 (2023): 1305-1308.