# Predicting-Detectability-of-Transiting-Exoplanets

This repository contains a Python script that models the signal-to-noise ratio and minimum detectable exoplanet transit depth for a small-aperture telescope. 
This is a user-friendly astronomy tool geared towards amateurs interested in studying transiting exoplanets.

To use this code, edit all of the data in the "INPUTS" section. The rest of the code is ready to run. This code uses the commonly cited "CCD Equation" to calculate signal-to-noise ratio, which is converted to minimum detectable transit depth in parts-per-thousand. Its accuracy has been tested using both amateur and research grade equipment.

The model accounts for:
- Telescope aperture and efficiency
- Sky background
- Read noise and dark current
- PSF concentration and saturation limits
- Exposure time and stellar magnitude

Requirements:
- Python 3.x
- numpy
- matplotlib

The model will output a graph showing minimum detectable transit depth on a transit depth vs. magnitude graph. Any targets that fall above the blue line and outside of the red bar should be detectable with your equipment. 
