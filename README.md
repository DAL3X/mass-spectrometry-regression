# Basics
A Matlab script to analyse data captured by a [OmniStar GSD 350 O](https://www.pfeiffer-vacuum.com/de/produkte/messung-analyse/analysegeraete/gasanalyse/gasanalyse-im-druckbereich-bis-1-000-hpa/omnistar-gsd-350-o/) for gas analysis.
It was developed as a part of scientific study about the thermal degradation of solid electrolyte interfaces (SEI) in lithium-ion batteries.

# Usage and Outputs
The script was used to track the temperature and following substances:

![substances](https://github.com/DAL3X/mass-spectrometry-regression/blob/master/pictures/interest.png)

After reading in the data the script will auto detect heating intervals and apply regression to each interval respectively for every substance.
The results are saved together with the mean squared error and [R squared value](https://en.wikipedia.org/wiki/Coefficient_of_determination).
The regression results will also be visualized for instant feedback:

![regression](https://github.com/DAL3X/mass-spectrometry-regression/blob/master/pictures/example.png)

The user can also force the usage of specific base functions for arbitrary intervals per substance.
