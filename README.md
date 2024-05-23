# Global Patterns in structural habitat complexity.
Briefly, the scripts contained in the respository follow the analysis sequence presented in 'Cant J., Schiettekatte N., Madin E.M.P., Madin J. & Dornelas M. (In prep). Unifying assessments of global habitat complexity', which entailed: generating global maps of selected complexity indices; assessing the scale ensitivity of these complexity indices; quantifying pairwise geometric relationships; evaluating ecosystem and land use type associated complexity characteristics; and visualising these quantified relationships and patterns. The global maps and data files created and used as part of these analyses have been archived seperatly and can be accessed at (**Insert Link**).

For further information please do not hesitate to contact **James Cant** at *james.cant91@gmail.com*.
---

## Script Details:

***1. Build Global DEM.R***
This script takes a series of raster tiles downloaded from [GEBCO](https://www.gebco.net/) detailing spatial topographical and bathymetric patterns, to form a single global elevation map (referred to as a Global Digital Elevation Model or DEM).

***2. Resolution Sensitivity.R***
This script tests the sensitivity of various complexity measures to changes in the scale at which they are computed by repeatedly calculating each measure for a series of successively smaller pixel cells selected from the same region of our global DEM. 

***3. Estimating Global Complexity.R***
Using the DEM generated above this script calculates global patterns in the geometic structural measures of *height range* (*ΔH*), *rugosity* (*R*), and *fractal dimension* (*D*), at a resolution of 1.87km. Please see the manuscript associated with this Github repository and also Torres-Pilliza et al. (2020) *Nat. Ecol. Evol.* **4** 1495-1501 for further details as to the definitions of these selected complexity measures. 

***4. Analysising Complexity patterns.R***
Using the maps generated showcasing global-scale patterns in the complexity measures of *height range*, *rugosity*, and *fractal dimension*, this script evaluates the pairwise relationships between these complexity variables, and how surface complexities differ between marine and terrestrial environments.

***5. Visualising Complexity.R***
Finally, this script can be used to generate plots showing (1) the global distribution of *height range*, *rugosity*, and *fractal dimension* estimates, (2) the relative *rugosity* and *fractal dimension* characteristics of selected global features, and (3) how estimates of *rugosity* and *fractal dimension* differ across different ecosystem types and land use classifications. This script also contains a small analysis quantifying the relationship between *rugosity* and *fractal dimension* across different ecosystem types.

### Enjoy!!
