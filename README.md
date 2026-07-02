# Global Patterns in structural habitat complexity.
Briefly, the scripts contained in the repository follow the analysis sequence presented in 'Cant J., Schiettekatte N., Seijmonsbergen A.C., Rijsdijk K.F, Madin E.M.P., Madin J. & Dornelas M. (in review). Landscape structural complexity influences the global distribution of ecosystems and people', which entailed: generating global maps of selected complexity indices; assessing the scale sensitivity of these complexity indices; quantifying pairwise geometric relationships; evaluating ecosystem and land use type associated complexity characteristics; visualising these quantified relationships and patterns; and exploring their associated biotic and abiotic impacts.

For further information, please do not hesitate to contact **James Cant** at *james.cant91@gmail.com*.

## Script Details:

***1. 00_RunScript.R:*** This script centralises the package installation and file directory pathway designation for the whole analysis pipeline.

***2. 01_Build Global DEM.R:*** This script takes a series of raster tiles downloaded from [GEBCO](https://www.gebco.net/) that detail spatial topographic and bathymetric patterns and generates a single global elevation map (referred to as a Global Digital Elevation Model, or DEM).

***3. 02_Estimating Global Complexity.R:*** Using the DEM generated above, this script calculates global patterns in the geometric structural measures of *height range* (*ΔH*), *rugosity* (*R*), and *fractal dimension* (*D*), at a resolution of 1.87km. Please see the manuscript associated with this GitHub repository and also Torres-Pilliza et al. (2020) *Nat. Ecol. Evol.* **4** 1495-1501 for further details as to the definitions of these selected complexity measures. 

***4. 03_Analysing Complexity patterns.R:*** Using the maps generated showcasing global-scale patterns in the complexity measures of *height range*, *rugosity*, and *fractal dimension*, this script evaluates the pairwise relationships between these complexity variables, and how surface complexities differ between marine and terrestrial environments.

***5. 04_Ecosystem type and Land use comparisons.R:*** This script can be used to generate plots showing (1) the global distribution of *height range*, *rugosity*, and *fractal dimension* estimates, (2) the relative *rugosity* and *fractal dimension* characteristics of selected global features, and (3) how estimates of *rugosity* and *fractal dimension* differ across different ecosystem types and land use classifications.

***6. 05_Estimating Global Climate Variability.R:*** The analyses presented within the script do not currently feature in the manuscript itself, but they comprise extracting climatology data from WorldClim for quantifying decadal patterns in mean monthly standard deviation in temperature across the globe. The output derived from this script is a map showing spatial patterns in the standard deviation of monthly temperature variability across the globe.

***7. 06_Evaluating Complexity associations.R:*** This script centres around analyses exploring the association between each of the measures of geometric complexity (*height range, rugosity, and fractal dimension*) and (1) *geodiversity* and (2) *human population density*. 

***8. 07_Island Complexity Analysis.R:*** This script centres around analyses evaluating the extent to which including the geometric complexity measures of height range, rugosity, and fractal dimension adds value to assessments of island-species area relationships, and how this approach can help to reveal key mechanistic drivers of global biodiversity patterns. 

***9. 08_Resolution Sensitivity.R:*** Finally, this script tests the sensitivity of various complexity measures to changes in the scale at which they are computed by repeatedly calculating each measure for a series of successively smaller pixel cells selected from the same region of our global DEM. 

### Enjoy!!
