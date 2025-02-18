# resistance

This project runs the analysis associated with insecticide
susceptibility. The analysis has 5 main steps:

1.  Read in some real insecticide susceptibility bioassay data, a nice
    clean subset of the the bioassay data from Moyes at al. (2019)
    <https://doi.org/10.5061/dryad.dn4676s>, restricted to Kenya,
    pyrethroids, and all records of members of the Anopheles gambiae
    complex together.

2.  Read in some environmental covariates for Kenya, and plot the IR
    data over it.

3.  Crop the spatial areas down to the region with sufficient IR data

4.  Fit and plot a spatial-only model

5.  Fit and plot a model with covariates and a spatial smooth
