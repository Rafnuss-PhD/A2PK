# Area-to-point kriging (A2PK)
[![DOI](https://zenodo.org/badge/105257986.svg)](https://zenodo.org/badge/latestdoi/105257986)

This MATLAS package allows you to generate fine-scale stocastic Gaussian simulation based on a known coarse scale image of the same domain. This is done through area-to-point kriging, a co-kriging estimation where the coase-scale is considered as a secondary variable so that the cross-covariance can be computed easily. 

## Area-to-point kriging

The main function is ``A2PK.m``, and an example script is provided in [``script.mlx``](https://rafnuss-phd.github.io/A2PK/script).

### Examples
The following two scripts give a step-by-step instruction to run area-to-point kriging:
- [``A2PK _gaussian.mlx``](https://rafnuss-phd.github.io/A2PK/examples/A2PK_gaussian) is the easiest place to start. It will guide you through the computation of a simulation based on a synthetic example.
- [``A2PK_cond_gaussian.mlx``](https://rafnuss-phd.github.io/A2PK/examples/A2PK_cond_gaussian) expends the previous example with a conditional case.

### Application to Eletrical Resistivity Tomography (ERT)
During my PhD, I applied A2PK to the simulation of a fine-scale scale electrical conductivity field based on the smooth result of a deterministic inversion of an eletrical resistivity tomography (ERT). The code used for this can be found in [the ERT folder](https://github.com/Rafnuss-PhD/A2PK/tree/master/ERT), and, in particular, the published script [``script_ERT.mlx``](https://Rafnuss-phd.github.io/A2PK/ERT/script_ERT). 

This work has lead to a paper which is currently under revision. I'll provide the link as soons as it's available. 

### Application to Hydraulic tomography (HT)
We also start working on an adaptation of the similar technique for hydraulic tomography, refers to [the HT folder](https://github.com/Rafnuss-PhD/A2PK/tree/master/HT), and more precisely to [the script used to generate the dataset and the result](https://github.com/Rafnuss-PhD/A2PK/blob/master/HT/script_elec_cond.m). Despite promising result, I did not carried on this work as I was finishing my PhD.


## References
- Nussbaumer, R., Mariethoz, G., Linde N., & Holliger, K. (2018). Simulation of fine-scale electrical conductivity fields using tomograms and area-to-point kriging, *GeoEnv2018, Belfast*. DOI: [10.13140/RG.2.2.31712.12801](https://www.doi.org/10.13140/RG.2.2.31712.12801)
- Kyriakidis, Phaedon C. 2004. â€œA Geostatistical Framework for Area-to-Point Spatial Interpolation.*Geographical Analysis 36(3):259*. (http://doi.wiley.com/10.1353/geo.2004.0009).
