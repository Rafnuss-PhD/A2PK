# A2PK
[![DOI](https://zenodo.org/badge/105257986.svg)](https://zenodo.org/badge/latestdoi/105257986)

This MATLAS package allows you to generate stocastic Gaussin simulation at a fine-scale based on a coarse scale image of the same domain. This is done through area-to-point kriging, a co-kriging estimation where the coase-scale is view as a secondary variable which allows the construction of the cross-covariance easily. 

## Area-to-point kriging

The main function to run area-to-point kriging is ``A2PK.m``. The script file [``script.mlx``](https://rafnuss-phd.github.io/A2PK/script) shows you how to use the function.


### Examples
The scripts below are generated with Matlab live-script and allow your to see step by step what is A2PK and how to use it.
- [``A2PK _gaussian.mlx``](https://rafnuss-phd.github.io/A2PK/examples/A2PK_gaussian) Live Script is the easiest place to start as this guide you through the computation of the simulation.
- [``A2PK_cond_gaussian.mlx``](https://rafnuss-phd.github.io/A2PK/examples/A2PK_cond_gaussian) expends the previous Live Script with the conditional case.

### Application to Eletrical Resistivity Tomography (ERT)
We applied A2PK to the simulation of a fine-scale scale electrical conductivity field based on the smooth result of an inversion based on ERT. Refers to [the ERT folder](https://github.com/Rafnuss-PhD/A2PK/tree/master/ERT) for all the codes, and more precisely to the published script [``script_ERT.mlx``](https://Rafnuss-phd.github.io/A2PK/ERT/script_ERT). 

This work has lead to a paper which is currently under revision. Contact me if you would like a pdf of the draft. 

### Application to Hydraulic tomography (HT)
Refers to [the HT folder](https://github.com/Rafnuss-PhD/A2PK/tree/master/HT), and more precisely to [the script used to generate the dataset and the result](https://github.com/Rafnuss-PhD/A2PK/blob/master/HT/script_elec_cond.m). I'll try to make an html page of this code soon. 
Despite promising result, we did not carried on this reasearch because of lack of time.

## Reference

- Nussbaumer, R., Mariethoz, G., Linde N., & Holliger, K. (2018). Simulation of fine-scale electrical conductivity fields using tomograms and area-to-point kriging, *GeoEnv2018, Belfast*. DOI: [10.13140/RG.2.2.31712.12801](https://www.doi.org/10.13140/RG.2.2.31712.12801)
- Kyriakidis, Phaedon C. 2004. â€œA Geostatistical Framework for Area-to-Point Spatial Interpolation.*Geographical Analysis 36(3):259*. (http://doi.wiley.com/10.1353/geo.2004.0009).
