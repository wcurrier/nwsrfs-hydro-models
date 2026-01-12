## This Fork
  
This fork includes the following updates from the original project:

### Overarching Goal
- Enable SAC-SMA/SNOW17 to have **restart (warm-state) capability**, allowing Unit Hydrograph (UH) and Lag-K routing to continue from where they left off.
- **Modified files**:
  - `R/sac-snow-uh.R`: added restart functionality (highlighted changes are in `R/README.md`).
  - Fortran source code: `sac-snow.f90`, `lagk_run.f90`, `duamel.f` (see comments in each file for changes).
- **Spin-up functionality** is maintained if no restart data are provided.
- Currently, only saves all states (SAC-SMA/SNOW17, Lag-K, UH) at the last time step.
  *Future work*: optionally save all time steps for more detailed hindcast evaluation.

### Testing
- `tests/test_all.R` was used to validate the code by splitting the model halfway through simulation and checking if the simulated streamflow matched.
- Individual tests were run on SAC-SMA/SNOW17, UH, and Lag-K components. Minor differences seem to occasionally appear.
- Examples include full run (no restart/spin-up) and restart run (split at midpoints and restarted). Reproducible with `rfchydromodels/tests/test_all.R`.

### Notes
- Only tested in R using the provided example datasets.
- Python version (`py-rfchydromodels`) not tested successfully due to Unit Hydrograph issues.
- Images or plots can be included using R Markdown code chunks, for example:

### Example Figures

#### Basin SAKW1
![Restart run for basin SAKW1](https://raw.githubusercontent.com/wcurrier/nwsrfs-hydro-models/main/rfchydromodels/tests/figures/restart_zoom_week_SAKW1.png)  
**Figure 1.** Example of SAC-SMA/SNOW17 restart run for basin SAKW1 – Includes Lag-K routing.

#### Basin WCHW1
![Restart run for basin WCHW1](https://raw.githubusercontent.com/wcurrier/nwsrfs-hydro-models/main/rfchydromodels/tests/figures/restart_zoom_week_WCHW1.png)  
**Figure 2.** Example of SAC-SMA/SNOW17 restart run for basin WCHW1 – Does not include Lag-K routing.


# NWRFC Operational Hydrology Models 

## Overview

The Northwest River Forecast Center (NWRFC) utilizes the National Weather Service River Forecasting System (NWSRFS) to provide timely information related to flooding, water supply, drought, recreation, navigation, and environmental flows. Originally developed in the late 1970s, NWSRFS remains a core component of the NWS Community Hydrologic Prediction System (CHPS). The system includes a suite of models that simulate soil moisture, snow accumulation and melt, flow routing, channel loss, and consumptive water use. For additional details on each model, see [this link](https://www.weather.gov/owp/oh_hrl_nwsrfs_users_manual_htm_xrfsdocpdf) .

To support hydrologic model calibration and development, NWRFC has created FORTRAN 90 wrappers to execute the original NWSRFS FORTRAN 77 source code. This repository contains the original FORTRAN 77 model code which has been verified to be is functionally equivalent to the Java-based implementation used in CHPS. The wrapped suite of available models includes SAC-SMA, SNOW17, Unit Hydrograph, LAGK, CHANLOSS, and CONS_USE.

Also included in this repository are Python and R packages that compile and interact with the FORTRAN 90 wrappers. These tools are intended to facilitate coupling the hydrologic models with modern optimization packages, supporting model calibration and evaluation.

**Languages:** R, Python, FORTRAN 77, and FORTRAN 90

**Compiler:** A FORTRAN compiler is required to install this package. This package has been tested with [gfortran](https://gcc.gnu.org/wiki/GFortran). See [this page](https://cran.r-project.org/bin/macosx/tools/) for a simple installation option on macOS

**Known OS Compatibility:** macOS and Red Hat OS (will likley work on any modern linux distro). Windows compatibility through WSL. 

**Time Step Compatibility:** This package and its wrappers have been tested only with a 6-hour time step. Use with other time steps may require additional configuration or validation.

## Installation

### R Package Installation

Install the R package from within R using the following command:

```R
devtools::install_github('NOAA-NWRFC/nwsrfs-hydro-models',subdir='rfchydromodels')
```   

or from the command line:

```bash
git clone https://github.com/NOAA-NWRFC/nwsrfs-hydro-models.git
cd nwsrfs-hydro-models
R CMD INSTALL rfchydromodels
```

See the documentation `?rfchydromodels` and `?sac_snow_uh` for examples. 

### Python Package Installation

**Tested Python Version:** 3.10.3\
**Package Dependencies:**  numpy, pandas\
**Dependencies:** numpy, pandas
numpy's `f2py` is used to compile the source code and FORTRAN wrappers. To compile the FORTRAN source:

```bash
git clone https://github.com/NOAA-NWRFC/nwsrfs-hydro-models.git
cd nwsrfs-hydro-models/py-rfchydromodels/utilities
make
```
See `nwsrfs-hydro-models/py-rfchydromodels/run_example.py` for example code demonstrating how to execute the NWSRFS models.

*Note:  An equivalent Python version of the R package is planned for a future release of this repository.*

## Credits and References

Please cite the following work when using this tool:

Walters, G., Bracken, C., et al., "A comprehensive calibration framework for the Northwest River Forecast Center." Unpublished manuscript, Submitted 2025, [Preprint](https://eartharxiv.org/repository/view/8993/)

If adapting this code, please credit this repository as the original source. 

### NWSRFS References

* Burnash, Robert J. C., et al. A generalized streamflow simulation system : conceptual modeling for digital computers. , National Weather Service, 1973
* Anderson, Eric. Snow Accumulation and Ablation Model. National Oceanic and Atmospheric Administration, 2006
* Linsley, R.K., et al. Hydrology for Engineers, McGraw-Hill series in water resources and environmental engineering. McGraw-Hill, 1982
* NOAA. Consumptive Use Operation. National Oceanic and Atmospheric Administration, 2005

## Acknowledgment

Guidance on compiling and running NWSRFS code was informed by work from Andy Wood ([andywood@ucar.edu](mailto:andywood@ucar.edu)) and collaborators. See: [NWS_hydro_models](https://github.com/NCAR/NWS_hydro_models/) GitHub repository

## Legal Disclaimer

This is a scientific product and does not represent official communication from NOAA or the U.S. Department of Commerce. All code is provided "as is."

See full disclaimer: [NOAA GitHub Policy](https://github.com/NOAAGov/Information)
 \
 \
 \
<img src="https://www.weather.gov/bundles/templating/images/header/header.png" alt="NWS-NOAA Banner">

[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [National Weather Service](https://www.weather.gov/) | [Northwest River Forecast Center](https://www.nwrfc.noaa.gov/rfc/)
