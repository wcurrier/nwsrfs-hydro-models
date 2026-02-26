
# This Fork

This fork extends the original project to support restartable and snapshot-based hydrologic modeling workflows, with a focus on reforecasting and ensemble applications.

## Overarching Goal

* Build upon the prior update that enabled saving model states at the final timestep.
* Extend this capability to allow:
   * Saving model states at user-defined timesteps
   * Restarting (warm starting) the model at arbitrary points after spin-up
   * Running simulations from minimal forcing inputs (precipitation and temperature only)
* Preserve original model behavior when restart/snapshot functionality is not used.

At a high level, this fork enables a workflow where a model can be:
1. Run once through spin-up,
2. Snapshotted at one or more timesteps, and
3. Restarted cleanly from those snapshots for continuation, branching, or reforecasting.

Additionally, forcing data can be provided in a minimal format (year, month, day, hour, map_mm, mat_degc) and the remaining required fields (ptps, pet_mm, etd_mm) are computed automatically.

## Scope of Updates

Updates span both model source code and R interfaces:

* Fortran source updates:
   * `sac_snow.f90` – SNOW17 and SAC-SMA restart and snapshot support
   * `duamel.f` – Unit Hydrograph restart capability
   * `lagk_run.f90` – partial restart and snapshot support
   * `flag7.f` – supporting changes for restart workflows

* R interface updates:
   * `R/sac-snow-uh.R` – extended to support snapshot saving, warm starts, and automatic forcing preparation
   * `R/sac-snow-snapshot-functions.R` – new helper functions for saving, restoring, and managing model state snapshots for reforecasting workflows
   * `R/prepare_forcing.R` – new functions for computing derived forcing fields from minimal inputs (see below)

These changes are designed to be backward compatible: if restart/snapshot features are not enabled, model behavior is equivalent to prior versions.

## Forcing Preparation (`prepare_forcing`)

A new `prepare_forcing()` function allows simulations to be run from minimal forcing files containing only precipitation and temperature. The following fields are computed automatically when missing:

* **ptps** (fraction of precipitation as snow): Computed via the `rsnwelev()` Fortran routine using the rain-snow elevation line method. Requires an area-elevation table (`ae_tbl`) and the `pxtemp`/`talr` parameters.
* **pet_mm** (potential evapotranspiration): Computed via the Hargreaves-Samani equation using sub-daily temperature aggregated to daily tave/tmax/tmin.
* **etd_mm** (evapotranspiration demand): Computed as `pet_mm * peadj`, where `peadj` is the monthly crop coefficient (`peadj_01`…`peadj_12`) interpolated between months centered on the 16th, matching the Fortran `fa_ts` interpolation scheme.

### Minimal forcing format

| Column   | Description                |
|----------|----------------------------|
| year     | Year                       |
| month    | Month (1–12)               |
| day      | Day of month               |
| hour     | Hour (0, 6, 12, 18 for 6h) |
| map_mm   | Mean areal precipitation   |
| mat_degc | Mean areal temperature (°C)|

### Required parameters for forcing preparation

| Parameter      | Source | Description |
|----------------|--------|-------------|
| pxtemp         | pars   | Rain/snow threshold temperature (°C) |
| talr           | pars   | Lapse rate (°C/100m, **positive** = temp decreases with elevation) |
| elev           | pars   | Reference elevation for temperature time series (m) |
| alat           | pars   | Latitude (decimal degrees) |
| peadj_01–12    | pars   | Monthly ET crop coefficients per zone |
| ae_tbl         | input  | Area-elevation curve (quantile + elevation per zone) |

### Integration

`prepare_forcing()` is called automatically within `sac_snow_uh_lagk_states_with_snapshots()` and `sac_snow_uh_lagk_states()` when the `ae_tbl` argument is provided. It can also be called independently for offline forcing preparation.

### Parameter sign conventions

* **talr** must be **positive** when passed to `rsnwelev`. The EX42 Fortran subroutine computes the rain-snow line elevation as `RSL = TAELEV + (PTA - PXTEMP) * (100 / TALR)`. A negative talr inverts the rain-snow partitioning.

## Restart and Snapshot Notes

* SAC-SMA and SNOW17 can be restarted cleanly from saved states.
* Unit Hydrograph routing supports warm starts via saved flow history.
* Lag-K routing does not yet restart cleanly in all cases and may introduce small inconsistencies. This component should be treated as experimental for restart-based workflows.

## Tests and Examples

Several test and diagnostic scripts are included to demonstrate and validate behavior:

* `tests/validate_restart_no_lagk.R` – Demonstrates how a model run can be split, restarted, and continued using saved states (excluding Lag-K). This is the recommended reference example.
* `tests/diagnose_restart_divergence.R` – Diagnostic tools for identifying restart-related divergence.
* `tests/example_snapshot_usage.R` – Illustrates snapshot saving and restoration for reforecast-style workflows.
* `tests/test_prepare_forcing.R` – Validates that minimal forcing (auto-computed ptps/pet/etd) produces identical results to pre-computed forcing through the full simulation and reforecast pipeline.

## Current Limitations

* Restart functionality has been tested primarily through the R interface.
* Lag-K restart behavior is incomplete and may introduce errors.
* Python interfaces have not yet been fully validated with the new snapshot and restart features.
* The Hargreaves-Samani PET computed by `prepare_forcing()` uses raw `mat_degc` rather than forcing-adjusted temperature (`mat_fa`), so PET/ETD values will differ slightly from those produced by the full `fa_ts` Fortran pipeline. This is an expected simplification.
* WARNING: For the test periods, results are functionally equivalent to the original code when `uh()` is run with `start_of_timestep = FALSE` and `backfill = FALSE`. These options control whether output flows are shifted by one timestep to account for forcing data labeled at the beginning of the timestep. The correct convention for this behavior is unclear. In this fork, the timestep-shifting logic has been removed, which may be incorrect.

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
