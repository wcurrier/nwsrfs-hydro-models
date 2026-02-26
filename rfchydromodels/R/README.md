# sac-snow-uh.R Updates

This directory contains updates to the SAC-SMA / SNOW17 hydrologic modeling and routing workflow, extending prior restart functionality to support snapshot-based model state saving and restoration for reforecasting and ensemble applications.

These updates build on an earlier enhancement that allowed model states to be saved at the final timestep, and now enable saving and restoring states at user-defined timesteps after spin-up.

## Summary of Capabilities

At a high level, the workflow now supports:

* Warm-start (restart) hydrologic simulations
* Snapshotting internal model states at arbitrary timesteps
* Restarting or branching simulations from saved snapshots
* Running simulations from minimal forcing inputs (precipitation and temperature only)
* Preserving original model behavior when restart/snapshot features are not used

## Updates to Core Functions

### 1. `sac_snow()`

The `sac_snow()` function supports state persistence and restoration:

* SAC-SMA and SNOW17 internal states can be:
  * Loaded at the beginning of a run (warm start)
  * Saved at user-defined timesteps or at the end of a run
* Enables continuation of simulations without loss of hydrologic memory

This functionality supports segmented runs, reforecasts, and ensemble branching workflows.

### 2. New Snapshot Utilities: `sac-snow-snapshot-functions.R`

A new helper script provides utilities for snapshot-based workflows, including:

* Saving model state objects at arbitrary timesteps
* Restoring saved states for restart or continuation
* Managing restart objects for reforecast-style applications

These functions are intended to simplify snapshot and restart logic outside the core model execution functions.

### 3. Forcing Preparation: `prepare_forcing.R`

A new module for computing derived forcing fields from minimal inputs. This enables simulations to run from forcing files containing only precipitation (`map_mm`) and temperature (`mat_degc`).

#### Functions

* **`prepare_forcing(forcing, pars, dt_hours, ae_tbl)`** — Entry point. Computes missing `ptps`, `pet_mm`, and `etd_mm` columns from minimal forcing. Requires an area-elevation table (`ae_tbl`). Called automatically within the main simulation functions when `ae_tbl` is provided, or can be called independently for offline forcing preparation.

* **`compute_pet_etd(forcing, pars, dt_hours)`** — Internal function. Aggregates sub-daily temperature to daily tave/tmax/tmin, computes daily PET via the Hargreaves-Samani equation, distributes evenly across sub-daily timesteps, and applies monthly `peadj` crop coefficients to produce `etd_mm`.

* **`pet_hs_daily(lat, jday, tave, tmax, tmin)`** — Hargreaves-Samani daily PET (mm/day). Negative values are clamped to zero.

* **`interpolate_peadj(peadj_monthly, year, month, day_of_month, dt_hours)`** — Interpolates monthly `peadj_01`…`peadj_12` values between months, centered on the 16th of each month, matching the Fortran `fa_ts` interpolation scheme. Handles leap years.

#### How ptps is computed

The fraction of precipitation as snow (`ptps`) is computed by the existing `rsnwelev()` function, which calls the Fortran EX42 subroutine. EX42 calculates a rain-snow line elevation from temperature, lapse rate, and a threshold temperature, then maps it to a snow fraction using an area-elevation curve.

**Important**: The `talr` (lapse rate) parameter must be **positive** when passed to `rsnwelev`. EX42 computes `RSL = TAELEV + (PTA - PXTEMP) * (100 / TALR)`. A positive `talr` means temperature decreases with elevation (standard atmospheric lapse rate). A negative `talr` inverts the rain-snow partitioning, producing incorrect ptps values (e.g., all snow at warm temperatures).

#### How etd is computed

Evapotranspiration demand is computed as `etd = pet * peadj`, where:
* `pet` is the Hargreaves-Samani PET computed from raw `mat_degc`
* `peadj` is the monthly crop coefficient interpolated between months centered on the 16th

The interpolation matches the Fortran `fa_ts` logic exactly. However, since `prepare_forcing` uses raw temperature rather than forcing-adjusted temperature (`mat_fa`), PET/ETD values will differ slightly from those produced by the full FA pipeline. This is an expected and intentional simplification for workflows that bypass the forcing adjustment step.

### 4. High-Level Workflow: `sac_snow_uh_lagk_states()`

A high-level orchestration function that executes the full hydrologic and routing model chain with optional restart, snapshot, and forcing preparation support.

#### Model Chain Execution Order

1. Forcing preparation (if `ae_tbl` provided): compute missing `ptps`, `pet_mm`, `etd_mm`
2. SAC-SMA + SNOW17 rainfall–runoff modeling
3. Unit Hydrograph (UH) routing
4. Lag-K routing with upstream tributary inflows
5. Channel loss routing

#### Outputs (when enabled)

* Routed flow time series
* SAC-SMA / SNOW17 state snapshots or time series
* UH restart states (`qprev`)
* Lag-K restart states (C-array)

These outputs enable deterministic continuation or branching of simulations.

### 5. Reforecast: `run_reforecast_from_snapshot()`

Updated to accept an `ae_tbl` parameter, which is passed through to `sac_snow_uh_lagk_states()` for automatic forcing preparation during reforecast runs.

## Unit Hydrograph Routing (`uh()`)

### Summary of Updates

* Added support for warm-start UH routing using `qprev` restart states
* Optionally returns UH restart states for continuation across runs
* Core routing logic and parameters remain unchanged
* Output time-shifting logic has been removed
* UH output is aligned directly with input timestep indices

## Lag-K Routing (`lagk()`)

### Summary of Updates and Limitations

* Added partial support for warm-start Lag-K routing using C-array states
* Optionally returns Lag-K restart data
* Corrected use of upstream tributary inputs
* Restart behavior for Lag-K is not yet fully reliable and may introduce inconsistencies when restarting mid-simulation
* Lag-K restart functionality should currently be treated as experimental.

Please let me know if you have any questions.
