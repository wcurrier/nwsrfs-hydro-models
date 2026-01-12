# SAC-SNOW-UH Routing Updates

This directory contains updates to the **SAC-SMA / SNOW17 hydrologic modeling and routing workflow**, with a focus on adding **restart (warm-state) capability** and improving the robustness of **Unit Hydrograph (UH)** and **Lag-K routing**.

This outlines the changes to the sac-snow-uh.R functions

---

## Summary of Changes

### 1. Updates to `sac_snow()`

The `sac_snow()` function has been updated to support **model state persistence**:

* SAC-SMA and SNOW17 model states can now be:

  * Loaded at the beginning of a run (if provided)
  * Saved at the end of a run for restart capability
* Enables warm-start hydrologic simulations across segmented or restart-based runs

This allows hydrologic simulations to be paused and resumed without loss of internal model state.

---

### 2. New Function: `sac_snow_uh_lagk_states()`

A new high-level workflow function has been added:

```r
sac_snow_uh_lagk_states()
```

This function executes the full hydrologic and routing model chain with optional warm-state (restart) support and returns both routed flows and restart objects.

#### Model Chain Execution Order

1. **SAC-SMA + SNOW17** rainfall–runoff modeling
2. **Unit Hydrograph (UH)** routing of SAC-SMA/SNOW17 output
3. **Lag-K routing** with upstream tributary inflows
4. **Channel loss routing** using:

   * UH-routed flow
   * Lag-K routed flow

#### Outputs

* Routed flow time series
* SAC-SMA / SNOW17 state time series
* UH restart states (`qprev`)
* Lag-K restart states (C-array)

These outputs enable **exact continuation of simulations across multiple runs** without loss of routing history.

---

## Unit Hydrograph Routing (`uh()`) Updates

### Changes Relative to Previous `uh()`

* Added support for **warm-starting UH routing** via `qprev` restart states
* Optionally returns UH restart states for continuation across runs
* Core routing behavior and parameter handling remain unchanged
* Removed output time shifting for start-of-timestep forcing
* UH output is now aligned directly with the input timestep index

---

## Lag-K Routing (`lagk()`) Updates

### Changes Relative to Previous `lagk()`

* Added support for **warm-starting Lag-K routing** via C-array restart states
* Optionally returns Lag-K restart data for continuation across runs
* Fixed use of `uptribs` instead of the previously undefined `upflow`
* Core routing behavior and parameter handling remain unchanged

---

**Author:** WRC
**Date:** January 2026

