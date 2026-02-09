# sac-snow-uh.R Updates

This directory contains updates to the **SAC-SMA / SNOW17 hydrologic modeling and routing workflow**, extending prior restart functionality to support **snapshot-based model state saving and restoration** for reforecasting and ensemble applications.

These updates build on an earlier enhancement that allowed model states to be saved at the **final timestep**, and now enable **saving and restoring states at user-defined timesteps after spin-up**.

---

## Summary of Capabilities

At a high level, the workflow now supports:

- Warm-start (restart) hydrologic simulations
- Snapshotting internal model states at arbitrary timesteps
- Restarting or branching simulations from saved snapshots
- Preserving original model behavior when restart/snapshot features are not used

---

## Updates to Core Functions

### 1. `sac_snow()`

The `sac_snow()` function supports **state persistence and restoration**:

- SAC-SMA and SNOW17 internal states can be:
  - Loaded at the beginning of a run (warm start)
  - Saved at user-defined timesteps or at the end of a run
- Enables continuation of simulations without loss of hydrologic memory

This functionality supports segmented runs, reforecasts, and ensemble branching workflows.

---

### 2. New Snapshot Utilities: `sac-snow-snapshot-functions.R`

A new helper script provides utilities for **snapshot-based workflows**, including:

- Saving model state objects at arbitrary timesteps
- Restoring saved states for restart or continuation
- Managing restart objects for reforecast-style applications

These functions are intended to simplify snapshot and restart logic outside the core model execution functions.

---

### 3. High-Level Workflow: `sac_snow_uh_lagk_states()`

A high-level orchestration function is provided:

`sac_snow_uh_lagk_states()`

This function executes the full hydrologic and routing model chain with optional restart and snapshot support.

## Model Chain Execution Order

1. SAC-SMA + SNOW17 rainfall–runoff modeling
2. Unit Hydrograph (UH) routing
3. Lag-K routing with upstream tributary inflows
4. Channel loss routing

## Outputs (when enabled)
* Routed flow time series
* SAC-SMA / SNOW17 state snapshots or time series
* UH restart states (qprev)
* Lag-K restart states (C-array)

These outputs enable deterministic continuation or branching of simulations.

# Unit Hydrograph Routing (uh())

##Summary of Updates

* Added support for warm-start UH routing using qprev restart states
* Optionally returns UH restart states for continuation across runs
* Core routing logic and parameters remain unchanged
* Output time-shifting logic has been removed
* UH output is aligned directly with input timestep indices

# Lag-K Routing (lagk())

## Summary of Updates and Limitations

* Added partial support for warm-start Lag-K routing using C-array states
* Optionally returns Lag-K restart data
* Corrected use of upstream tributary inputs
* **Restart behavior for Lag-K is not yet fully reliable** and may introduce
inconsistencies when restarting mid-simulation

Lag-K restart functionality should currently be treated as **experimental**.

