# ============================================================================
# CONSOLIDATED STATE SNAPSHOT FUNCTIONALITY FOR REFORECASTING
# ============================================================================
#
# This file contains all the snapshot-related functions needed for saving
# and restoring model states for reforecasting applications.
#
# Key design decisions:
# 1. Single uh() function - no separate uh_with_timeseries()
# 2. UH snapshots reconstructed from TCI history (no Fortran changes needed)
# 3. Lag-K snapshots obtained from c_array_ts saved by Fortran at each timestep
# 4. All snapshots stored in a single RDS file per basin
# 5. Snapshots are saved incrementally to avoid memory exhaustion
#
# ============================================================================


# ============================================================================
# SNAPSHOT EXTRACTION FUNCTIONS
# ============================================================================

#' Extract UH snapshot at a specific timestep
#'
#' Reconstructs the qprev array needed for UH warm start from the
#' TCI time series. The qprev array contains the last M TCI values
#' that are used in the convolution.
#'
#' @param states_df data frame of model states (must contain tci_* columns)
#' @param idx timestep index to extract snapshot at
#' @param pars parameter data frame
#' @param m UH convolution length (default 1000)
#'
#' @return list of UH states for each zone
#' @export
extract_uh_snapshot <- function(states_df, idx, pars, m = 1000) {

  # Count zones from tci columns
  tci_cols <- grep("^tci_", names(states_df), value = TRUE)
  n_zones <- length(tci_cols)

  if (n_zones == 0) {
    stop("No tci_* columns found in states_df")
  }

  uh_snapshot <- list()

  for (z in 1:n_zones) {
    # Get shape and scale parameters
    shape <- pars[pars$name == "unit_shape", ]$value[z]
    toc_gis <- pars[pars$name == "unit_toc", ]$value[z]
    toc_adj <- pars[pars$name == "unit_toc_adj", ]$value[z]

    if (is.na(toc_gis) | is.na(toc_adj)) {
      scale <- pars[pars$name == "unit_scale", ]$value[z]
    } else {
      toc <- toc_gis * toc_adj
      scale <- uh2p_get_scale(shape, toc, 1)
    }

    # Extract TCI values for this zone
    tci_col <- paste0("tci_", z)
    tci_vec <- states_df[[tci_col]]

    # Build qprev: the last M TCI values ending at timestep idx
    # qprev[1] = oldest value (idx - M + 1)
    # qprev[M] = newest value (idx)
    qprev <- numeric(m)

    for (j in 1:m) {
      tci_idx <- idx - m + j
      if (tci_idx >= 1 && tci_idx <= length(tci_vec)) {
        qprev[j] <- tci_vec[tci_idx]
      } else {
        qprev[j] <- 0.0  # Pad with zeros for indices before start
      }
    }

    uh_snapshot[[z]] <- list(
      qprev = qprev,
      shape = shape,
      scale = scale
    )
  }

  return(uh_snapshot)
}


#' Extract SAC-SMA and SNOW17 states at a specific timestep
#'
#' This function extracts the exact model states at a given timestep index
#' using the actual cs array from the Fortran output (cs_ts), not reconstructed values.
#'
#' @param states_df data frame of model states (from sac_snow with return_states=TRUE)
#' @param idx timestep index to extract snapshot at
#' @param n_zones number of zones
#' @param cs_ts 3D array of cs values at each timestep (19 x sim_length x n_zones)
#' @param taprev_ts matrix of taprev values at each timestep (sim_length x n_zones)
#'
#' @return list of SAC/SNOW states for each zone
#' @export
extract_sac_snow_snapshot <- function(states_df, idx, n_zones, cs_ts, taprev_ts) {

  snapshot <- list()

  for (z in 1:n_zones) {
    zone_snapshot <- list(
      # SAC-SMA states - directly from states_df
      uztwc = states_df[[paste0("uztwc_", z)]][idx],
      uzfwc = states_df[[paste0("uzfwc_", z)]][idx],
      lztwc = states_df[[paste0("lztwc_", z)]][idx],
      lzfsc = states_df[[paste0("lzfsc_", z)]][idx],
      lzfpc = states_df[[paste0("lzfpc_", z)]][idx],
      adimc = states_df[[paste0("adimc_", z)]][idx],

      # SNOW17 simple states - directly from states_df
      swe = states_df[[paste0("swe_", z)]][idx],
      neghs = states_df[[paste0("neghs_", z)]][idx],
      liqw = states_df[[paste0("liqw_", z)]][idx],

      # Temperature carryover - from taprev_ts (EXACT value from Fortran)
      taprev = taprev_ts[idx, z],

      # SNOW17 cs array - EXACT values from Fortran cs_ts output
      # cs_ts is dimensioned as (19, sim_length, n_zones)
      cs = cs_ts[, idx, z]
    )

    snapshot[[paste0("zone_", z)]] <- zone_snapshot
  }

  return(snapshot)
}


#' Extract Lag-K C-array snapshot at a specific timestep
#'
#' Extracts the exact Lag-K C-array state from the c_array_ts 3D array
#' that was saved by Fortran at every timestep during the simulation.
#'
#' @param c_array_ts 3D array of C-array values (100 x sim_length x n_uptribs)
#' @param idx timestep index
#' @param n_uptribs number of upstream tributaries
#'
#' @return matrix of C-array values (100 x n_uptribs) at the specified timestep
#' @export
extract_lagk_snapshot <- function(c_array_ts, idx, n_uptribs) {
  # c_array_ts is dimensioned as (100, sim_length, n_uptribs)
  # Extract the C-array at timestep idx for all uptribs
  c_array_at_idx <- matrix(0, nrow = 100, ncol = n_uptribs)

  for (u in 1:n_uptribs) {
    c_array_at_idx[, u] <- c_array_ts[, idx, u]
  }

  return(c_array_at_idx)
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Determine which timestep indices should have snapshots saved
#'
#' @param forcing_dt POSIXct vector of simulation datetimes
#' @param snapshot_interval integer (every N timesteps), "daily", or "all"
#' @param snapshot_hour hour of day for daily snapshots (0-23)
#' @param dt_hours model timestep in hours
#'
#' @return integer vector of indices to save
#' @export
determine_snapshot_indices <- function(
    forcing_dt,
    snapshot_interval,
    snapshot_hour = 0,
    dt_hours = 6
) {

  if (is.character(snapshot_interval) && snapshot_interval == "daily") {
    # Save at specific hour each day
    hour_of_day <- as.integer(format(forcing_dt, "%H"))

    # Find indices matching the target hour
    indices <- which(hour_of_day == snapshot_hour)

    # If target hour doesn't exist in data, find closest
    if (length(indices) == 0) {
      warning(
        sprintf(
          "Hour %d not found in data. Using closest available hour.",
          snapshot_hour
        )
      )

      # Find the hour closest to target for each day
      dates <- as.Date(forcing_dt)
      unique_dates <- unique(dates)

      indices <- sapply(unique_dates, function(d) {
        day_indices <- which(dates == d)
        day_hours <- hour_of_day[day_indices]
        closest_idx <- day_indices[which.min(abs(day_hours - snapshot_hour))]
        return(closest_idx)
      })
    }

    return(sort(indices))

  } else if (is.numeric(snapshot_interval)) {
    # Save every N timesteps
    indices <- seq(snapshot_interval, length(forcing_dt), by = snapshot_interval)
    return(indices)

  } else if (snapshot_interval == "all") {
    # Save every timestep
    return(seq_along(forcing_dt))

  } else {
    stop("snapshot_interval must be 'daily', 'all', or a positive integer")
  }
}


# ============================================================================
# MAIN SNAPSHOT FUNCTION
# ============================================================================

#' Run SAC-SMA/SNOW17/UH/Lag-K with state snapshots at specified intervals
#'
#' This is the main function for running simulations with snapshot saving.
#' It runs the complete model chain and saves exact model states at
#' user-specified intervals for reforecasting applications.
#'
#' Snapshots are written to disk in chunks to avoid memory exhaustion.
#'
#' @param dt_hours timestep in hours
#' @param forcing list of forcing data frames (one per zone)
#' @param uptribs list of upstream tributary flow data frames (or NULL)
#' @param pars parameter data frame
#' @param restart_file optional SAC-SMA/SNOW17 restart CSV
#' @param lagk_restart_file optional Lag-K restart RDS file
#' @param uh_restart_file optional UH restart RDS file
#' @param debug_components logical; include intermediate flow components in output
#' @param save_snapshots logical; whether to save state snapshots
#' @param snapshot_file path to save snapshots (RDS file)
#' @param snapshot_interval "daily", "all", or integer (every N timesteps)
#' @param snapshot_hour integer; if snapshot_interval="daily", save at this hour (0-23)
#' @param snapshot_chunk_size integer; number of snapshots to accumulate before
#'        flushing to disk (default 2000). Smaller values use less memory.
#'
#' @return A list containing:
#' \describe{
#'   \item{flow_cfs}{Total routed flow after channel loss}
#'   \item{states}{Data frame of SAC-SMA/SNOW17 state time series}
#'   \item{restart}{SAC-SMA/SNOW17 restart states}
#'   \item{lagk_restart}{Lag-K restart C-array}
#'   \item{uh_restart}{UH restart states}
#'   \item{snapshot_info}{Information about saved snapshots}
#' }
#' @export
sac_snow_uh_lagk_states_with_snapshots <- function(
    dt_hours,
    forcing,
    uptribs = NULL,
    pars,
    restart_file = NULL,
    lagk_restart_file = NULL,
    uh_restart_file = NULL,
    debug_components = FALSE,
    save_snapshots = FALSE,
    snapshot_file = NULL,
    snapshot_interval = "daily",
    snapshot_hour = 0,
    snapshot_chunk_size = 2000,
    ae_tbl = NULL) {

  # ---- Input validation ----
  if (save_snapshots && is.null(snapshot_file)) {
    stop("snapshot_file must be provided when save_snapshots = TRUE")
  }

  if (save_snapshots) {
    snapshot_dir <- dirname(snapshot_file)
    if (snapshot_dir != "." && !dir.exists(snapshot_dir)) {
      dir.create(snapshot_dir, recursive = TRUE)
      cat("Created snapshot directory:", snapshot_dir, "\n")
    }
  }

  # Auto-apply rsnwelev if pxtemp is in pars and ae_tbl provided
  if (!is.null(ae_tbl) && any(pars$name == "pxtemp")) {
    forcing <- rsnwelev(forcing, pars, ae_tbl)
  }
    
  # ---- Build datetime vector ----
  forcing_dt <- as.POSIXct(
    ISOdatetime(
      forcing[[1]]$year,
      forcing[[1]]$month,
      forcing[[1]]$day,
      forcing[[1]]$hour,
      0, 0,
      tz = "UTC"
    ),
    tz = "UTC"
  )

  sim_length <- length(forcing_dt)
  n_zones <- length(forcing)

  cat("=== Simulation Configuration ===\n")
  cat("  Timesteps:", sim_length, "\n")
  cat("  Period:", format(min(forcing_dt)), "to", format(max(forcing_dt)), "\n")
  cat("  Zones:", n_zones, "\n")
  cat("  Upstream tributaries:", if(is.null(uptribs)) 0 else length(uptribs), "\n")

  # ---- Determine snapshot indices ----
  save_indices <- NULL
  if (save_snapshots) {
    save_indices <- determine_snapshot_indices(
      forcing_dt,
      snapshot_interval,
      snapshot_hour,
      dt_hours
    )
    cat("  Snapshots to save:", length(save_indices), "\n")
  }
  cat("\n")

  # ---- SAC-SMA + SNOW17 ----
  cat("=== Running SAC-SMA + SNOW17 ===\n")

  # When saving snapshots, we need save_restart=TRUE to get cs_ts and taprev_ts
  sac_out <- sac_snow(
    dt_hours,
    forcing,
    pars,
    return_states = TRUE,
    save_restart = save_snapshots,  # Enable to get cs_ts/taprev_ts when saving snapshots
    restart_file = restart_file
  )

  # Extract outputs based on whether save_restart was enabled
  if (save_snapshots) {
    # sac_out is a list with states, restart, cs_ts, taprev_ts
    states_df <- sac_out$states
    restart_vals <- sac_out$restart
    cs_ts <- sac_out$cs_ts       # 3D array: 19 x sim_length x n_zones
    taprev_ts <- sac_out$taprev_ts  # matrix: sim_length x n_zones
  } else {
    # sac_out is just the states data frame
    states_df <- sac_out
    restart_vals <- NULL
    cs_ts <- NULL
    taprev_ts <- NULL
  }

  cat("  Complete.\n\n")

  # ---- Extract TCI ----
  tci_cols <- grep("^tci_", names(states_df), value = TRUE)
  if (length(tci_cols) == 0) {
    stop("No tci_* columns found in sac_snow output")
  }
  tci <- as.matrix(states_df[, tci_cols, drop = FALSE])

  # ---- UH routing ----
  cat("=== Running UH Routing ===\n")
  uh_restart <- NULL
  if (!is.null(uh_restart_file) && file.exists(uh_restart_file)) {
    cat("  Loading UH restart from:", uh_restart_file, "\n")
    uh_restart <- readRDS(uh_restart_file)
  }

  uh_result <- uh(
    dt_hours,
    tci,
    pars,
    return_states = TRUE,
    uh_restart = uh_restart
  )

  flow_uh <- uh_result$flow
  uh_restart_out <- uh_result$restart
  cat("  Complete.\n\n")

  # ---- Lag-K routing ----
  lagk_flow <- 0
  lagk_restart <- NULL
  c_array_ts <- NULL
  n_uptribs <- 0

  if (!is.null(uptribs) && length(uptribs) > 0) {
    cat("=== Running Lag-K Routing ===\n")
    n_uptribs <- length(uptribs)

    # Load initial C-array if provided
    initial_c_array <- NULL
    if (!is.null(lagk_restart_file) && file.exists(lagk_restart_file)) {
      cat("  Loading Lag-K restart from:", lagk_restart_file, "\n")
      initial_c_array <- readRDS(lagk_restart_file)
    }

    # Run Lag-K with save_states=TRUE if we need snapshots
    lagk_result <- lagk(
      dt_hours,
      uptribs,
      pars,
      sum_routes = TRUE,
      return_states = TRUE,
      save_states = save_snapshots,  # Save c_array at every timestep for snapshots
      restart_c_array = initial_c_array
    )

    if (save_snapshots) {
      # lagk_result is a list with flow, c_array, c_array_ts
      lagk_flow <- lagk_result$flow
      lagk_restart <- lagk_result$c_array
      c_array_ts <- lagk_result$c_array_ts  # 3D: 100 x sim_length x n_uptribs
    } else {
      # lagk_result is a data frame with c_array attribute
      routed_cols <- grep("^routed_", names(lagk_result), value = TRUE)
      if (length(routed_cols) > 0) {
        lagk_flow <- rowSums(as.matrix(lagk_result[, routed_cols, drop = FALSE]))
      }
      lagk_restart <- attr(lagk_result, "c_array")
    }

    cat("  Complete.\n\n")
  } else {
    cat("=== Skipping Lag-K (no upstream tributaries) ===\n\n")
  }

  # ---- Channel loss ----
  cat("=== Applying Channel Loss ===\n")
  total_flow_cfs <- chanloss(
    flow_uh + lagk_flow,
    forcing,
    dt_hours,
    pars
  )
  cat("  Complete.\n\n")

  # ---- Build and save snapshots (chunked to manage memory) ----
  snapshot_info <- NULL

  if (save_snapshots && length(save_indices) > 0) {
    cat("=== Saving State Snapshots ===\n")

    # Force a GC + small allocation to surface any heap corruption
    # from the Fortran calls BEFORE entering the snapshot loop.
    # If the heap is corrupted, this will crash here with a clear context.
    gc(verbose = FALSE)
    .heap_test <- vector("list", 100)
    rm(.heap_test)
    cat("  Heap integrity check passed.\n")

    n_snapshots <- length(save_indices)

    # --- Chunked snapshot saving ---
    # Instead of building one giant list in memory (which caused heap corruption),
    # we process snapshots in chunks and merge them at the end.
    #
    # Strategy:
    #   1. Process save_indices in chunks of snapshot_chunk_size
    #   2. Each chunk produces a temporary RDS file
    #   3. After all chunks, merge into final file
    #   4. Clean up temp files

    chunk_files <- character(0)
    chunk_start <- 1

    while (chunk_start <= n_snapshots) {
      chunk_end <- min(chunk_start + snapshot_chunk_size - 1, n_snapshots)
      chunk_indices <- save_indices[chunk_start:chunk_end]

      # Build snapshot list for this chunk
      chunk_snapshots <- vector("list", length(chunk_indices))
      chunk_names <- character(length(chunk_indices))

      for (j in seq_along(chunk_indices)) {
        idx <- chunk_indices[j]
        timestamp <- format(forcing_dt[idx], "%Y-%m-%d %H:%M:%S")
        chunk_names[j] <- timestamp

        # SAC/SNOW snapshot - using exact cs_ts values from Fortran
        sac_snapshot <- extract_sac_snow_snapshot(
          states_df,
          idx,
          n_zones,
          cs_ts,
          taprev_ts
        )

        # UH snapshot (reconstructed from TCI - this is exact)
        uh_snapshot <- extract_uh_snapshot(
          states_df,
          idx,
          pars
        )

        # Lag-K snapshot - using exact c_array_ts values from Fortran
        lagk_snapshot <- NULL
        if (!is.null(c_array_ts) && n_uptribs > 0) {
          lagk_snapshot <- extract_lagk_snapshot(c_array_ts, idx, n_uptribs)
        }

        # Combine
        chunk_snapshots[[j]] <- list(
          datetime = forcing_dt[idx],
          timestep_index = idx,
          sac_snow = sac_snapshot,
          uh = uh_snapshot,
          lagk = lagk_snapshot
        )
      }

      names(chunk_snapshots) <- chunk_names

      # Save chunk to temp file
      chunk_file <- tempfile(
        pattern = sprintf("snapshot_chunk_%05d_", chunk_start),
        fileext = ".rds"
      )
      saveRDS(chunk_snapshots, chunk_file)
      chunk_files <- c(chunk_files, chunk_file)

      # Progress
      cat(sprintf("  Saved chunk %d-%d / %d (%.1f%%)\n",
                  chunk_start, chunk_end, n_snapshots,
                  100 * chunk_end / n_snapshots))

      # Free chunk memory
      rm(chunk_snapshots, chunk_names)

      chunk_start <- chunk_end + 1
    }

    # Merge all chunks into final file
    cat("  Merging", length(chunk_files), "chunks into final snapshot file...\n")
    all_snapshots <- list()
    for (cf in chunk_files) {
      chunk_data <- readRDS(cf)
      # Append to master list
      for (nm in names(chunk_data)) {
        all_snapshots[[nm]] <- chunk_data[[nm]]
      }
      rm(chunk_data)
      file.remove(cf)  # Clean up temp file
    }

    # Save final merged file
    saveRDS(all_snapshots, snapshot_file)
    file_size_mb <- file.info(snapshot_file)$size / 1024 / 1024

    cat(sprintf("  Saved %d snapshots to %s (%.2f MB)\n",
                length(all_snapshots), snapshot_file, file_size_mb))

    snapshot_info <- list(
      file = snapshot_file,
      n_snapshots = length(all_snapshots),
      dates = names(all_snapshots),
      first_date = names(all_snapshots)[1],
      last_date = names(all_snapshots)[length(all_snapshots)],
      file_size_mb = file_size_mb
    )

    rm(all_snapshots)
  }

  cat("\n=== Simulation Complete ===\n")

  # ---- Free large snapshot-only arrays before building result ----
  if (!is.null(cs_ts)) rm(cs_ts)
  if (!is.null(taprev_ts)) rm(taprev_ts)
  if (!is.null(c_array_ts)) rm(c_array_ts)

  # ---- Build result ----
  result <- list(
    flow_cfs = total_flow_cfs,
    states = states_df
  )

  if (debug_components) {
    result$flow_uh <- flow_uh
    result$flow_lagk <- lagk_flow
    result$flow_before_chanloss <- flow_uh + lagk_flow
  }

  if (!is.null(restart_vals)) result$restart <- restart_vals
  if (!is.null(lagk_restart)) result$lagk_restart <- lagk_restart
  if (!is.null(uh_restart_out)) result$uh_restart <- uh_restart_out
  if (!is.null(snapshot_info)) result$snapshot_info <- snapshot_info

  return(result)
}


# ============================================================================
# SNAPSHOT LOADING AND CONVERSION FUNCTIONS
# ============================================================================

#' Load a specific snapshot from the snapshot file
#'
#' @param snapshot_file path to snapshot RDS file
#' @param snapshot_date character or POSIXct; the datetime to load
#'
#' @return snapshot list with sac_snow, uh, and lagk states
#' @export
load_snapshot <- function(snapshot_file, snapshot_date) {

  if (!file.exists(snapshot_file)) {
    stop("Snapshot file not found: ", snapshot_file)
  }

  snapshots <- readRDS(snapshot_file)

  # Convert date to character format used in snapshots
  if (inherits(snapshot_date, "POSIXct") || inherits(snapshot_date, "POSIXt")) {
    snapshot_key <- format(snapshot_date, "%Y-%m-%d %H:%M:%S")
  } else {
    snapshot_key <- as.character(snapshot_date)
  }

  if (!snapshot_key %in% names(snapshots)) {
    available_dates <- names(snapshots)
    stop(
      sprintf(
        "Snapshot date '%s' not found.\nFirst 5 available: %s\nLast 5 available: %s\nTotal: %d",
        snapshot_key,
        paste(head(available_dates, 5), collapse = ", "),
        paste(tail(available_dates, 5), collapse = ", "),
        length(available_dates)
      )
    )
  }

  return(snapshots[[snapshot_key]])
}


#' List all available snapshot dates in a snapshot file
#'
#' @param snapshot_file path to snapshot RDS file
#'
#' @return character vector of available snapshot datetimes
#' @export
list_snapshot_dates <- function(snapshot_file) {

  if (!file.exists(snapshot_file)) {
    stop("Snapshot file not found: ", snapshot_file)
  }

  snapshots <- readRDS(snapshot_file)
  return(names(snapshots))
}


#' Convert a snapshot to restart file format
#'
#' @param snapshot snapshot list from load_snapshot()
#' @param output_dir directory to save restart files
#' @param basin_name basin identifier for file naming
#'
#' @return list of created file paths
#' @export
snapshot_to_restart_files <- function(snapshot, output_dir, basin_name) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  created_files <- list()

  # ---- SAC/SNOW restart CSV ----
  n_zones <- length(snapshot$sac_snow)
  cat("Converting snapshot with", n_zones, "zones to restart files\n")

  restart_df <- data.frame(
    zone = 1:n_zones,
    uztwc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$uztwc),
    uzfwc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$uzfwc),
    lztwc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$lztwc),
    lzfsc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$lzfsc),
    lzfpc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$lzfpc),
    adimc = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$adimc),
    taprev = sapply(1:n_zones, function(z) snapshot$sac_snow[[z]]$taprev)
  )

  # Add SNOW17 cs array columns
  if (!is.null(snapshot$sac_snow[[1]]$cs)) {
    for (i in 1:19) {
      restart_df[[paste0("cs_", i)]] <- sapply(
        1:n_zones,
        function(z) snapshot$sac_snow[[z]]$cs[i]
      )
    }
  }

  restart_csv <- file.path(output_dir, paste0("restart_", basin_name, ".csv"))
  write.csv(restart_df, restart_csv, row.names = FALSE)
  created_files$sac_snow <- restart_csv
  cat("  Created SAC/SNOW restart:", restart_csv, "\n")

  # ---- UH restart RDS ----
  if (!is.null(snapshot$uh)) {
    uh_rds <- file.path(output_dir, paste0("uh_restart_", basin_name, ".rds"))
    saveRDS(snapshot$uh, uh_rds)
    created_files$uh <- uh_rds
    cat("  Created UH restart:", uh_rds, "\n")
  }

  # ---- Lag-K restart RDS ----
  if (!is.null(snapshot$lagk)) {
    lagk_rds <- file.path(output_dir, paste0("lagk_restart_", basin_name, ".rds"))
    saveRDS(snapshot$lagk, lagk_rds)
    created_files$lagk <- lagk_rds
    cat("  Created Lag-K restart:", lagk_rds, "\n")
  } else {
    cat("  NOTE: No Lag-K state in snapshot (may be NULL if no upstream tribs)\n")
  }

  return(created_files)
}


#' Run a reforecast from a saved snapshot
#'
#' @param snapshot_file path to snapshot RDS file
#' @param snapshot_date datetime to load snapshot from
#' @param forcing list of forcing data frames for forecast period
#' @param uptribs list of upstream tributary flows (or NULL)
#' @param pars parameter data frame
#' @param basin_name basin identifier
#' @param dt_hours timestep in hours (default 6)
#' @param temp_dir directory for temporary restart files
#'
#' @return results from sac_snow_uh_lagk_states
#' @export
run_reforecast_from_snapshot <- function(
    snapshot_file,
    snapshot_date,
    forcing,
    uptribs = NULL,
    pars,
    basin_name,
    dt_hours = 6,
    temp_dir = tempdir()
) {

  cat("=== Running Reforecast ===\n")
  cat("  Loading snapshot from:", snapshot_date, "\n")
  snapshot <- load_snapshot(snapshot_file, snapshot_date)

  cat("  Converting snapshot to restart files...\n")
  restart_files <- snapshot_to_restart_files(snapshot, temp_dir, basin_name)

  cat("  Running forecast simulation...\n")
  result <- sac_snow_uh_lagk_states(
    dt_hours = dt_hours,
    forcing = forcing,
    uptribs = uptribs,
    pars = pars,
    restart_file = restart_files$sac_snow,
    uh_restart_file = restart_files$uh,
    lagk_restart_file = restart_files$lagk
  )

  cat("  Forecast complete!\n")
  return(result)
}


#' Get snapshot summary statistics
#'
#' @param snapshot_file path to snapshot RDS file
#'
#' @return data frame with summary statistics
#' @export
snapshot_summary <- function(snapshot_file) {

  if (!file.exists(snapshot_file)) {
    stop("Snapshot file not found: ", snapshot_file)
  }

  snapshots <- readRDS(snapshot_file)

  dates <- names(snapshots)
  n_snapshots <- length(dates)
  first_snap <- snapshots[[1]]
  n_zones <- length(first_snap$sac_snow)
  file_size_mb <- file.info(snapshot_file)$size / 1024 / 1024

  # Check for Lag-K states
  has_lagk <- !is.null(first_snap$lagk)
  n_with_lagk <- sum(sapply(snapshots, function(s) !is.null(s$lagk)))

  cat("=== Snapshot File Summary ===\n")
  cat("  File:", snapshot_file, "\n")
  cat("  Size:", sprintf("%.2f MB", file_size_mb), "\n")
  cat("  Number of snapshots:", n_snapshots, "\n")
  cat("  Number of zones:", n_zones, "\n")
  cat("  Date range:", dates[1], "to", dates[n_snapshots], "\n")
  cat("  Has UH states:", !is.null(first_snap$uh), "\n")
  cat("  Snapshots with Lag-K states:", n_with_lagk, "/", n_snapshots, "\n")

  invisible(list(
    file = snapshot_file,
    file_size_mb = file_size_mb,
    n_snapshots = n_snapshots,
    n_zones = n_zones,
    first_date = dates[1],
    last_date = dates[n_snapshots],
    has_uh = !is.null(first_snap$uh),
    n_with_lagk = n_with_lagk
  ))
}


