#!/usr/bin/env Rscript
# ============================================================================
# RESTART DIVERGENCE DIAGNOSTIC
# ============================================================================
# Tests each component individually to isolate restart error source.
# ============================================================================

library(rfchydromodels)
source("../R/sac-snow-snapshot-functions.R")
source("../R/sac-snow-uh.R")

data(forcingSAKW1)
data(parsSAKW1)
data(upflowSAKW1)

cat("======================================================================\n")
cat("RESTART DIVERGENCE DIAGNOSTIC\n")
cat("======================================================================\n\n")

sim_length <- nrow(forcingSAKW1[[1]])
n_zones <- length(forcingSAKW1)
has_uptribs <- !is.null(upflowSAKW1)
restart_step <- 3000

cat("Sim length:", sim_length, "\n")
cat("Zones:", n_zones, "\n")
cat("Has uptribs:", has_uptribs, "\n")
cat("Restart step:", restart_step, "\n\n")

# Helper
compare_flows <- function(label, full_vec, restart_vec, n_show = 20) {
  n <- min(length(full_vec), length(restart_vec))
  full_vec <- full_vec[1:n]
  restart_vec <- restart_vec[1:n]
  diffs <- restart_vec - full_vec
  abs_diffs <- abs(diffs)

  cat(sprintf("  Max abs diff:  %.6f\n", max(abs_diffs)))
  cat(sprintf("  RMSE:          %.6f\n", sqrt(mean(diffs^2))))
  cat(sprintf("  Max diff step: %d\n", which.max(abs_diffs)))

  cat(sprintf("\n  First %d timesteps:\n", min(n_show, n)))
  for (i in 1:min(n_show, n)) {
    cat(sprintf("    Step %3d: full=%12.4f  restart=%12.4f  diff=%12.6f\n",
                i, full_vec[i], restart_vec[i], diffs[i]))
  }

  if (n >= 100) {
    rmse_first50 <- sqrt(mean(diffs[1:50]^2))
    rmse_last50 <- sqrt(mean(tail(diffs, 50)^2))
    cat(sprintf("\n  RMSE first 50: %.6f  |  RMSE last 50: %.6f\n", rmse_first50, rmse_last50))
    if (rmse_first50 > rmse_last50 * 5) {
      cat("  >> Error concentrated at START and decays (spinup artifact)\n")
    }
  }

  return(max(abs_diffs))
}

# Helper: write proper restart CSV combining SAC states + SNOW17 cs + taprev
write_restart_csv_from_snapshot <- function(states_df, cs_ts, taprev_ts, step, n_zones, filepath) {
  # SAC states from states_df at the given timestep
  uztwc <- sapply(1:n_zones, function(z) states_df[[paste0("uztwc_", z)]][step])
  uzfwc <- sapply(1:n_zones, function(z) states_df[[paste0("uzfwc_", z)]][step])
  lztwc <- sapply(1:n_zones, function(z) states_df[[paste0("lztwc_", z)]][step])
  lzfsc <- sapply(1:n_zones, function(z) states_df[[paste0("lzfsc_", z)]][step])
  lzfpc <- sapply(1:n_zones, function(z) states_df[[paste0("lzfpc_", z)]][step])
  adimc <- sapply(1:n_zones, function(z) states_df[[paste0("adimc_", z)]][step])

  # SNOW17 carryover from cs_ts (19 x sim_length x n_zones)
  cs_snapshot <- cs_ts[, step, ]  # 19 x n_zones (or vector if 1 zone)
  if (is.null(dim(cs_snapshot))) cs_snapshot <- matrix(cs_snapshot, ncol = 1)

  # taprev from taprev_ts (sim_length x n_zones)
  taprev <- taprev_ts[step, ]

  cat(sprintf("  Restart at step %d:\n", step))
  for (z in 1:n_zones) {
    cat(sprintf("    Zone %d: uztwc=%.2f lztwc=%.2f swe(cs5)=%.2f taprev=%.2f\n",
                z, uztwc[z], lztwc[z], cs_snapshot[5, z], taprev[z]))
  }

  # Build CSV in the format sac_snow() expects
  df <- data.frame(
    uztwc = uztwc,
    uzfwc = uzfwc,
    lztwc = lztwc,
    lzfsc = lzfsc,
    lzfpc = lzfpc,
    adimc = adimc,
    taprev = taprev
  )
  for (j in 1:19) {
    df[[paste0("cs_", j)]] <- cs_snapshot[j, ]
  }

  write.csv(df, filepath, row.names = FALSE)
  cat("  Wrote restart CSV:", filepath, "\n\n")
}

# ======================================================================
# STEP 1: Full continuous run with state saving
# ======================================================================
cat("--- STEP 1: Full continuous SAC-SMA/SNOW17 run ---\n")

full_sac <- sac_snow(6, forcingSAKW1, parsSAKW1,
                     return_states = TRUE, save_restart = TRUE)

cat("  Output names:", paste(names(full_sac), collapse=", "), "\n")
cat("  cs_ts dim:", paste(dim(full_sac$cs_ts), collapse=" x "), "\n\n")

# Write restart CSV from snapshot
restart_csv <- tempfile(fileext = ".csv")
write_restart_csv_from_snapshot(
  full_sac$states, full_sac$cs_ts, full_sac$taprev_ts,
  restart_step, n_zones, restart_csv
)

# Verify the CSV
cat("  Restart CSV contents:\n")
restart_check <- read.csv(restart_csv)
print(restart_check[, 1:7])
cat("\n")

# ======================================================================
# TEST A: SAC-SMA/SNOW17 alone
# ======================================================================
cat("======================================================================\n")
cat("TEST A: SAC-SMA/SNOW17 alone (isolate model state restart)\n")
cat("======================================================================\n")

# Subset forcing from restart_step+1 onward (restart_step state is the
# END state after processing that timestep, so next input is step+1)
forcing_restart <- lapply(forcingSAKW1, function(x) x[(restart_step + 1):sim_length, ])
n_compare <- sim_length - restart_step

restart_sac <- sac_snow(6, forcing_restart, parsSAKW1,
                        return_states = TRUE, save_restart = FALSE,
                        restart_file = restart_csv)

# Get TCI for comparison
tci_cols <- grep("^tci_", names(full_sac$states), value = TRUE)
full_tci_total <- rowSums(full_sac$states[, tci_cols, drop = FALSE])

if (is.data.frame(restart_sac)) {
  restart_tci_cols <- grep("^tci_", names(restart_sac), value = TRUE)
  restart_tci_total <- rowSums(restart_sac[, restart_tci_cols, drop = FALSE])
} else {
  restart_tci_total <- rowSums(restart_sac)
}

# Compare: full run from step restart_step+1 onward vs restart
full_seg_a <- full_tci_total[(restart_step + 1):sim_length]
restart_seg_a <- restart_tci_total[1:n_compare]

sac_err <- compare_flows("SAC", full_seg_a, restart_seg_a)

if (sac_err < 0.001) {
  cat("\n  >> VERDICT: SAC-SMA/SNOW17 restart is PERFECT\n\n")
} else {
  cat("\n  >> VERDICT: SAC-SMA/SNOW17 restart DIVERGES\n\n")
}

# ======================================================================
# TEST B: UH convolution
# ======================================================================
cat("======================================================================\n")
cat("TEST B: UH convolution (isolate unit hydrograph restart)\n")
cat("======================================================================\n")

# Full UH on continuous TCI
full_tci_matrix <- as.matrix(full_sac$states[, tci_cols])
full_uh <- uh(6, full_tci_matrix, parsSAKW1, sum_zones = TRUE, return_states = TRUE)

cat("  Full UH flow length:", length(full_uh$flow), "\n")

# Build UH restart from TCI history (manual qprev extraction)
# The UH convolution needs the last M timesteps of TCI as its state.
# extract_uh_snapshot expects states_df with tci columns - let's build that.
uh_snapshot <- tryCatch({
  # Build a minimal "states_df" that has tci columns for extract_uh_snapshot
  fake_states <- full_sac$states[, tci_cols, drop = FALSE]
  extract_uh_snapshot(fake_states, restart_step, parsSAKW1)
}, error = function(e) {
  cat("  extract_uh_snapshot error:", e$message, "\n")
  cat("  Building UH restart manually from TCI history...\n")

  # Manual construction: qprev = last M TCI values before restart_step
  m <- 1000
  restart_states <- list()
  for (z in 1:n_zones) {
    shape <- parsSAKW1[parsSAKW1$name == "unit_shape", ]$value[z]
    toc_gis <- parsSAKW1[parsSAKW1$name == "unit_toc", ]$value[z]
    toc_adj <- parsSAKW1[parsSAKW1$name == "unit_toc_adj", ]$value[z]
    if (is.na(toc_gis) | is.na(toc_adj)) {
      scale <- parsSAKW1[parsSAKW1$name == "unit_scale", ]$value[z]
    } else {
      scale <- uh2p_get_scale(shape, toc_gis * toc_adj, 1)
    }

    # qprev: the last m TCI values ending at restart_step
    start_idx <- max(1, restart_step - m + 1)
    tci_hist <- full_tci_matrix[start_idx:restart_step, z]
    # Pad with leading zeros if not enough history
    if (length(tci_hist) < m) {
      tci_hist <- c(rep(0, m - length(tci_hist)), tci_hist)
    }
    # Reverse so qprev[1] = most recent
    restart_states[[z]] <- list(
      qprev = rev(tci_hist),
      shape = shape,
      scale = scale
    )
  }
  restart_states
})

if (!is.null(uh_snapshot)) {
  cat("  UH snapshot zones:", length(uh_snapshot), "\n")
  cat("  Zone 1 qprev non-zero:", sum(uh_snapshot[[1]]$qprev != 0), "\n")
}

# Restart TCI for UH input (from sac restart)
if (is.data.frame(restart_sac)) {
  restart_tci_cols <- grep("^tci_", names(restart_sac), value = TRUE)
  restart_tci_matrix <- as.matrix(restart_sac[, restart_tci_cols])
} else {
  restart_tci_matrix <- restart_sac[1:n_compare, ]
}

# UH WITH restart
restart_uh_with <- tryCatch(
  uh(6, restart_tci_matrix, parsSAKW1, sum_zones = TRUE,
     return_states = FALSE, uh_restart = uh_snapshot),
  error = function(e) { cat("  UH with restart error:", e$message, "\n"); NULL }
)

# UH WITHOUT restart (cold start)
restart_uh_without <- uh(6, restart_tci_matrix, parsSAKW1, sum_zones = TRUE,
                         return_states = FALSE, uh_restart = NULL)

full_uh_seg <- full_uh$flow[(restart_step + 1):sim_length]

get_flow <- function(x) if (is.list(x) && "flow" %in% names(x)) x$flow else x

if (!is.null(restart_uh_with)) {
  cat("\n  WITH UH restart:\n")
  uh_with_err <- compare_flows("UH-with", full_uh_seg, get_flow(restart_uh_with)[1:n_compare])
} else {
  uh_with_err <- NA
  cat("\n  UH restart not available\n")
}

cat("\n  WITHOUT UH restart (cold start):\n")
uh_without_err <- compare_flows("UH-cold", full_uh_seg, get_flow(restart_uh_without)[1:n_compare])

if (!is.na(uh_with_err) && uh_with_err < 0.01) {
  cat("\n  >> VERDICT: UH restart is PERFECT\n\n")
} else if (!is.na(uh_with_err) && uh_with_err < uh_without_err * 0.1) {
  cat("\n  >> VERDICT: UH restart HELPS significantly\n\n")
} else {
  cat("\n  >> VERDICT: UH restart not effective\n\n")
}

# ======================================================================
# TEST C: Lag-K alone
# ======================================================================
if (has_uptribs) {
  cat("======================================================================\n")
  cat("TEST C: Lag-K alone (isolate C-array restart)\n")
  cat("======================================================================\n")

  full_lagk <- lagk(6, upflowSAKW1, parsSAKW1,
                    sum_routes = TRUE, return_states = TRUE,
                    save_states = TRUE)

  full_lagk_flow <- full_lagk$flow
  c_array_ts <- full_lagk$c_array_ts
  cat("  c_array_ts dim:", paste(dim(c_array_ts), collapse=" x "), "\n")

  n_uptribs <- length(upflowSAKW1)
  c_restart <- matrix(c_array_ts[, restart_step, ], nrow = 100, ncol = n_uptribs)
  cat("  c_restart non-zero:", sum(c_restart != 0), "of", length(c_restart), "\n")

  # Restart from step+1
  uptribs_restart <- lapply(upflowSAKW1, function(x) x[(restart_step + 1):sim_length, ])

  restart_lagk <- lagk(6, uptribs_restart, parsSAKW1,
                       sum_routes = TRUE, return_states = FALSE,
                       save_states = FALSE,
                       restart_c_array = c_restart)

  full_seg_c <- full_lagk_flow[(restart_step + 1):sim_length]
  restart_seg_c <- get_flow(restart_lagk)[1:n_compare]

  cat("\n")
  lagk_err <- compare_flows("LagK", full_seg_c, restart_seg_c)

  if (lagk_err < 0.01) {
    cat("\n  >> VERDICT: Lag-K restart is PERFECT\n\n")
  } else {
    cat("\n  >> VERDICT: Lag-K restart DIVERGES\n\n")
  }
} else {
  lagk_err <- 0
}

# ======================================================================
# SUMMARY
# ======================================================================
cat("======================================================================\n")
cat("SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("  A. SAC-SMA/SNOW17 (TCI):  max diff = %.6f\n", sac_err))
if (!is.na(uh_with_err)) {
  cat(sprintf("  B. UH (with restart):     max diff = %.6f\n", uh_with_err))
}
cat(sprintf("     UH (cold start):       max diff = %.6f\n", uh_without_err))
if (has_uptribs) {
  cat(sprintf("  C. Lag-K:                 max diff = %.4f cfs\n", lagk_err))
}

cat("\nDiagnosis:\n")
if (sac_err < 0.001) cat("  * SAC-SMA/SNOW17 restart is PERFECT\n")
if (sac_err >= 0.001) cat("  * SAC-SMA/SNOW17 restart DIVERGES\n")
if (uh_without_err > 1.0) cat("  * UH cold-start causes large error (missing convolution history)\n")
if (!is.na(uh_with_err) && uh_with_err < 0.01) {
  cat("  * UH restart FIXES the UH error\n")
} else if (!is.na(uh_with_err) && uh_with_err < uh_without_err * 0.5) {
  cat("  * UH restart helps but doesn't fully fix error\n")
}
if (has_uptribs && lagk_err < 0.01) cat("  * Lag-K restart is PERFECT\n")
if (has_uptribs && lagk_err >= 0.01) cat("  * Lag-K restart DIVERGES\n")
cat("\n")

# Cleanup
file.remove(restart_csv)