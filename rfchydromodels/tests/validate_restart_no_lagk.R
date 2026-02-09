#!/usr/bin/env Rscript
# ============================================================================
# RESTART VALIDATION: SAC-SMA/SNOW17 + UH + Channel Loss (no Lag-K)
# ============================================================================
# 1. Run a full cold-start continuous simulation, saving daily snapshots.
# 2. Pick 4 random dates (after spinup, before end).
# 3. Warm-start from each snapshot and run to end of forcing period.
# 4. Compare warm-start flow to cold-start flow in 4-panel figures.
# ============================================================================

library(rfchydromodels)
source("../R/sac-snow-snapshot-functions.R")
source("../R/sac-snow-uh.R")

set.seed(42)  # reproducible random dates

data(forcingSAKW1)
data(parsSAKW1)

dt_hours <- 6
sim_length <- nrow(forcingSAKW1[[1]])
n_zones <- length(forcingSAKW1)

# Build datetime vector
forcing_dt <- as.POSIXct(
  ISOdatetime(
    forcingSAKW1[[1]]$year,
    forcingSAKW1[[1]]$month,
    forcingSAKW1[[1]]$day,
    forcingSAKW1[[1]]$hour,
    0, 0, tz = "UTC"
  ), tz = "UTC"
)

cat("======================================================================\n")
cat("STEP 1: Full cold-start continuous simulation with daily snapshots\n")
cat("======================================================================\n")
cat("  Timesteps:", sim_length, "\n")
cat("  Period:", format(min(forcing_dt)), "to", format(max(forcing_dt)), "\n")
cat("  Zones:", n_zones, "\n\n")

snapshot_file <- "snapshots_validation.rds"

full_run <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours = dt_hours,
  forcing = forcingSAKW1,
  uptribs = NULL,           # No Lag-K
  pars = parsSAKW1,
  save_snapshots = TRUE,
  snapshot_file = snapshot_file,
  snapshot_interval = "daily",
  snapshot_hour = 0,
  snapshot_chunk_size = 2000
)

full_flow <- full_run$flow_cfs
cat("\n  Full run flow range:", sprintf("%.2f to %.2f cfs\n", min(full_flow), max(full_flow)))

# ======================================================================
# STEP 2: Select 4 random snapshot dates (after 1-year spinup)
# ======================================================================
cat("\n======================================================================\n")
cat("STEP 2: Select 4 random snapshot dates for warm-start tests\n")
cat("======================================================================\n")

available_dates <- list_snapshot_dates(snapshot_file)
available_dt <- as.POSIXct(available_dates, tz = "UTC")

# Spinup: skip first 365 days. Buffer: stop 180 days before end.
spinup_end <- min(forcing_dt) + 365 * 86400
buffer_start <- max(forcing_dt) - 180 * 86400

valid_mask <- available_dt > spinup_end & available_dt < buffer_start
valid_dates <- available_dates[valid_mask]
valid_dt <- available_dt[valid_mask]

cat("  Available snapshots:", length(available_dates), "\n")
cat("  Valid (post-spinup, pre-buffer):", length(valid_dates), "\n")

# Pick 4 random dates
pick_idx <- sort(sample(length(valid_dates), 4))
test_dates <- valid_dates[pick_idx]
test_dt <- valid_dt[pick_idx]

for (i in 1:4) {
  cat(sprintf("  Test %d: %s\n", i, test_dates[i]))
}

# ======================================================================
# STEP 3: Run warm-start reforecasts from each snapshot
# ======================================================================
cat("\n======================================================================\n")
cat("STEP 3: Run warm-start reforecasts\n")
cat("======================================================================\n")

reforecast_results <- list()

for (i in 1:4) {
  cat(sprintf("\n--- Reforecast %d from %s ---\n", i, test_dates[i]))

  # Find the timestep index for this snapshot date
  snap_idx <- which(forcing_dt == test_dt[i])
  if (length(snap_idx) == 0) {
    # Find closest
    snap_idx <- which.min(abs(forcing_dt - test_dt[i]))
  }
  cat("  Snapshot timestep index:", snap_idx, "\n")

  # Subset forcing from snap_idx+1 onward
  forcing_forecast <- lapply(forcingSAKW1, function(x) x[(snap_idx + 1):sim_length, ])
  cat("  Forecast timesteps:", nrow(forcing_forecast[[1]]), "\n")

  # Run reforecast using the snapshot machinery
  result <- run_reforecast_from_snapshot(
    snapshot_file = snapshot_file,
    snapshot_date = test_dates[i],
    forcing = forcing_forecast,
    uptribs = NULL,   # No Lag-K
    pars = parsSAKW1,
    basin_name = "SAKW1",
    dt_hours = dt_hours
  )

  reforecast_results[[i]] <- list(
    date = test_dates[i],
    dt = test_dt[i],
    snap_idx = snap_idx,
    flow = result$flow_cfs,
    n_steps = length(result$flow_cfs)
  )

  # Quick comparison
  full_seg <- full_flow[(snap_idx + 1):sim_length]
  restart_seg <- result$flow_cfs[1:length(full_seg)]
  max_diff <- max(abs(full_seg - restart_seg))
  rmse <- sqrt(mean((full_seg - restart_seg)^2))
  cat(sprintf("  Max diff: %.6f cfs  |  RMSE: %.6f cfs\n", max_diff, rmse))
}

# ======================================================================
# STEP 4: Create comparison figures
# ======================================================================
cat("\n======================================================================\n")
cat("STEP 4: Creating comparison figures\n")
cat("======================================================================\n")

if (!dir.exists("figures")) dir.create("figures")

# --- Figure 1: 4-panel comparison (flow + difference) ---
png("figures/restart_validation_4panel.png", width = 1800, height = 2400, res = 150)
par(mfrow = c(4, 2), mar = c(3, 4, 2, 1), oma = c(2, 0, 2, 0))

for (i in 1:4) {
  r <- reforecast_results[[i]]
  snap_idx <- r$snap_idx

  # Time vectors
  full_time <- forcing_dt[(snap_idx + 1):sim_length]
  full_seg <- full_flow[(snap_idx + 1):sim_length]
  restart_seg <- r$flow[1:length(full_seg)]
  diff_seg <- restart_seg - full_seg

  # Include a bit of context before the restart point
  context_start <- max(1, snap_idx - 4 * 30)  # ~30 days before
  context_time <- forcing_dt[context_start:sim_length]
  context_flow <- full_flow[context_start:sim_length]

  # Limit display to ~1 year after restart for readability
  display_end <- min(length(full_seg), 4 * 365)  # ~1 year of 6-hr steps

  # Left panel: Flow comparison
  plot_time <- full_time[1:display_end]
  plot_full <- full_seg[1:display_end]
  plot_restart <- restart_seg[1:display_end]

  # Also show context
  ctx_end_idx <- snap_idx - context_start + 1
  ctx_time <- forcing_dt[context_start:snap_idx]
  ctx_flow <- full_flow[context_start:snap_idx]

  ylim <- range(c(plot_full, plot_restart, ctx_flow), na.rm = TRUE)

  plot(ctx_time, ctx_flow, type = "l", col = "gray70",
       xlim = range(c(ctx_time, plot_time)),
       ylim = ylim,
       xlab = "", ylab = "Flow (cfs)",
       main = sprintf("Restart from %s", format(r$dt, "%Y-%m-%d")))
  lines(plot_time, plot_full, col = "black", lwd = 1.5)
  lines(plot_time, plot_restart, col = "red", lwd = 1, lty = 2)
  abline(v = as.numeric(r$dt), col = "blue", lty = 3, lwd = 2)
  legend("topright",
         legend = c("Cold start (reference)", "Warm start", "Restart point"),
         col = c("black", "red", "blue"),
         lty = c(1, 2, 3), lwd = c(1.5, 1, 2),
         cex = 0.7, bg = "white")

  # Right panel: Difference
  plot_diff <- diff_seg[1:display_end]
  plot(plot_time, plot_diff, type = "l", col = "darkred",
       xlab = "", ylab = "Diff (cfs)",
       main = sprintf("Difference (max=%.4f, RMSE=%.4f)",
                       max(abs(diff_seg)), sqrt(mean(diff_seg^2))))
  abline(h = 0, col = "gray50", lty = 2)
  abline(v = as.numeric(r$dt), col = "blue", lty = 3, lwd = 2)
}

mtext("Restart Validation: SAC-SMA/SNOW17 + UH (no Lag-K)",
      outer = TRUE, line = 0.5, cex = 1.2, font = 2)
dev.off()

cat("  Saved: figures/restart_validation_4panel.png\n")

# --- Figure 2: Summary statistics ---
png("figures/restart_validation_summary.png", width = 1200, height = 800, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in 1:4) {
  r <- reforecast_results[[i]]
  snap_idx <- r$snap_idx

  full_seg <- full_flow[(snap_idx + 1):sim_length]
  restart_seg <- r$flow[1:length(full_seg)]
  diff_seg <- restart_seg - full_seg

  # Rolling RMSE (window = 4 steps per day * 7 days = 28 steps)
  window <- 28
  n <- length(diff_seg)
  if (n > window) {
    rolling_rmse <- numeric(n - window + 1)
    for (j in 1:(n - window + 1)) {
      rolling_rmse[j] <- sqrt(mean(diff_seg[j:(j + window - 1)]^2))
    }
    days_after <- (1:length(rolling_rmse)) * dt_hours / 24

    plot(days_after, rolling_rmse, type = "l", col = "darkblue",
         xlab = "Days after restart", ylab = "7-day rolling RMSE (cfs)",
         main = sprintf("Restart %s", format(r$dt, "%Y-%m-%d")),
         log = ifelse(max(rolling_rmse) > 0.01, "y", ""))
    abline(h = 0.01, col = "green", lty = 2)
    legend("topright", "0.01 cfs threshold", col = "green", lty = 2, cex = 0.7)
  }
}

dev.off()
cat("  Saved: figures/restart_validation_summary.png\n")

# ======================================================================
# SUMMARY TABLE
# ======================================================================
cat("\n======================================================================\n")
cat("SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("%-22s  %12s  %12s  %12s\n",
            "Restart Date", "Max Diff", "RMSE", "Steps"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (i in 1:4) {
  r <- reforecast_results[[i]]
  snap_idx <- r$snap_idx
  full_seg <- full_flow[(snap_idx + 1):sim_length]
  restart_seg <- r$flow[1:length(full_seg)]
  diff_seg <- restart_seg - full_seg

  cat(sprintf("%-22s  %12.6f  %12.6f  %12d\n",
              format(r$dt, "%Y-%m-%d %H:%M"),
              max(abs(diff_seg)),
              sqrt(mean(diff_seg^2)),
              length(diff_seg)))
}

cat("\nDone. Check figures/ for plots.\n")

# Cleanup
file.remove(snapshot_file)
