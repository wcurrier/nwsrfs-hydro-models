#!/usr/bin/env Rscript

# example_snapshot_usage.R
# 
# Demonstrates how to use the state snapshot functionality for reforecasting
# 
# This script shows:
# 1. Running a simulation with daily snapshots saved
# 2. Loading a specific snapshot
# 3. Running a reforecast from that snapshot
# 4. Comparing reforecast with original run


suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(ggplot2)
  load_all(".")  # Load your package
})

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(ggplot2)
  load_all(".")  # Load your package
})

# Source the snapshot functions
# source("sac_snapshot_functions.R")

# ============================================================================
# CONFIGURATION
# ============================================================================

basin <- "SAKW1"
dt_hours <- 6

# Load data
data(forcingSAKW1)
data(parsSAKW1)
data(upflowSAKW1)

forcing <- forcingSAKW1
pars <- parsSAKW1
uptribs <- upflowSAKW1

# Create output directories
dir.create("snapshots", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("reforecast_states", showWarnings = FALSE)

# ============================================================================
# EXAMPLE 1: RUN SIMULATION WITH DAILY SNAPSHOTS
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 1: Running simulation with daily snapshots\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

snapshot_file <- "snapshots/state_snapshots_SAKW1.rds"

# Run with snapshots saved daily at 00Z (hour 0)
result_with_snapshots <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours = dt_hours,
  forcing = forcing,
  uptribs = uptribs,
  pars = pars,
  save_snapshots = TRUE,
  snapshot_file = snapshot_file,
  snapshot_interval = "daily",  # Save once per day
  snapshot_hour = 0              # At 00Z (midnight UTC)
)

cat("\nSnapshot information:\n")
print(result_with_snapshots$snapshot_info)

# ============================================================================
# EXAMPLE 2: LIST AVAILABLE SNAPSHOTS
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 2: Listing available snapshots\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

available_dates <- list_snapshot_dates(snapshot_file)
cat("Total snapshots available:", length(available_dates), "\n")
cat("\nFirst 10 snapshot dates:\n")
print(head(available_dates, 10))
cat("\nLast 10 snapshot dates:\n")
print(tail(available_dates, 10))

# ============================================================================
# EXAMPLE 3: LOAD A SPECIFIC SNAPSHOT
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 3: Loading a specific snapshot\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

# Choose a date in the middle of the simulation
reforecast_date <- available_dates[length(available_dates) %/% 2]
cat("Loading snapshot from:", reforecast_date, "\n")

snapshot <- load_snapshot(snapshot_file, reforecast_date)

cat("\nSnapshot structure:\n")
cat("  - Datetime:", format(snapshot$datetime), "\n")
cat("  - Timestep index:", snapshot$timestep_index, "\n")
cat("  - Number of zones:", length(snapshot$sac_snow), "\n")
cat("\nSAC-SMA states for zone 1:\n")
print(unlist(snapshot$sac_snow$zone_1[1:6]))

# ============================================================================
# EXAMPLE 4: RUN REFORECAST FROM SNAPSHOT
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 4: Running reforecast from snapshot\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

# Determine the forcing period for the reforecast
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

reforecast_dt <- as.POSIXct(reforecast_date, tz = "UTC")
reforecast_idx <- which(forcing_dt == reforecast_dt)[1]

cat("Reforecast starting at index:", reforecast_idx, "\n")
cat("Reforecast datetime:", format(forcing_dt[reforecast_idx]), "\n")

cat("starting subsetting\n")

# Extract forcing for reforecast period (from snapshot date to end)
forcing_reforecast <- lapply(forcing, function(df) {
  df[reforecast_idx:nrow(df), ]
})
cat("Subsetted forcing reforecast\n")

uptribs_reforecast <- lapply(uptribs, function(df) {
  df[reforecast_idx:nrow(df), ]
})
cat("Subsetted uptribs reforecast\n")

# Run reforecast using the convenience function
reforecast_result <- run_reforecast_from_snapshot(
  snapshot_file = snapshot_file,
  snapshot_date = reforecast_date,
  forcing = forcing_reforecast,
  uptribs = uptribs_reforecast,
  pars = pars,
  basin_name = basin,
  temp_dir = "reforecast_states"
)

cat("Reforecast simulation complete!\n")
cat("Reforecast length:", length(reforecast_result$flow_cfs), "timesteps\n")

# ============================================================================
# EXAMPLE 5: COMPARE REFORECAST WITH ORIGINAL RUN
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 5: Comparing reforecast with original run\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

# Extract corresponding portion from original run
original_flow <- result_with_snapshots$flow_cfs[reforecast_idx:length(forcing_dt)]
reforecast_flow <- reforecast_result$flow_cfs

# Check lengths match
min_len <- min(length(original_flow), length(reforecast_flow))
original_flow <- original_flow[1:min_len]
reforecast_flow <- reforecast_flow[1:min_len]

# Calculate difference
flow_diff <- reforecast_flow - original_flow
max_diff <- max(abs(flow_diff))
rmse <- sqrt(mean(flow_diff^2))

cat("\nComparison statistics:\n")
cat("  Maximum absolute difference:", max_diff, "cfs\n")
cat("  RMSE:", rmse, "cfs\n")
cat("  Mean difference:", mean(flow_diff), "cfs\n")

# Create comparison plot
reforecast_dt_vec <- forcing_dt[reforecast_idx:(reforecast_idx + min_len - 1)]

plot_df <- data.frame(
  datetime = reforecast_dt_vec,
  original = original_flow,
  reforecast = reforecast_flow,
  difference = flow_diff
)

# Plot first 30 days
plot_window <- plot_df[1:min(nrow(plot_df), 30 * 24 / dt_hours), ]

p1 <- ggplot(plot_window) +
  geom_line(aes(datetime, original, color = "Original run"), size = 0.8) +
  geom_line(aes(datetime, reforecast, color = "Reforecast"), 
            size = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Original run" = "blue", "Reforecast" = "red")) +
  labs(
    title = paste("Reforecast comparison:", basin),
    subtitle = paste("Starting from", reforecast_date),
    x = "Date",
    y = "Flow (cfs)",
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p2 <- ggplot(plot_window) +
  geom_line(aes(datetime, difference), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Difference (Reforecast - Original)",
    x = "Date",
    y = "Difference (cfs)"
  ) +
  theme_minimal()

# Save plots
ggsave(
  "figures/reforecast_comparison.png",
  gridExtra::grid.arrange(p1, p2, ncol = 1),
  width = 12,
  height = 8
)

cat("\nPlot saved to: figures/reforecast_comparison.png\n")

# ============================================================================
# EXAMPLE 6: SAVE SNAPSHOTS AT DIFFERENT INTERVALS
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 6: Saving snapshots at different intervals\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

# Example 6a: Save every 4 timesteps (every 24 hours if dt_hours=6)
cat("6a. Saving every 4 timesteps (24 hours)...\n")
result_4step <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours = dt_hours,
  forcing = forcing,
  uptribs = uptribs,
  pars = pars,
  save_snapshots = TRUE,
  snapshot_file = "snapshots/snapshots_every_24hr_SAKW1.rds",
  snapshot_interval = 4  # Every 4 timesteps
)
cat("Saved", result_4step$snapshot_info$n_snapshots, "snapshots\n")

# Example 6b: Save every timestep (for high-frequency reforecasting)
# WARNING: This can create very large files!
cat("\n6b. Saving every timestep (WARNING: large file)...\n")
cat("Skipping this example to avoid large files.\n")
cat("To enable, uncomment the code below:\n")
cat("
# result_all <- sac_snow_uh_lagk_states_with_snapshots(
#   dt_hours = dt_hours,
#   forcing = forcing,
#   uptribs = uptribs,
#   pars = pars,
#   save_snapshots = TRUE,
#   snapshot_file = 'snapshots/snapshots_all_SAKW1.rds',
#   snapshot_interval = 'all'
# )
\n")

# ============================================================================
# EXAMPLE 7: BATCH REFORECASTING
# ============================================================================

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nEXAMPLE 7: Batch reforecasting (pseudo-operational)\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

cat("Demonstrating how to generate multiple reforecasts...\n\n")

# Get all available snapshot dates
all_snapshots <- list_snapshot_dates(snapshot_file)

# Select every 10th snapshot for reforecasting (to save time)
reforecast_dates <- all_snapshots[seq(10, length(all_snapshots), by = 10)]

cat("Will generate", length(reforecast_dates), "reforecasts\n")
cat("Reforecast dates:\n")
print(reforecast_dates)

# In a real operational setting, you would loop through these dates
# and generate forecasts for each one:
#
# for (fcst_date in reforecast_dates) {
#   
#   # Load appropriate forcing data for forecast period
#   # (In operations, this might be from NWP models)
#   
#   # Run reforecast
#   fcst <- run_reforecast_from_snapshot(
#     snapshot_file = snapshot_file,
#     snapshot_date = fcst_date,
#     forcing = forecast_forcing,
#     uptribs = forecast_uptribs,
#     pars = pars,
#     basin_name = basin,
#     temp_dir = tempdir()
#   )
#   
#   # Save forecast output
#   # saveRDS(fcst, paste0("forecasts/forecast_", fcst_date, ".rds"))
# }

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\nAll examples complete!\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n\n")

cat("Summary of created files:\n")
cat("  - snapshots/state_snapshots_SAKW1.rds (daily snapshots)\n")
cat("  - snapshots/snapshots_every_24hr_SAKW1.rds (4-timestep snapshots)\n")
cat("  - figures/reforecast_comparison.png\n")
cat("  - reforecast_states/ (temporary restart files)\n")
