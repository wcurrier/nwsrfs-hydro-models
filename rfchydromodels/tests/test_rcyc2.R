#!/usr/bin/env Rscript
# test_RCYC2.R
#
# 3-zone RCYC2 watershed test (no Lag-K):
#   1. Run POR simulation with daily snapshots at 12Z
#   2. Run reforecast from 2000-12-27 using MM forcing
#   3. Compare reforecast against continuous run and observed flow
#   4. Plot results

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
  library(ggplot2)
  load_all(".")
})

# ============================================================================
# CONFIGURATION
# ============================================================================

basin    <- "RCYC2"
dt_hours <- 6
n_zones  <- 3
data_dir <- "/home/wcurrier/SacSnow"
out_dir  <- "test_RCYC2_output"
snap_dir <- file.path(out_dir, "snapshots")
fig_dir  <- file.path(out_dir, "figures")

for (d in c(out_dir, snap_dir, fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

zones <- paste0(basin, "-", 1:n_zones)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("=== Loading data ===\n")

# Parameters
pars <- read.csv(file.path(data_dir, "parsRCYC2.csv"), stringsAsFactors = FALSE)
cat("  Parameters:", nrow(pars), "rows\n")

# Area-elevation curve
ae_tbl <- read.csv(file.path(data_dir, "area_elev_curve_RCYC2.csv"), check.names = FALSE)
cat("  Area-elevation curve:", nrow(ae_tbl), "rows,", ncol(ae_tbl) - 1, "zones\n")

# Observed daily streamflow
obs_file <- "/home/wcurrier/src/nwsrfs-hydro-autocalibration/runs/CBRFC/RCYC2/flow_daily_RCYC2.csv"
obs_raw <- fread(obs_file)
obs_daily <- data.frame(
  date     = as.Date(ISOdate(obs_raw$year, obs_raw$month, obs_raw$day)),
  flow_cfs = obs_raw$flow_cfs
)
cat("  Observed flow:", nrow(obs_daily), "days (",
    format(min(obs_daily$date)), "to", format(max(obs_daily$date)), ")\n")

# POR forcing (minimal: year, month, day, hour, map_mm, mat_degc)
forcing <- list()
for (z in zones) {
  f <- fread(file.path(data_dir, paste0("forcing_por_", z, ".csv")))
  forcing[[z]] <- as.data.frame(f)
  cat("  Forcing", z, ":", nrow(f), "timesteps\n")
}

# Reforecast forcing (different column names: precip_mm, t2m_degc)
forcing_fcst <- list()
for (z in zones) {
  f <- fread(file.path(data_dir, paste0("forcing_20001227_MM_", z, ".csv")))
  # Rename to standard column names
  names(f)[names(f) == "precip_mm"] <- "map_mm"
  names(f)[names(f) == "t2m_degc"]  <- "mat_degc"
  forcing_fcst[[z]] <- as.data.frame(f)
  cat("  Forecast forcing", z, ":", nrow(f), "timesteps\n")
}

cat("\n")

# ============================================================================
# BUILD DATETIME VECTORS (used throughout)
# ============================================================================

forcing_dt <- as.POSIXct(
  ISOdatetime(forcing[[1]]$year, forcing[[1]]$month, forcing[[1]]$day,
              forcing[[1]]$hour, 0, 0, tz = "UTC"),
  tz = "UTC"
)

refc_dt <- as.POSIXct(
  ISOdatetime(forcing_fcst[[1]]$year, forcing_fcst[[1]]$month,
              forcing_fcst[[1]]$day, forcing_fcst[[1]]$hour,
              0, 0, tz = "UTC"),
  tz = "UTC"
)

# ============================================================================
# STEP 1: RUN POR SIMULATION WITH DAILY SNAPSHOTS AT 12Z
# ============================================================================

cat("=== Step 1: POR simulation with daily snapshots at 12Z ===\n")

snap_file <- file.path(snap_dir, "snapshots_RCYC2.rds")

result_por <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours        = dt_hours,
  forcing         = forcing,
  uptribs         = NULL,
  pars            = pars,
  ae_tbl          = ae_tbl,
  save_snapshots  = TRUE,
  snapshot_file   = snap_file,
  snapshot_interval = "daily",
  snapshot_hour   = 12
)

cat("  POR run complete.\n")
cat("  Snapshots saved:", result_por$snapshot_info$n_snapshots, "\n")
cat("  Flow range:", round(min(result_por$flow_cfs), 2), "to",
    round(max(result_por$flow_cfs), 2), "cfs\n\n")

# Aggregate POR simulation to daily mean for comparison with obs
por_daily <- data.frame(
  date     = as.Date(forcing_dt),
  flow_cfs = result_por$flow_cfs
)
por_daily <- aggregate(flow_cfs ~ date, data = por_daily, FUN = mean)
names(por_daily)[2] <- "sim_cfs"

# Merge with observed
por_vs_obs <- merge(por_daily, obs_daily, by = "date", suffixes = c("_sim", "_obs"))
por_vs_obs$diff <- por_vs_obs$sim_cfs - por_vs_obs$flow_cfs

cat("  POR vs Obs comparison:\n")
cat("    Overlapping days:", nrow(por_vs_obs), "\n")
nse_por <- 1 - sum((por_vs_obs$sim_cfs - por_vs_obs$flow_cfs)^2) /
               sum((por_vs_obs$flow_cfs - mean(por_vs_obs$flow_cfs))^2)
cat("    NSE:", round(nse_por, 4), "\n")
cat("    Bias:", round(mean(por_vs_obs$diff), 2), "cfs\n\n")

# ============================================================================
# STEP 2: REFORECAST FROM 2000-12-27 12:00 UTC
# ============================================================================

cat("=== Step 2: Reforecast from 2000-12-27 ===\n")

# Snapshot date — snapshots are at 12Z
snapshot_date <- "2000-12-27 12:00:00"

# Verify snapshot exists
available <- list_snapshot_dates(snap_file)
if (!snapshot_date %in% available) {
  cat("  Exact date not found, searching available snapshots...\n")
  matches <- grep("2000-12-27", available, value = TRUE)
  if (length(matches) > 0) {
    snapshot_date <- matches[1]
    cat("  Using:", snapshot_date, "\n")
  } else {
    stop("No snapshot found for 2000-12-27")
  }
}

cat("  Snapshot date:", snapshot_date, "\n")

# Run reforecast
refc_result <- run_reforecast_from_snapshot(
  snapshot_file = snap_file,
  snapshot_date = snapshot_date,
  forcing       = forcing_fcst,
  uptribs       = NULL,
  pars          = pars,
  basin_name    = basin,
  dt_hours      = dt_hours,
  ae_tbl        = ae_tbl
)

cat("  Reforecast complete:", length(refc_result$flow_cfs), "timesteps\n")
cat("  Reforecast flow range:", round(min(refc_result$flow_cfs), 2), "to",
    round(max(refc_result$flow_cfs), 2), "cfs\n\n")

# Aggregate reforecast to daily
refc_daily <- data.frame(
  date     = as.Date(refc_dt),
  flow_cfs = refc_result$flow_cfs
)
refc_daily <- aggregate(flow_cfs ~ date, data = refc_daily, FUN = mean)
names(refc_daily)[2] <- "refc_cfs"

# ============================================================================
# STEP 3: COMPARE REFORECAST VS POR (sub-daily)
# ============================================================================

cat("=== Step 3: Comparing reforecast vs continuous POR ===\n")

overlap_start <- max(min(refc_dt), min(forcing_dt))
overlap_end   <- min(max(refc_dt), max(forcing_dt))

por_idx  <- which(forcing_dt >= overlap_start & forcing_dt <= overlap_end)
refc_idx <- which(refc_dt >= overlap_start & refc_dt <= overlap_end)

if (length(por_idx) > 0 && length(refc_idx) > 0) {
  min_len <- min(length(por_idx), length(refc_idx))
  por_flow  <- result_por$flow_cfs[por_idx[1:min_len]]
  refc_flow <- refc_result$flow_cfs[refc_idx[1:min_len]]
  overlap_dt <- forcing_dt[por_idx[1:min_len]]

  flow_diff <- refc_flow - por_flow
  cat("  Overlap period:", format(overlap_start), "to", format(overlap_end), "\n")
  cat("  Overlap timesteps:", min_len, "\n")
  cat("  Max |difference|:", round(max(abs(flow_diff)), 4), "cfs\n")
  cat("  RMSE:", round(sqrt(mean(flow_diff^2)), 4), "cfs\n\n")
} else {
  cat("  No overlapping period found between POR and reforecast.\n")
  overlap_dt <- NULL
}

# ============================================================================
# STEP 4: PLOTS
# ============================================================================

cat("=== Step 4: Generating figures ===\n")

# --- Figure 1: Full POR — simulated vs observed (daily) ---
df_por_plot <- merge(por_daily, obs_daily, by = "date", all.x = TRUE)

p1 <- ggplot(df_por_plot, aes(x = date)) +
  geom_line(aes(y = flow_cfs, color = "Observed"), linewidth = 0.3) +
  geom_line(aes(y = sim_cfs, color = "Simulated"), linewidth = 0.3) +
  scale_color_manual(values = c("Observed" = "black", "Simulated" = "steelblue")) +
  labs(title = paste(basin, "— POR: Simulated vs Observed (daily mean)"),
       subtitle = sprintf("NSE = %.3f, Bias = %.2f cfs", nse_por, mean(por_vs_obs$diff)),
       x = NULL, y = "Flow (cfs)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "01_por_sim_vs_obs.png"), p1, width = 14, height = 5, dpi = 150)

# --- Figure 2: Zoomed POR around reforecast period ---
zoom_start <- as.Date("2000-10-01")
zoom_end   <- as.Date("2001-03-31")
df_zoom <- df_por_plot[df_por_plot$date >= zoom_start & df_por_plot$date <= zoom_end, ]

if (nrow(df_zoom) > 0) {
  p2 <- ggplot(df_zoom, aes(x = date)) +
    geom_line(aes(y = flow_cfs, color = "Observed"), linewidth = 0.5) +
    geom_line(aes(y = sim_cfs, color = "Simulated"), linewidth = 0.5) +
    geom_vline(xintercept = as.Date("2000-12-27"), linetype = "dashed", color = "red", linewidth = 0.4) +
    annotate("text", x = as.Date("2000-12-27"),
             y = max(df_zoom$sim_cfs, df_zoom$flow_cfs, na.rm = TRUE) * 0.95,
             label = "Snapshot", hjust = -0.1, color = "red", size = 3) +
    scale_color_manual(values = c("Observed" = "black", "Simulated" = "steelblue")) +
    labs(title = paste(basin, "— WY2001 Context"),
         x = NULL, y = "Flow (cfs)", color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "02_por_zoom_wy2001.png"), p2, width = 10, height = 5, dpi = 150)
}

# --- Figure 3: Reforecast vs POR vs Observed (daily) ---
df_refc_comp <- merge(refc_daily, obs_daily, by = "date", all.x = TRUE)
df_refc_comp <- merge(df_refc_comp, por_daily, by = "date", all.x = TRUE)

p3 <- ggplot(df_refc_comp, aes(x = date)) +
  geom_line(aes(y = flow_cfs, color = "Observed"), linewidth = 0.6) +
  geom_line(aes(y = sim_cfs, color = "POR Simulation"), linewidth = 0.5) +
  geom_line(aes(y = refc_cfs, color = "Reforecast (MM)"), linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("Observed" = "black",
                                "POR Simulation" = "steelblue",
                                "Reforecast (MM)" = "firebrick")) +
  labs(title = paste(basin, "— Reforecast vs POR vs Observed"),
       subtitle = paste("Reforecast from", snapshot_date, "using MM forcing"),
       x = NULL, y = "Daily Mean Flow (cfs)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "03_reforecast_vs_por_vs_obs.png"), p3, width = 10, height = 5, dpi = 150)

# --- Figure 4: Reforecast vs Observed only ---
refc_vs_obs <- df_refc_comp[!is.na(df_refc_comp$flow_cfs) & !is.na(df_refc_comp$refc_cfs), ]
if (nrow(refc_vs_obs) > 0) {
  nse_refc <- 1 - sum((refc_vs_obs$refc_cfs - refc_vs_obs$flow_cfs)^2) /
                  sum((refc_vs_obs$flow_cfs - mean(refc_vs_obs$flow_cfs))^2)
  bias_refc <- mean(refc_vs_obs$refc_cfs - refc_vs_obs$flow_cfs)

  p4 <- ggplot(refc_vs_obs, aes(x = date)) +
    geom_line(aes(y = flow_cfs, color = "Observed"), linewidth = 0.6) +
    geom_line(aes(y = refc_cfs, color = "Reforecast (MM)"), linewidth = 0.6) +
    scale_color_manual(values = c("Observed" = "black", "Reforecast (MM)" = "firebrick")) +
    labs(title = paste(basin, "— Reforecast vs Observed"),
         subtitle = sprintf("NSE = %.3f, Bias = %.2f cfs", nse_refc, bias_refc),
         x = NULL, y = "Daily Mean Flow (cfs)", color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "04_reforecast_vs_obs.png"), p4, width = 10, height = 5, dpi = 150)
} else {
  cat("  No overlapping obs data for reforecast period.\n")
  nse_refc <- NA
  bias_refc <- NA
}

cat("  Figures saved to:", fig_dir, "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY — ", basin, "\n")
cat(strrep("=", 70), "\n")
cat("POR simulation:\n")
cat("  Timesteps:     ", length(result_por$flow_cfs), "\n")
cat("  Period:        ", format(min(forcing_dt)), "to", format(max(forcing_dt)), "\n")
cat("  Snapshots:     ", result_por$snapshot_info$n_snapshots, "(daily at 12Z)\n")
cat("  vs Observed:    NSE =", round(nse_por, 4), ", Bias =", round(mean(por_vs_obs$diff), 2), "cfs\n")
cat("Reforecast:\n")
cat("  From snapshot: ", snapshot_date, "\n")
cat("  Timesteps:     ", length(refc_result$flow_cfs), "\n")
cat("  Period:        ", format(min(refc_dt)), "to", format(max(refc_dt)), "\n")
if (!is.na(nse_refc)) {
  cat("  vs Observed:    NSE =", round(nse_refc, 4), ", Bias =", round(bias_refc, 2), "cfs\n")
}
cat("Outputs:\n")
cat("  Snapshots:     ", snap_file, "\n")
cat("  Figures:       ", fig_dir, "\n")
cat(strrep("=", 70), "\n")
