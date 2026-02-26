#!/usr/bin/env Rscript
# test_prepare_forcing.R
#
# Verifies that minimal forcing (year, month, day, hour, map_mm, mat_degc)
# produces identical results to pre-computed forcing when both use the same
# pet/etd/ptps derivation method (prepare_forcing).
#
# Workflow:
#   1. Load original forcingSAKW1, strip ptps/pet_mm/etd_mm to create minimal forcing
#   2. Run prepare_forcing() offline to compute those columns -> "reference" forcing
#   3. Write CSVs for both versions
#   4. Run sac_snow_uh_lagk_states_with_snapshots with reference forcing (pre-computed)
#   5. Run sac_snow_uh_lagk_states_with_snapshots with minimal forcing (auto-computed)
#   6. Run reforecast from snapshot for both
#   7. Compare and plot

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(ggplot2)
  load_all(".")
})

# ============================================================================
# CONFIGURATION
# ============================================================================

basin       <- "SAKW1"
dt_hours    <- 6
out_dir     <- "test_prepare_forcing_output"
snap_dir    <- file.path(out_dir, "snapshots")
csv_dir     <- file.path(out_dir, "csv")
fig_dir     <- file.path(out_dir, "figures")
restart_dir <- file.path(out_dir, "restarts")

for (d in c(out_dir, snap_dir, csv_dir, fig_dir, restart_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Load package data
data(forcingSAKW1)
data(parsSAKW1)
data(upflowSAKW1)

forcing_orig <- forcingSAKW1
pars         <- parsSAKW1
uptribs      <- upflowSAKW1
zone_names   <- names(forcing_orig)

# Add pxtemp and talr parameters (needed for rsnwelev ptps computation)
# These would normally be in the parameter file; added here for testing.
extra_pars <- data.frame(
  p_name = c("pxtemp_SAKW1-1", "pxtemp_SAKW1-2", "talr_SAKW1-1", "talr_SAKW1-2"),
  name   = c("pxtemp",         "pxtemp",          "talr",         "talr"),
  type   = c("snow",           "snow",             "snow",         "snow"),
  zone   = c("SAKW1-1",        "SAKW1-2",          "SAKW1-1",      "SAKW1-2"),
  value  = c(-0.867040346419604, 2.9642471707163, 0.711034803854563, 0.745278555068483),
  stringsAsFactors = FALSE
)
pars <- rbind(as.data.frame(pars), extra_pars)

# ============================================================================
# NOTE: ae_tbl for ptps — required for prepare_forcing
# ============================================================================
# Create synthetic area-elevation curve for SAKW1 testing.
# SAKW1-1 elev = 762m, SAKW1-2 elev = 1599m (from pars)
# Curve spans realistic elevation ranges around those values.

# ============================================================================
# NOTE: ae_tbl for ptps — required for prepare_forcing
# ============================================================================
# Real area-elevation curve for SAKW1

ae_tbl <- data.frame(
  quantile = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
               0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
  `SAKW1-1` = c(305, 366, 427, 488, 518, 549, 580, 610, 641, 671,
                701, 732, 762, 793, 823, 854, 914, 975, 1067, 1189, 1372),
  `SAKW1-2` = c(610, 762, 914, 1006, 1067, 1128, 1189, 1250, 1311, 1372,
                1433, 1494, 1555, 1616, 1676, 1737, 1829, 1920, 2042, 2195, 2438),
  check.names = FALSE
)
cat("  ae_tbl loaded (SAKW1):\n")
print(head(ae_tbl))

# ============================================================================
# STEP 1: CREATE MINIMAL FORCING (strip computed columns)
# ============================================================================

cat("=== Step 1: Creating minimal forcing ===\n")

forcing_minimal <- lapply(forcing_orig, function(df) {
  keep_cols <- c("year", "month", "day", "hour", "map_mm", "mat_degc")
  df[, keep_cols, drop = FALSE]
})

cat("  Minimal forcing columns:", paste(names(forcing_minimal[[1]]), collapse = ", "), "\n\n")

# ============================================================================
# STEP 2: RUN prepare_forcing() OFFLINE TO BUILD REFERENCE FORCING
# ============================================================================

cat("=== Step 2: Running prepare_forcing() offline ===\n")

forcing_reference <- prepare_forcing(
  forcing  = lapply(forcing_minimal, as.data.frame),  # copy so minimal stays clean
  pars     = pars,
  dt_hours = dt_hours,
  ae_tbl   = ae_tbl
)

cat("  Reference forcing columns:", paste(names(forcing_reference[[1]]), collapse = ", "), "\n\n")

# ============================================================================
# STEP 3: WRITE CSVs
# ============================================================================

cat("=== Step 3: Writing CSV files ===\n")

for (i in seq_along(zone_names)) {
  z <- zone_names[i]

  # Minimal (what you'd distribute / store)
  write.csv(
    forcing_minimal[[i]],
    file.path(csv_dir, sprintf("forcing_minimal_%s.csv", z)),
    row.names = FALSE
  )

  # Reference (pre-computed pet/etd, same method as prepare_forcing)
  write.csv(
    forcing_reference[[i]],
    file.path(csv_dir, sprintf("forcing_reference_%s.csv", z)),
    row.names = FALSE
  )

  # Original (for comparison with the old pet/etd values)
  write.csv(
    forcing_orig[[i]],
    file.path(csv_dir, sprintf("forcing_original_%s.csv", z)),
    row.names = FALSE
  )
}

cat("  CSVs written to:", csv_dir, "\n\n")

# ============================================================================
# STEP 4: RUN WITH REFERENCE FORCING (pre-computed pet/etd)
# ============================================================================

cat("=== Step 4: Simulation with REFERENCE forcing (pre-computed) ===\n")

snap_file_ref <- file.path(snap_dir, "snapshots_reference.rds")

result_ref <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours        = dt_hours,
  forcing         = forcing_reference,
  uptribs         = uptribs,
  pars            = pars,
  ae_tbl          = ae_tbl,
  save_snapshots  = TRUE,
  snapshot_file   = snap_file_ref,
  snapshot_interval = "daily",
  snapshot_hour   = 0
)

cat("  Reference run complete. Snapshots:", result_ref$snapshot_info$n_snapshots, "\n\n")

# ============================================================================
# STEP 5: RUN WITH MINIMAL FORCING (auto-computed inside pipeline)
# ============================================================================

cat("=== Step 5: Simulation with MINIMAL forcing (auto-computed) ===\n")

snap_file_min <- file.path(snap_dir, "snapshots_minimal.rds")

# Convert to plain data.frame so prepare_forcing can add columns
forcing_minimal_df <- lapply(forcing_minimal, as.data.frame)

result_min <- sac_snow_uh_lagk_states_with_snapshots(
  dt_hours        = dt_hours,
  forcing         = forcing_minimal_df,
  uptribs         = uptribs,
  pars            = pars,
  ae_tbl          = ae_tbl,
  save_snapshots  = TRUE,
  snapshot_file   = snap_file_min,
  snapshot_interval = "daily",
  snapshot_hour   = 0
)

cat("  Minimal run complete. Snapshots:", result_min$snapshot_info$n_snapshots, "\n\n")

# ============================================================================
# STEP 6: COMPARE FULL SIMULATION FLOWS
# ============================================================================

cat("=== Step 6: Comparing full simulation flows ===\n")

flow_diff <- result_ref$flow_cfs - result_min$flow_cfs
cat("  Max absolute difference:", max(abs(flow_diff)), "\n")
cat("  Mean absolute difference:", mean(abs(flow_diff)), "\n")
cat("  Flows identical:", all.equal(result_ref$flow_cfs, result_min$flow_cfs), "\n\n")

# Build datetime vector for plotting
forcing_dt <- as.POSIXct(
  ISOdatetime(
    forcing_orig[[1]]$year,
    forcing_orig[[1]]$month,
    forcing_orig[[1]]$day,
    forcing_orig[[1]]$hour,
    0, 0, tz = "UTC"
  ),
  tz = "UTC"
)

# ============================================================================
# STEP 7: REFORECAST FROM SNAPSHOT — BOTH VERSIONS
# ============================================================================

cat("=== Step 7: Running reforecasts from snapshot ===\n")

# Pick a reforecast date near the middle
available_dates <- list_snapshot_dates(snap_file_ref)
reforecast_date <- available_dates[length(available_dates) %/% 2]
reforecast_dt   <- as.POSIXct(reforecast_date, tz = "UTC")
reforecast_idx  <- which(forcing_dt == reforecast_dt)[1]

cat("  Reforecast date:", reforecast_date, "(index", reforecast_idx, ")\n")

# Subset forcing and uptribs from reforecast date onward
forcing_ref_sub <- lapply(forcing_reference, function(df) df[reforecast_idx:nrow(df), ])
forcing_min_sub <- lapply(forcing_minimal_df, function(df) df[reforecast_idx:nrow(df), ])
uptribs_sub     <- lapply(uptribs, function(df) df[reforecast_idx:nrow(df), ])

# Reforecast with reference forcing
cat("  Running reforecast (reference)...\n")
refc_ref <- run_reforecast_from_snapshot(
  snapshot_file = snap_file_ref,
  snapshot_date = reforecast_date,
  forcing       = forcing_ref_sub,
  uptribs       = uptribs_sub,
  pars          = pars,
  basin_name    = basin,
  dt_hours      = dt_hours,
  temp_dir      = restart_dir,
  ae_tbl        = ae_tbl
)

# Reforecast with minimal forcing
cat("  Running reforecast (minimal)...\n")
refc_min <- run_reforecast_from_snapshot(
  snapshot_file = snap_file_min,
  snapshot_date = reforecast_date,
  forcing       = forcing_min_sub,
  uptribs       = uptribs_sub,
  pars          = pars,
  basin_name    = basin,
  dt_hours      = dt_hours,
  temp_dir      = restart_dir,
  ae_tbl        = ae_tbl
)

refc_diff <- refc_ref$flow_cfs - refc_min$flow_cfs
cat("  Reforecast max absolute difference:", max(abs(refc_diff)), "\n")
cat("  Reforecast flows identical:", all.equal(refc_ref$flow_cfs, refc_min$flow_cfs), "\n\n")

# ============================================================================
# STEP 8: PLOTS
# ============================================================================

cat("=== Step 8: Generating figures ===\n")

# --- Figure 1: Full simulation comparison ---
df_full <- data.frame(
  datetime  = forcing_dt,
  reference = result_ref$flow_cfs,
  minimal   = result_min$flow_cfs
)

p1 <- ggplot(df_full, aes(x = datetime)) +
  geom_line(aes(y = reference, color = "Reference (pre-computed)"), linewidth = 0.4) +
  geom_line(aes(y = minimal,   color = "Minimal (auto-computed)"), linewidth = 0.4, linetype = "dashed") +
  scale_color_manual(values = c("Reference (pre-computed)" = "steelblue",
                                "Minimal (auto-computed)"   = "firebrick")) +
  labs(title = "Full Simulation: Reference vs Minimal Forcing",
       subtitle = sprintf("Max |diff| = %.6f cfs", max(abs(flow_diff))),
       x = NULL, y = "Flow (cfs)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "01_full_simulation_comparison.png"), p1,
       width = 12, height = 5, dpi = 150)

# --- Figure 2: Difference time series ---
df_diff <- data.frame(datetime = forcing_dt, diff_cfs = flow_diff)

p2 <- ggplot(df_diff, aes(x = datetime, y = diff_cfs)) +
  geom_line(linewidth = 0.3, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Flow Difference (Reference - Minimal)",
       x = NULL, y = "Difference (cfs)") +
  theme_minimal()

ggsave(file.path(fig_dir, "02_flow_difference.png"), p2,
       width = 12, height = 4, dpi = 150)

# --- Figure 3: Reforecast comparison ---
refc_dt <- forcing_dt[reforecast_idx:length(forcing_dt)]

df_refc <- data.frame(
  datetime  = refc_dt,
  reference = refc_ref$flow_cfs,
  minimal   = refc_min$flow_cfs
)

# Also overlay the original continuous run for context
df_refc$continuous <- result_ref$flow_cfs[reforecast_idx:length(result_ref$flow_cfs)]

p3 <- ggplot(df_refc, aes(x = datetime)) +
  geom_line(aes(y = continuous, color = "Continuous run"), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = reference,  color = "Reforecast (reference)"), linewidth = 0.4) +
  geom_line(aes(y = minimal,    color = "Reforecast (minimal)"), linewidth = 0.4, linetype = "dashed") +
  scale_color_manual(values = c("Continuous run"            = "grey50",
                                "Reforecast (reference)"    = "steelblue",
                                "Reforecast (minimal)"      = "firebrick")) +
  labs(title = sprintf("Reforecast from %s", reforecast_date),
       subtitle = sprintf("Max |diff| = %.6f cfs", max(abs(refc_diff))),
       x = NULL, y = "Flow (cfs)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "03_reforecast_comparison.png"), p3,
       width = 12, height = 5, dpi = 150)

# --- Figure 4: Computed vs original pet_mm / etd_mm (sanity check) ---
# Compare what prepare_forcing computed vs what was in the original .rda
df_pet <- data.frame(
  datetime     = forcing_dt,
  pet_computed = forcing_reference[[1]]$pet_mm,
  pet_original = forcing_orig[[1]]$pet_mm,
  etd_computed = forcing_reference[[1]]$etd_mm,
  etd_original = forcing_orig[[1]]$etd_mm
)

p4a <- ggplot(df_pet[1:min(1460, nrow(df_pet)), ], aes(x = datetime)) +
  geom_line(aes(y = pet_original, color = "Original pet_mm"), linewidth = 0.4) +
  geom_line(aes(y = pet_computed,  color = "Computed pet_mm"), linewidth = 0.4, linetype = "dashed") +
  labs(title = "PET Comparison: Original vs Hargreaves-Samani (Zone 1, first year)",
       x = NULL, y = "pet_mm", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

p4b <- ggplot(df_pet[1:min(1460, nrow(df_pet)), ], aes(x = datetime)) +
  geom_line(aes(y = etd_original, color = "Original etd_mm"), linewidth = 0.4) +
  geom_line(aes(y = etd_computed,  color = "Computed etd_mm"), linewidth = 0.4, linetype = "dashed") +
  labs(title = "ETD Comparison: Original vs Computed (Zone 1, first year)",
       x = NULL, y = "etd_mm", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "04a_pet_comparison.png"), p4a, width = 12, height = 4, dpi = 150)
ggsave(file.path(fig_dir, "04b_etd_comparison.png"), p4b, width = 12, height = 4, dpi = 150)

# --- Figure 4c: ptps comparison (rsnwelev-computed vs original) ---
df_ptps <- data.frame(
  datetime     = forcing_dt,
  ptps_computed = forcing_reference[[1]]$ptps,
  ptps_original = forcing_orig[[1]]$ptps,
  mat_degc      = forcing_orig[[1]]$mat_degc
)

p4c <- ggplot(df_ptps[1:min(1460, nrow(df_ptps)), ], aes(x = datetime)) +
  geom_line(aes(y = ptps_original, color = "Original ptps"), linewidth = 0.4) +
  geom_line(aes(y = ptps_computed,  color = "Computed ptps (rsnwelev)"), linewidth = 0.4, linetype = "dashed") +
  labs(title = "PTPS Comparison: Original vs rsnwelev-computed (Zone 1, first year)",
       subtitle = "Note: synthetic ae_tbl — values will differ from original but ref vs minimal should match",
       x = NULL, y = "ptps (fraction snow)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "04c_ptps_comparison.png"), p4c, width = 12, height = 4, dpi = 150)

# Quick sanity check: ptps at warm vs cold temps
cat("\n  PTPS sanity check (computed, zone 1):\n")
warm_idx <- which(forcing_reference[[1]]$mat_degc > 10)[1:3]
cold_idx <- which(forcing_reference[[1]]$mat_degc < -5)[1:3]
cat("    Warm temps (>10C): mat =", forcing_reference[[1]]$mat_degc[warm_idx],
    "-> ptps =", forcing_reference[[1]]$ptps[warm_idx], "\n")
cat("    Cold temps (<-5C): mat =", forcing_reference[[1]]$mat_degc[cold_idx],
    "-> ptps =", forcing_reference[[1]]$ptps[cold_idx], "\n")

cat("\n  Figures saved to:", fig_dir, "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")
cat("Full simulation:\n")
cat("  Max |flow diff|:  ", max(abs(flow_diff)), "cfs\n")
cat("  Flows identical:  ", isTRUE(all.equal(result_ref$flow_cfs, result_min$flow_cfs)), "\n")
cat("Reforecast:\n")
cat("  Max |flow diff|:  ", max(abs(refc_diff)), "cfs\n")
cat("  Flows identical:  ", isTRUE(all.equal(refc_ref$flow_cfs, refc_min$flow_cfs)), "\n")
cat("Outputs:\n")
cat("  CSVs:     ", csv_dir, "\n")
cat("  Snapshots:", snap_dir, "\n")
cat("  Figures:  ", fig_dir, "\n")
cat(strrep("=", 70), "\n")
