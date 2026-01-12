#!/usr/bin/env Rscript

# test_restart_general.R - Verify restart capability for multiple basins

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  load_all(".")
})

# -----------------------------
# User settings
# -----------------------------
basins <- c("SAKW1", "WCHW1", "WGCM8")
dt_hours <- 6
split_dt <- as.POSIXct("2001-03-31 18:00", tz = "UTC")  # Adjust as needed
days_window <- 7

# Create output directories
if(!dir.exists("states")) dir.create("states")
if(!dir.exists("figures")) dir.create("figures")

# -----------------------------
# Basin loop
# -----------------------------
for(basin in basins) {
  cat("\n==============================\n")
  cat("Processing basin:", basin, "\n")
  cat("==============================\n")
  
  # Load forcing and parameters
  forcing_file <- paste0("../data/forcing", basin, ".rda")
  pars_file <- paste0("../data/pars", basin, ".rda")
  
  load(forcing_file)
  load(pars_file)
  
  forcing_list <- get(paste0("forcing", basin))
  pars <- get(paste0("pars", basin))
    
  # -----------------------------
  # Initialize upflow list
  # -----------------------------
  upflow_var <- paste0("upflow", basin)
  has_upflow <- exists(upflow_var) && !is.null(get(upflow_var))
  upflow_list <- if(has_upflow) get(upflow_var) else NULL
  
  # Forcing datetimes
  forcing_dt <- as.POSIXct(ISOdatetime(
    forcing_list[[1]]$year,
    forcing_list[[1]]$month,
    forcing_list[[1]]$day,
    forcing_list[[1]]$hour,
    0, 0,
    tz = "UTC"
  ))
  
  cat("Total forcing period:", min(forcing_dt), "to", max(forcing_dt), "\n")
  cat("Total timesteps:", length(forcing_dt), "\n")
  
  split_idx <- which(forcing_dt == split_dt)[1]
  cat("Split at index:", split_idx, "datetime:", forcing_dt[split_idx], "\n")
  
  # ------------------------------
  # RUN 1: Start to split
  # ------------------------------
  cat("\n===== RUN 1: Start to split date =====\n")
  
  forcing1 <- lapply(forcing_list, function(df) df[1:split_idx, ])
  upflow1 <- if(has_upflow) lapply(upflow_list, function(df) df[1:split_idx, ]) else NULL

                                   
  out1 <- sac_snow_uh_lagk_states(
      dt_hours,
      forcing1,
      upflow1,                      # NULL for non-upflow basins
      pars
  )
  
  # Save restart files
  restart_csv <- paste0("states/restart_", basin, ".csv")
    restart1 <- out1$restart
    n_zones <- length(restart1$uztwc)

    # Build restart dataframe
    restart_df <- data.frame(
      zone = 1:n_zones,
      uztwc = restart1$uztwc,
      uzfwc = restart1$uzfwc,
      lztwc = restart1$lztwc,
      lzfsc = restart1$lzfsc,
      lzfpc = restart1$lzfpc,
      adimc = restart1$adimc,
      taprev = restart1$taprev
    )

    # Add snow states
    for (i in 1:19) {
      restart_df[[paste0("cs_", i)]] <- restart1$cs[i, ]
    }

  # SAC states
  restart_csv <- paste0("states/restart_", basin, ".csv")
  # Build restart_df manually here
  write.csv(restart_df, restart_csv, row.names = FALSE)

  # UH restart (always)
  uh_rds <- paste0("states/uh_restart_", basin, ".rds")
  saveRDS(out1$uh_restart, uh_rds)
  cat("Saved", uh_rds, "\n")

  # Lag-K restart only if upflow exists
  if(has_upflow && !is.null(out1$lagk_restart)) {
      lagk_rds <- paste0("states/lagk_restart_", basin, ".rds")
      saveRDS(out1$lagk_restart, lagk_rds)
      cat("Saved", lagk_rds, "\n")
  } else {
      lagk_rds <- NULL
  }
                                   
  # ------------------------------
  # RUN 2: Split to end (warm start)
  # ------------------------------
  cat("\n===== RUN 2: Split date to end (WARM START) =====\n")
  
  forcing2 <- lapply(forcing_list, function(df) df[(split_idx+1):nrow(df), ])
  upflow2 <- if(has_upflow) lapply(upflow_list, function(df) df[(split_idx+1):nrow(df), ]) else NULL
  
  restart_file <- restart_csv
  uh_restart_file <- if(has_upflow) uh_rds else NULL
  lagk_restart_file <- if(has_upflow) lagk_rds else NULL
  
  out2 <- sac_snow_uh_lagk_states(
    dt_hours,
    forcing2,
    upflow2,                      # NULL if no upflow
    pars,
    restart_file = restart_csv,
    lagk_restart_file = lagk_rds, # NULL if no upflow
    uh_restart_file = uh_rds      # always load
  )
  
  # ------------------------------
  # RUN 3: Full period (cold start)
  # ------------------------------
  cat("\n===== RUN 3: Full period (COLD START) =====\n")
  
  out3 <- sac_snow_uh_lagk_states(
    dt_hours,
    forcing_list,
    upflow_list,
    pars
  )
  
  # ------------------------------
  # TCI comparison
  # ------------------------------
  tci_file <- paste0("../data/tci", basin, ".rda")
  load(tci_file)
  tci_ref <- get(paste0("tci", basin))
  
  cat("\nTCI comparison for basin", basin, "\n")
  cat("Run1 final TCI (zone1):", tail(out1$states$tci_1, 3), "\n")
  cat("Run2 initial TCI (zone1):", head(out2$states$tci_1, 3), "\n")
  cat("Run3 TCI at split:", out3$states$tci_1[split_idx:(split_idx+2)], "\n")
  
  # ------------------------------
  # Flow comparison
  # ------------------------------
  flow_combined <- c(out1$flow_cfs, out2$flow_cfs)
  flow_full <- out3$flow_cfs
  
  restart_idx <- length(out1$flow_cfs)
  
  plot_df <- data.frame(
    datetime = forcing_dt,
    full_run = flow_full,
    combined = flow_combined,
    difference = flow_combined - flow_full
  )
  
  restart_time <- plot_df$datetime[restart_idx]
  
  plot_zoom <- plot_df %>%
    filter(datetime >= restart_time - days_window*86400 &
             datetime <= restart_time + days_window*86400)
  
  p_flow <- ggplot(plot_zoom) +
    geom_line(aes(datetime, full_run, color="Full run"), size=0.8) +
    geom_line(aes(datetime, combined, color="Restart run"), size=0.8, linetype="dashed") +
    geom_vline(xintercept=restart_time, linetype="dotted") +
    scale_color_manual(values=c("Full run"="blue", "Restart run"="red")) +
    labs(title=paste("Flow comparison around restart:", basin),
         x="Date", y="Flow (cfs)", color="") +
    theme_minimal()
  
  p_diff <- ggplot(plot_zoom) +
    geom_line(aes(datetime, difference), size=0.8) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=restart_time, linetype="dotted") +
    labs(title=paste("Flow difference (Restart − Full):", basin),
         x="Date", y="Difference (cfs)") +
    theme_minimal()
  
  plot_file <- paste0("figures/restart_zoom_week_", basin, ".png")
  ggsave(plot_file, grid.arrange(p_flow, p_diff, ncol=1), width=10, height=8)
  cat("Saved", plot_file, "\n\n")
}