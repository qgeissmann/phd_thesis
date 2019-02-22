rm(list=ls());gc()

METADATA <- "metadata.csv"
CACHE <- "./cache/"
RESULT_DIR <- "/data/ethoscope_results"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"
SD_START_END <- hours(c(0,12))
REF_SLEEP_WINDOW <- hours(c(0, 3)) - hours(12)
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(1)
PDF_NAME <- "L_phase_dsd.pdf"

source("../trunk.R")

stat_rebound_dt <- make_plot()
stat_rebound_dt <- stat_rebound_dt[dt[, .(interactions = sum(interactions)), by=id]]
stat_rebound_dt[sdi==10, round(mean_cl_boot(180 * quiet_rebound_day3h_diff),3), by="sex"]
stat_rebound_dt[sdi==10, round(mean_cl_boot(interactions),2), by="sex"]
stat_rebound_dt[, .N, by="sex,treatment"]


