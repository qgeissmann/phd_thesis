rm(list=ls());gc()

METADATA <- "metadata.csv"
CACHE <- "./cache/"
RESULT_DIR <- "/data/ethoscope_results"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"
SD_START_END <- hours(c(0,8))
REF_SLEEP_WINDOW <- hours(c(8, 11)) - hours(24)
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(1)
PDF_NAME <- "zt0-8_phase_dsd.pdf"


source("../trunk.R")
make_plot()


