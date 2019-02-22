rm(list=ls());gc()

METADATA <- "metadata.csv"
CACHE <- "./cache/"
RESULT_DIR <- "/data/ethoscope_results"
REMOTE_DIR <- "ftp://nas.lab.gilest.ro/auto_generated_data/ethoscope_results/"
SD_START_END <- hours(c(20,24))
REF_SLEEP_WINDOW <- hours(c(0, 3))
TARGET_SLEEP_WINDOW <- REF_SLEEP_WINDOW + days(1)
PDF_NAME <- "time_window_dsd.pdf"
source("../trunk.R")
make_plot()

