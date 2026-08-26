
# Load data
# Uses `input_file` and `log_dir` variable set by run_msstats.Rmd from params

msstats_format <- FragPipetoMSstatsFormat(
  input = input_file,
  useUniquePeptide = TRUE,
  removeFewMeasurements = TRUE,
  removeProtein_with1Feature = FALSE,
  summaryforMultipleRows = max,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = file.path(log_dir, "MSstats_import.log")
)