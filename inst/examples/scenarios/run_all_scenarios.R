# =============================================================================
# Run every PredictProR scenario script and print one PASS / SKIP / FAIL table
# =============================================================================
#   source(system.file("examples", "scenarios", "run_all_scenarios.R", package = "PredictProR"))
#
# Each script is plain user code (full model_execute() calls); every call is
# recorded with log_result(). SKIP = optional software not available
# (ASReml-R, a Python backend, Java for Beagle). A full run with every backend
# takes roughly 20-30 minutes on a desktop.
# To run only some scripts:  only <- c("01", "03")  before sourcing this file.
# =============================================================================
scenario_dir <- system.file("examples", "scenarios", package = "PredictProR")
if (!nzchar(scenario_dir)) scenario_dir <- "."

scripts <- sort(list.files(scenario_dir, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
scripts <- scripts[!grepl("^00_", basename(scripts))]
if (exists("only", inherits = FALSE)) scripts <- scripts[substr(basename(scripts), 1, 2) %in% only]

# a fresh log for this run
assign("scenario_log", data.frame(scenario = character(), status = character(), seconds = numeric(),
                                  note = character(), stringsAsFactors = FALSE), envir = globalenv())
script_of_row <- character()

for (s in scripts) {
  cat("\n\n#################################################################\n")
  cat("#", basename(s), "\n")
  cat("#################################################################\n")
  before <- nrow(get("scenario_log", envir = globalenv()))
  ok <- tryCatch({ source(s, local = new.env(), chdir = TRUE); TRUE }, error = function(e) {
    tab <- get("scenario_log", envir = globalenv())
    tab[nrow(tab) + 1L, ] <- list(basename(s), "FAIL", NA_real_, paste("script stopped:", conditionMessage(e)))
    assign("scenario_log", tab, envir = globalenv())
    FALSE
  })
  after <- nrow(get("scenario_log", envir = globalenv()))
  script_of_row <- c(script_of_row, rep(basename(s), after - before))
}

summary_table <- cbind(script = script_of_row, get("scenario_log", envir = globalenv()))
rownames(summary_table) <- NULL
old <- options(width = 250)
cat("\n\n======================= ALL SCENARIOS =======================\n")
print(summary_table, row.names = FALSE, right = FALSE)
cat(sprintf("\nTOTAL  PASS %d   SKIP %d   FAIL %d\n", sum(summary_table$status == "PASS"),
            sum(summary_table$status == "SKIP"), sum(summary_table$status == "FAIL")))
options(old)
invisible(summary_table)
