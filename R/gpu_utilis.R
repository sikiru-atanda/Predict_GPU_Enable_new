gp_is_expected_gpu_probe_error <- function(e) {
  msg <- tolower(conditionMessage(e) %||% "")
  any(grepl(
    paste(
      c(
        "no module named 'torch'",
        'no module named "torch"',
        "modulenotfounderror",
        "reticulate",
        "python shared library",
        "unable to locate python",
        "python specified in reticulate_python",
        "installation of python not found"
      ),
      collapse = "|"
    ),
    msg,
    perl = TRUE
  ))
}

detect_physical_num_gpus <- function() {
  nvidia_smi <- Sys.which("nvidia-smi")
  if (!nzchar(nvidia_smi)) {
    return(NA_integer_)
  }
  out <- tryCatch(
    system2(
      nvidia_smi,
      args = c("--query-gpu=index", "--format=csv,noheader,nounits"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  out <- trimws(out)
  out <- out[nzchar(out)]
  if (!length(out)) {
    return(NA_integer_)
  }
  as.integer(length(out))
}

# Detect number of GPUs
detect_num_gpus <- function(purpose = c("gp", "dl", "general")) {
  purpose <- match.arg(purpose)
  physical <- detect_physical_num_gpus()
  if (is.finite(physical) && physical > 0L) {
    return(as.integer(physical))
  }
  tryCatch({
    status <- gp_python_runtime_status(initialize = TRUE, purpose = purpose)
    as.integer(status$num_gpus %||% 0L)
  }, error = function(e) {
    if (!gp_is_expected_gpu_probe_error(e)) {
      logger::log_warn("GPU detection failed; assuming no GPUs: {conditionMessage(e)}")
    }
    0L
  })
}

get_gpu_usage <- function() {
  nvidia_smi <- Sys.which("nvidia-smi")
  if (!nzchar(nvidia_smi)) {
    return(NA_real_)
  }
  out <- tryCatch(
    system2(
      nvidia_smi,
      args = c("--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )
  values <- suppressWarnings(as.numeric(trimws(out)))
  if (!length(values) || !any(is.finite(values))) {
    return(NA_real_)
  }
  mean(values, na.rm = TRUE)
}

gp_gpu_busy_threshold <- function() {
  val <- suppressWarnings(as.numeric(Sys.getenv("PREDICTPRO_GPU_BUSY_THRESHOLD", "85")))
  if (!is.finite(val)) 85 else max(0, min(100, val))
}

gp_gpu_is_busy <- function(usage = get_gpu_usage(),
                           threshold = gp_gpu_busy_threshold()) {
  is.finite(usage) && is.finite(threshold) && as.numeric(usage) >= as.numeric(threshold)
}

# Ensure daemon cleanup
stop_daemons <- function() {
  if (requireNamespace("mirai", quietly = TRUE)) {
    tryCatch({
      mirai::daemons(0)
      logger::log_info("All mirai daemons stopped")
    }, error = function(e) {
      logger::log_error("Failed to stop mirai daemons: {conditionMessage(e)}")
    })
  }
}
