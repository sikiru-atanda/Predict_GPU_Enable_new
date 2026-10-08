
# Metrics registry for thread-safe aggregation
MetricsRegistry <- R6::R6Class("MetricsRegistry",
                           public = list(
                             queue_depth = 0,
                             task_times = numeric(0),
                             cpu_usage = numeric(0),
                             gpu_usage = numeric(0),
                             mem_usage = numeric(0),

                             initialize = function() {
                               self$queue_depth <- 0
                               self$task_times <- numeric(0)
                               self$cpu_usage <- numeric(0)
                               self$gpu_usage <- numeric(0)
                               self$mem_usage <- numeric(0)
                             },

                             update = function(queue_depth = NA_integer_, task_time = NA_real_, cpu = NA_real_, gpu = NA_real_, mem = NA_real_) {
                               if (!is.na(queue_depth)) self$queue_depth <- queue_depth
                               if (!is.na(task_time)) self$task_times <- c(self$task_times, task_time)
                               if (!is.na(cpu)) self$cpu_usage <- c(self$cpu_usage, cpu)
                               if (!is.na(gpu)) self$gpu_usage <- c(self$gpu_usage, gpu)
                               if (!is.na(mem)) self$mem_usage <- c(self$mem_usage, mem)
                               invisible(self)
                             },

                             get_metrics = function() {
                               list(
                                 queue_depth = self$queue_depth,
                                 avg_task_time = if (length(self$task_times) > 0) mean(self$task_times, na.rm = TRUE) else NA_real_,
                                 avg_cpu_usage = if (length(self$cpu_usage) > 0) mean(self$cpu_usage, na.rm = TRUE) else NA_real_,
                                 avg_gpu_usage = if (length(self$gpu_usage) > 0) mean(self$gpu_usage, na.rm = TRUE) else NA_real_,
                                 avg_mem_usage = if (length(self$mem_usage) > 0) mean(self$mem_usage, na.rm = TRUE) else NA_real_
                               )
                             }
                           )
)

# Global metrics registry
metrics_registry <- MetricsRegistry$new()

# Start Prometheus endpoint
# Opt-in (GP_PROMETHEUS_ENABLE=true) and bound to localhost unless
# GP_PROMETHEUS_HOST says otherwise: an unauthenticated endpoint must not be
# opened on every network interface by default.
start_metrics_endpoint <- function(port = as.integer(Sys.getenv("GP_PROMETHEUS_PORT", "9090")),
                                   host = Sys.getenv("GP_PROMETHEUS_HOST", "127.0.0.1")) {
  if (!tolower(Sys.getenv("GP_PROMETHEUS_ENABLE", "false")) %in% c("true", "1", "yes")) {
    return(NULL)
  }
  if (!requireNamespace("promises", quietly = TRUE) || !requireNamespace("httpuv", quietly = TRUE)) {
    logger::log_warn("promises or httpuv not available, skipping metrics endpoint")
    return(NULL)
  }
  tryCatch({
    server <- httpuv::startServer(host, port, list(
      call = function(req) {
        metrics <- metrics_registry$get_metrics()
        response <- paste0(
          "# HELP queue_depth Number of tasks in mirai queue\n",
          "# TYPE queue_depth gauge\n",
          "queue_depth ", metrics$queue_depth %||% 0, "\n",
          "# HELP avg_task_time Average task duration in seconds\n",
          "# TYPE avg_task_time gauge\n",
          "avg_task_time ", metrics$avg_task_time %||% 0, "\n",
          "# HELP avg_cpu_usage Average CPU usage percentage\n",
          "# TYPE avg_cpu_usage gauge\n",
          "avg_cpu_usage ", metrics$avg_cpu_usage %||% 0, "\n",
          "# HELP avg_gpu_usage Average GPU usage percentage\n",
          "# TYPE avg_gpu_usage gauge\n",
          "avg_gpu_usage ", metrics$avg_gpu_usage %||% 0, "\n",
          "# HELP avg_mem_usage Average memory usage in MB\n",
          "# TYPE avg_mem_usage gauge\n",
          "avg_mem_usage ", metrics$avg_mem_usage %||% 0, "\n"
        )
        list(status = 200L, headers = list("Content-Type" = "text/plain"), body = response)
      }
    ))
    logger::log_info("Prometheus endpoint started on port {port}")
    server
  }, error = function(e) {
    logger::log_error("Failed to start Prometheus endpoint: {conditionMessage(e)}")
    NULL
  })
}

# Stop Prometheus endpoint
stop_metrics_endpoint <- function(server) {
  if (!is.null(server)) {
    tryCatch({
      httpuv::stopServer(server)
      logger::log_info("Prometheus endpoint stopped")
    }, error = function(e) {
      logger::log_error("Failed to stop Prometheus endpoint: {conditionMessage(e)}")
    })
  }
}
