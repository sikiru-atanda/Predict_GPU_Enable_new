gp_set_seed <- function(seed) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    current_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (!is.integer(current_seed) || length(current_seed) <= 1L) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }

  set.seed(seed)
  invisible(NULL)
}

gp_model_task_seed <- function(base_seed = 123L, task_index = 1L) {
  base_seed <- suppressWarnings(as.numeric(base_seed))
  task_index <- suppressWarnings(as.numeric(task_index))
  if (length(base_seed) != 1L || !is.finite(base_seed)) base_seed <- 123
  if (length(task_index) != 1L || !is.finite(task_index) || task_index < 1) {
    stop("`task_index` must be a positive integer.", call. = FALSE)
  }
  # Calculate in double precision before the modulus to avoid integer
  # overflow, then return the scalar integer required by set.seed().  The
  # historical default remains 10123 for task 1, while a user-supplied
  # random_state now genuinely changes true-prediction Bayesian/ML/DL draws.
  seed <- (base_seed + floor(task_index) * 10000) %% .Machine$integer.max
  if (!is.finite(seed) || seed <= 0) seed <- 1
  as.integer(seed)
}

gp_clear_invalid_random_seed <- function(previous_seed = NULL, previous_seed_existed = FALSE) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    return(invisible(NULL))
  }

  current_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (!is.integer(current_seed) || length(current_seed) <= 1L) {
    if (isTRUE(previous_seed_existed) && is.integer(previous_seed) && length(previous_seed) > 1L) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
  invisible(NULL)
}

# Evaluate `expr` with a reproducible MCMC/model stream: Mersenne-Twister
# seeded with `seed` (the same in a sequential run and in any parallel
# worker), then restore the caller's RNG kind and stream. A NULL/NA seed
# evaluates `expr` without touching the RNG.
gp_with_pinned_seed <- function(seed, expr) {
  seed <- suppressWarnings(as.numeric(seed)[1L])
  if (!length(seed) || !is.finite(seed)) {
    return(expr)
  }
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(abs(seed) %% .Machine$integer.max))
  expr
}