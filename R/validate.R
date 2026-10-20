# Argument checks shared by the user-facing functions. Each helper returns
# its (possibly coerced) argument or stops with a message naming it.

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  x
}

.check_positive_int <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != round(x))
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  as.integer(x)
}

.check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop(sprintf("`%s` must be a single positive number.", name),
         call. = FALSE)
  as.numeric(x)
}

.check_alpha <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x >= 1)
    stop(sprintf("`%s` must be a single number in (0, 1).", name),
         call. = FALSE)
  as.numeric(x)
}

.check_seed <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(sprintf("`%s` must be NULL or a single number.", name),
         call. = FALSE)
  x
}

# Evaluate `expr` with the RNG seeded by `seed`, then restore the caller's
# RNG state, so that a seeded PDB calibration leaves the global stream
# untouched. `seed = NULL` evaluates `expr` on the current stream.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = env)
    else rm(".Random.seed", envir = env)
  }, add = TRUE)
  set.seed(seed)
  expr
}
