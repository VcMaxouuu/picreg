# =============================================================================
# Two-panel figure (grayscale, line-styles only):
#  (a) gradient of l_n at 0          -> sup-norm law depends on the censoring p0
#  (b) gradient of sqrt(l_n) at 0    -> p0 cancels (pivotal)
# =============================================================================
library(ggplot2)
library(patchwork)
set.seed(1)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

## ---------------------------- parameters ----------------------------------
p  <- 50
M  <- 5000
n0 <- 200
p0_grid <- c(0.3, 0.5, 0.7, 0.9)

## ----------------- fixed centered, standardized design --------------------
make_design <- function(n, p) {
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2, colMeans(X), "-")
  s <- sqrt(colMeans(X^2)); s[s == 0] <- 1
  sweep(X, 2, s, "/")
}

## ----- sup-norm of grad l_n at 0 (no root) and grad sqrt(l_n) at 0 ---------
stats_one <- function(X, y, delta) {
  n   <- nrow(X)
  ord <- order(y)
  Xs  <- X[ord, , drop = FALSE]; ds <- delta[ord]
  sizes  <- n:1
  suffix <- apply(Xs[n:1, , drop = FALSE], 2, cumsum)[n:1, , drop = FALSE]
  xbar   <- suffix / sizes
  num    <- max(abs(colSums(ds * (Xs - xbar)) / n))
  B_n    <- 2 * sqrt(sum(ds * log(sizes)) / n)
  c(noroot = sqrt(n) * num,
    root   = 2 * sqrt(n * log(n)) * (num / B_n))
}
sample_stats <- function(X, p0, M) {
  n <- nrow(X); crate <- (1 - p0) / p0
  nr <- numeric(M); rr <- numeric(M)
  for (m in seq_len(M)) {
    Tt <- rexp(n, 1); Cc <- rexp(n, crate)
    s <- stats_one(X, pmin(Tt, Cc), as.numeric(Tt <= Cc))
    nr[m] <- s["noroot"]; rr[m] <- s["root"]
  }
  list(noroot = nr, root = rr)
}
X0 <- make_design(n0, p)
S0 <- lapply(p0_grid, function(p0) sample_stats(X0, p0, M))

## ------------------------- long data frames -------------------------------
mk_long <- function(key)
  do.call(rbind, lapply(seq_along(S0), function(i)
    data.frame(value = as.numeric(S0[[i]][[key]]),
               p0 = factor(p0_grid[i], levels = p0_grid))))
d_no <- mk_long("noroot")
d_rt <- mk_long("root")

## ------------------------------ plotting ----------------------------------
lt_vals <- c("solid", "dashed", "dotted", "dotdash")
base_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "plain"))

# panel (a): no legend
p_a <- ggplot(d_no, aes(value, linetype = p0)) +
  geom_density(linewidth = 0.7) +
  scale_linetype_manual(values = lt_vals, name = expression(p[0])) +
  labs(title = "(a) Without normalisation: the law depends on censoring",
       x = expression(sqrt(n) ~ "||" * nabla * italic(l)[n](0) * "||"[infinity]),
       y = "density") +
  base_theme +
  theme(legend.position = "none")

# panel (b): classic boxed legend, inside top-right
p_b <- ggplot(d_rt, aes(value, linetype = p0)) +
  geom_density(linewidth = 0.7, key_glyph = "path") +
  scale_linetype_manual(values = lt_vals, name = expression(p[0])) +
  labs(title = "(b) Square-root normalisation: invariance to censoring",
       x = expression(2 * sqrt(n * log(n)) ~ "||" * nabla * sqrt(italic(l)[n])(0) * "||"[infinity]),
       y = "density") +
  base_theme +
  theme(legend.position      = c(0.98, 0.98),   # inside, top-right
        legend.justification = c(1, 1),
        legend.key.width     = unit(1.2, "cm"),
        legend.title         = element_text(size = 10),
        legend.text          = element_text(size = 9))

fig <- p_a + p_b      # no shared legend
ggsave(file.path(script_dir, "results/pivotality_graph.pdf"), fig,
       width = 9.5, height = 4.0)
cat("Figure saved to:", normalizePath(file.path(script_dir, "pivotality_graph.pdf")), "\n")