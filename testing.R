## Data
library('rrBLUP')
data('wheat', package = 'BGLR')
n_ps <- 100L
n_max_ts <- 200L
h2 <- 0.2

X <- wheat.X
n <- nrow(X)
m <- ncol(X)

alpha <- rnorm(m)
tbv <- drop(X %*% alpha)
sigma_e <- sqrt(var(tbv) * (1 - h2) / h2)
y <- rnorm(n = length(tbv), mean = tbv, sd = sigma_e)
G <- rrBLUP::A.mat(X = 2 * (X - 0.5))
ix_ps <- sort(sample.int(n = n, size = n_ps, replace = FALSE))
ix_ts <- sort(sample(x = setdiff(x = seq_len(n), y = ix_ps), size = n_max_ts, replace = FALSE))
tbv_ps <- tbv[ix_ps]

position <- rep(T, n_max_ts)

## fn <- function(particle) {
##   position <- particle$position
##   if (!any(position))
##     return (-Inf)

##   ix_ts_act <- ix_ts[position]
##   ix <- c(ix_ps, ix_ts_act)
##   y_ts <- y
##   y_ts[ix_ps] <- NA_real_
##   fm <- rrBLUP::mixed.solve(y = y_ts[ix], K = G[ix, ix])
##   cor(tbv_ps, fm$u[seq_len(n_ps)])
## }

fn <- function(position, ix_ts_, ix_ps_, y_, G_, tbv_ps_) {
  n_ps <- length(tbv_ps_)
  if (!any(position))
    return (-Inf)

  ix_ts_act <- ix_ts_[position]
  ix <- c(ix_ps_, ix_ts_act)
  y_ts <- y_
  y_ts[ix_ps_] <- NA_real_
  fm <- rrBLUP::mixed.solve(y = y_ts[ix], K = G_[ix, ix])
  cor(tbv_ps_, fm$u[seq_len(n_ps)])
}

library(BPSO)
fm <- bpsoptim(par = rep(FALSE, n_max_ts),
               fn = fn,
               ix_ts_= ix_ts, ix_ps_ = ix_ps, y_= y, G_= G, tbv_ps_ = tbv_ps,
               control = list('maxit' = 30L, REPORT = 1L, v.max = 5))

str(fm$stats)

## u <- BPSO_mod(fun = fun,
##               n = n_max_ts,
##               n_iter = 50,
##               print_each = 5L)

## u <- BPSO(fun = fun,
##               n = n_max_ts,
##               n_iter = 50)


## n <- 100L ## dimension of the problem

## ## Parameters
## n_particles <- floor(10.0 + 2.0 * sqrt(n))
## n_iter <- 50L
## omega <- 1.0 / (2.0 * log(2.0))
## phi_p <- phi_g <- 1.0 + log(2.0)
## prob <- 0.0

## ## Heterozygosity
## ## evaluate_particle <- function(particle) {
## ##   p <- colMeans(X[particle$position, , drop = FALSE])
## ##   sum(p * (1.0 - p))
## ## }

## ## genetic variance
## va <- as.vector(X %*% alpha)
## evaluate_particle <- function(particle) {
##   v <- va[particle$position]
##   if (length(v) < 2L)
##     return(-Inf)
##   var(v)
## }

## var(X%*%alpha)
## u <- BPSO(fun = evaluate_particle, n = n, n_iter = 5000)
## u <- BPSO_mod(fun = evaluate_particle, n = n, n_iter = 5000)

## u$position
## var(X[u$position,]%*%alpha)


## optim()
## ?optim

## var(X%*%alpha)
## var(X[swarm$best_position,]%*%alpha)

## vapply(X = swarm$particles, FUN = function(x) x$best_value, FUN.VALUE = numeric(1L))

## u <- sapply(swarm$particles, function(x) x$position)
## cor(u)


##     ## pos <- lapply(X = swarm$particles, FUN = function(x) x$position)
##     ## ##  best_pos <- lapply(X = swarm$particles, FUN = function(x) x$best_position)
##     ## ##  mat <- do.call(what = rbind, args = c(list(swarm$best_position), best_pos, pos))
##     ## mat <- do.call(rbind, pos)

##     ## p <- colMeans(x = mat)
##     ## mat <- mat[, p * (1.0 - p) > sqrt(.Machine$double.eps)]
##     ## pca <- prcomp(x = as.data.frame(t(mat)), center = TRUE, scale. = TRUE)
##     ## ## print(pca)
##     ## ## plot(pca)
##     ## ## summary(pca)
##     ## ## str(pca$rotation)
##     ## biplot(pca, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
