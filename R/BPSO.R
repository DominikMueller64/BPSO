initialize_swarm <- function(n, n_particles, prob = 0.5, fun) {
  swarm <- list()
  swarm$particles <- replicate(n = n_particles, simplify = FALSE, expr = {
    position <- best_position <- sample(x = c(FALSE, TRUE), size = n, replace = TRUE,
                                        prob = c(1.0 - prob, prob))
    velocity <- runif(n = n, min = 0, max = 1)
    value <- best_value <- NA_real_
    particle <- list('position' = position,
                     'velocity' = velocity,
                     'value' = NA_real_,
                     'best_position' = best_position,
                     'best_value' = NA_real_)

    particle$value <- particle$best_value <- fun(particle)
    particle
  })


  values <- vapply(X = swarm$particles, FUN = function(x) x$value, FUN.VALUE = numeric(1L))
  swarm$best_value <- max(values)
  swarm$best_position <- swarm$particles[[which.max(values)]]$position
  swarm
}


update_swarm <- function(swarm,  omega, phi_p, phi_g, fun) {

  update_particle <- function(particle) {
    n <- length(particle$position)
    r_p <- runif(n = n, min = 0, max = 1)
    r_g <- runif(n = n, min = 0, max = 1)
    velocity_star <- omega * particle$velocity +
      phi_p * r_p * (particle$best_position - particle$position) +
      phi_g * r_g * (swarm$best_position - particle$position)
    s <- 1.0 / (1.0 + exp(-velocity_star))
    r <- runif(n = n, min = 0, max = 1)
    logsub <- r < s

    ## Update position and value
    particle$position[logsub] <- TRUE
    particle$position[!logsub] <- FALSE
    particle$value <- fun(particle)

    ## Update best_position and best_value
    if (particle$value > particle$best_value) {
      particle$best_position <- particle$position
      particle$best_value <- particle$value
    }

    particle
  }

  for (i in seq_along(swarm$particles)) {
    swarm$particles[[i]] <- update_particle(particle = swarm$particles[[i]])
  }

  values <- vapply(X = swarm$particles, FUN = function(x) x$value, FUN.VALUE = numeric(1L))
  max_ix <- which.max(values)
  if (values[max_ix] > swarm$best_value) {
    swarm$best_value <- values[max_ix]
    swarm$best_position <- swarm$particles[[max_ix]]$position
  }

  swarm
}


#' @export
BPSO <- function(fun, n, n_iter, n_particles = floor(10.0 + 2.0 * sqrt(n)),
                 omega = 1.0 / (2.0 * log(2.0)), phi_p = 1.0 + log(2.0),
                 phi_g = 1.0 + log(2.0), verbose = TRUE) {

  ## Initialization
  swarm <- initialize_swarm(n = n, n_particles = n_particles, prob = 0.0, fun = fun)

  ## Evolution
  for (iter in seq_len(n_iter)) {

    if (verbose && (iter == 1L || iter %% 5L == 0L)) {
      cat(sprintf("iteration: %5d | best value = %.3f\n", iter, swarm$best_value))
    }

    swarm <- update_swarm(swarm = swarm, omega = omega, phi_p = phi_p,
                          phi_g = phi_g, fun = fun)
  }
  list('value' = swarm$best_value,
       'position' = swarm$best_position)
}


initialize_swarm_mod <- function(n, n_particles, fun) {
  swarm <- list()
  swarm$particles <- replicate(n = n_particles, simplify = FALSE, expr = {
    position <- best_position <- sample(x = c(FALSE, TRUE), size = n, replace = TRUE)
    ## velocity <- runif(n = n, min = 0, max = 1)
    velocity <- velocity_0 <- velocity_1 <- numeric(length = n)
    value <- best_value <- -Inf
    particle <- list('position' = position,
                     'velocity_0' = velocity_0,
                     'velocity_1' = velocity_1,
                     'velocity' = velocity,
                     'value' = value,
                     'best_position' = best_position,
                     'best_value' = best_value)

    particle$value <- particle$best_value <- fun(particle)
    particle
  })

  ## values <- vapply(X = swarm$particles, FUN = function(x) x$value, FUN.VALUE = numeric(1L))
  swarm$best_value <- -Inf
  ## swarm$best_position <- swarm$particles[[which.max(values)]]$position
  swarm$best_position <- logical(n)
  swarm
}


update_swarm_mod <- function(swarm,  omega, phi_p, phi_g, fun) {

  for (i in seq_along(swarm$particles)) {
    particle <- swarm$particles[[i]]
    particle$value <- fun(particle)
    n <- length(particle$position)
    ## Compare current performance to its best performance
    if (particle$value > particle$best_value) {
      particle$best_value <- particle$value
      particle$best_position <- particle$position
    }

    ## Compare current performance to globally best performance
    if (particle$value > swarm$best_value) {
      swarm$best_value <- particle$value
      swarm$best_position <- particle$position
    }

    ## Change velocities of the particle
    d11 <- d01 <- d12 <- d02 <- logical(length = n)
    c1r1 <- runif(n = n, min = 0, max = phi_p)
    c2r2 <- runif(n = n, min = 0, max = phi_g)
    pbp <- particle$best_position
    gbp <- swarm$best_position

    d11[pbp] <- c1r1[pbp]
    d01[pbp] <- -c1r1[pbp]
    d01[!pbp] <- c1r1[!pbp]
    d11[!pbp] <- -c1r1[!pbp]
    d12[gbp] <- c2r2[gbp]
    d02[gbp] <- -c2r2[gbp]
    d02[!gbp] <- c2r2[!gbp]
    d12[!gbp] <- -c2r2[!gbp]

    particle$velocity_1 <- omega * particle$velocity_1 + d11 + d12
    particle$velocity_0 <- omega * particle$velocity_0 + d01 + d02

    ## Change velocities of bit-change
    pp <- particle$position
    particle$velocity[!pp] <- particle$velocity_1[!pp]
    particle$velocity[pp] <- particle$velocity_0[pp]

    ## Update position and value
    s <- 1.0 / (1.0 + exp(-particle$velocity))
    logsub <- runif(n = n) < s
    particle$position[logsub] <- !particle$position[logsub]
    ## leave unchanged otherwise
    swarm$particles[[i]] <- particle
  }

  swarm
}


#' @export
BPSO_mod <- function(fun, n, n_iter, n_particles = floor(10.0 + 2.0 * sqrt(n)),
                     omega = 1.0 / (2.0 * log(2.0)), phi_p = 1.0 + log(2.0),
                     phi_g = 1.0 + log(2.0), verbose = TRUE) {

  ## Initialization
  swarm <- initialize_swarm_mod(n = n, n_particles = n_particles, fun = fun)

  ## Evolution
  for (iter in seq_len(n_iter)) {

    if (verbose && (iter == 1L || iter %% 5L == 0L)) {
      cat(sprintf("iteration: %5d | best value = %.3f\n", iter, swarm$best_value))
    }

    swarm <- update_swarm_mod(swarm = swarm, omega = omega, phi_p = phi_p,
                              phi_g = phi_g, fun = fun)
  }
  list('value' = swarm$best_value,
       'position' = swarm$best_position)
}



## ## Data
## data('wheat', package = 'BGLR')
## X <- wheat.X
## n <- nrow(X)
## ## n <- 100L ## dimension of the problem
## alpha <- rnorm(ncol(X))

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
