initialize_swarm <- function(par, fn, s) {
  n <- length(par)
  swarm <- list()
  swarm$particles <- replicate(n = s, simplify = FALSE, expr = {
    position <- best_position <- sample(x = c(FALSE, TRUE), size = n, replace = TRUE)
    velocity <- velocity_0 <- velocity_1 <- numeric(length = n)
    value <- best_value <- -Inf
    particle <- list('position' = position,
                     'velocity_0' = velocity_0,
                     'velocity_1' = velocity_1,
                     'velocity' = velocity,
                     'value' = value,
                     'best_position' = best_position,
                     'best_value' = best_value)

    particle$value <- particle$best_value <- fn(particle$position)
    particle
  })

  swarm$best_value <- -Inf
  swarm$best_position <- logical(n)
  swarm
}


update_swarm <- function(fn, swarm,  w, c_p, c_g) {

  for (i in seq_along(swarm$particles)) {
    particle <- swarm$particles[[i]]
    particle$value <- fn(particle$position)
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
    c1r1 <- runif(n = n, min = 0, max = c_p)
    c2r2 <- runif(n = n, min = 0, max = c_g)
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

    particle$velocity_1 <- w * particle$velocity_1 + d11 + d12
    particle$velocity_0 <- w * particle$velocity_0 + d01 + d02

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

#' @title Binary Particle Swarm Optimization
#'
#' @description This is an implementation of the binary particle swarm optimization as
#' described in Khanesar et al. (2007).
#'
#' @param par Vector with length defining the dimensionality of the
#'            optimization problem. Providing actual values of \code{par} are not
#'            necessary (\code{NA} is just fine). If values are provided, the first
#'            particle will be initialized to the position provided by \code{par}.
#'
#' @param fn A function to be minimized (or maximized), with first argument the vector
#'           of parameters over which minimization is to take place.  It should return
#'           a scalar result.
#'
#' @param \dots Further arguments to be passed to \code{fn}.
#'
#' @param control A list of control parameters. See \dQuote{Details}.
#'
#' @details
#' By default this function performs minimization using a binary particle swarm
#' algorithm, but it will maximize if \code{control$fnscale} is negative.
#'
#' The \code{control} argument is a list that can supply any of the following components:
#' \describe{
#'   \item{trace:}{Logical; if \code{TRUE}, tracing information on
#'     the progress of the optimization is produced. Defaults to \code{TRUE}.
#'   }
#'  \item{fnscale:}{An overall scaling to be applied to the value of \code{fn}
#'    during optimization. If negative, turns the problem into a maximization problem.
#'    Optimization is performed on \code{fn(par)/fnscale}. Defaults to \code{1}.
#'  }
#'  \item{maxit:}{
#'    The maximum number of iterations. Defaults to \code{1000}.
#'  }
#'  \item{REPORT:}{
#'    The frequency for reports if \code{control$trace} is
#'    positive. Defaults to \code{10}.
#'  }
#'  \item{trace.stats:}{Logical; if \code{TRUE} statistics at every
#'    reporting step are collected and returned. Defaults to \code{FALSE}.
#'  }
#'  \item{s:}{
#'    The swarm size. Defaults to \code{floor(10+2*sqrt(length(par)))}.
#'  }
#'  \item{w:}{
#'    The inertia weight. Defaults to \code{1/(2*log(2))}.
#'  }
#'  \item{c.p:}{
#'    The local exploration constant. Defaults to \code{1+log(2)}.
#'  }
#'  \item{c.g:}{
#'    The global exploration constant. Defaults to \code{1+log(2)}.
#'  }
#'  \item{v.max:}{
#'    The maximal (euclidean) length of the velocity vector. Defaults to \code{NA}
#'    which disables clamping of the velocity. Currently not implemented!
#'  }
#'  \item{maxit.stagnate:}{
#'    The maximum number of iterations without improvement.
#'    Defaults to \code{Inf}.
#'  }
#' }
#' \value{
#'  A list with components:
#'  \item{par}{
#'     The best set of parameters found.
#'  }
#'  \item{value}{
#'    The value of \code{fn} corresponding to \code{par}.
#'  }
#' }
#' If \code{trace.stats} is \code{TRUE} additionally the component:
#' \item{stats}{
#'   A list of statistics collected at every reporting step with the following components:
#'   \describe{
#'     \item{\code{it}}{The iteration number}
#'     \item{\code{par}}{The currently best set of parameters.}
#'     \item{\code{value}}{The corresponding value of \code{fn}.}
#'     \item{\code{f}}{The current values of \code{fn} of each particle.}
#'     \item{\code{x}}{The current particle positions as a matrix.}
#'   }
#' }
#'
#' \references{
#' Algorithm is taken from:
#'
#' Khanesar, M.A et al. (2007)
#' \url{http://zeus.inf.ucv.cl/~bcrawford/MII-748-METAHEURISTICAS/Papers_BinaryMH/A%20novel%20BPSO.pdf}
#'
#'  Default parameters are taken from:
#'  Cortez, P. (2014)
#'  \url{http://www.springer.com/de/book/9783319082622}
#' }
#'
#' @export
bpsoptim <- function(par, fn, ..., control = list()) {

  n <- length(par)
  con <- list(trace = TRUE,
              REPORT = 10L,
              'fnscale' = 1.0, maxit = 1000L,
              s = floor(10.0 + 2.0 * sqrt(n)),
              w = 1.0 / (2.0 * log(2.0)),
              c.p = 1.0 + log(2.0),
              c.g = 1.0 + log(2.0),
              v.max = NA_real_,
              maxit.stagnate = Inf,
              trace.stats = FALSE)

  names_con <- names(con)
  con[(names_control <- names(control))] <- control
  if (length(no_names <- names_control[!(names_control %in% names_con)]))
    warning("unknown names in control: ", paste(no_names, collapse = ", "))

  ## Re-scale function
  fn1 <- function(par) fn(par, ...) / con$fnscale

  ## Re-name parameters
  trace <- con$trace > 0L
  REPORT <- con$REPORT
  trace_stats <- con$trace.stats

  ## Check v.max
  if (!is.na(con$v.max)) {
    .NotYetUsed('v.max', error = FALSE)
  }

  ## Initialize swarm
  swarm <- initialize_swarm(par = par, fn = fn1, s = con$s)

  ## Evolve swarm
  if (trace_stats) {
    trace_stats_ret <- list()
    tsc <- 1L
  }

  it <- 1L
  it_stagnate <- 0L
  while(TRUE) {

    ## Update particles
    best_value <- swarm$best_value
    swarm <- update_swarm(fn = fn1, swarm = swarm, w = con$w, c_p = con$c.p, c_g = con$c.g)

    cond <- isTRUE(all.equal(target = best_value, current = swarm$best_value))
    if (cond) {
      it_stagnate <- it_stagnate + 1L
    } else {
      it_stagnat <- 0L
    }

    ## Print trace
    if (trace) {
      if (it %% REPORT == 0L) {

        tmp <- lapply(X = swarm$particles, FUN = function(x) 1.0 / (1.0 + exp(-x$velocity)))
        prob_bit_change <- mean(Reduce(f = `+`, x = tmp) / length(tmp))

        message(sprintf("iteration = %4d | value = %.3f | P(change) = %1.3f",
                        it, swarm$best_value, prob_bit_change))
      }
    }

    ## Store trace stats
    if (trace_stats) {
      if (it %% REPORT == 0L) {
        trace_stats_ret[[tsc]] <-
          list('iteration' = it,
               'par' = swarm$best_position,
               'value' = swarm$best_value,
               'f' = vapply(swarm$particles, function(x) x$value, FUN.VALUE = numeric(1L)),
               'x' = do.call(rbind, lapply(swarm$particles, function(x) x$position))
               )
        tsc <- tsc + 1L
      }
    }

    ## Check termination condition
    if (it == con$maxit || it_stagnate == con$maxit.stagnate) {
      break
    }
    it <- it + 1L
  }

  list('par' = swarm$best_position,
       'value' = swarm$best_value,
       'stats' = if (trace_stats) trace_stats_ret else NULL
       )
}

## initialize_swarm <- function(n, n_particles, prob = 0.5, fun) {
##   swarm <- list()
##   swarm$particles <- replicate(n = n_particles, simplify = FALSE, expr = {
##     position <- best_position <- sample(x = c(FALSE, TRUE), size = n, replace = TRUE,
##                                         prob = c(1.0 - prob, prob))
##     velocity <- runif(n = n, min = 0, max = 1)
##     value <- best_value <- NA_real_
##     particle <- list('position' = position,
##                      'velocity' = velocity,
##                      'value' = NA_real_,
##                      'best_position' = best_position,
##                      'best_value' = NA_real_)

##     particle$value <- particle$best_value <- fun(particle)
##     particle
##   })

##   values <- vapply(X = swarm$particles, FUN = function(x) x$value, FUN.VALUE = numeric(1L))
##   swarm$best_value <- max(values)
##   swarm$best_position <- swarm$particles[[which.max(values)]]$position
##   swarm
## }


## update_swarm <- function(swarm,  w, c_p, c_g, fun) {

##   update_particle <- function(particle) {
##     n <- length(particle$position)
##     r_p <- runif(n = n, min = 0, max = 1)
##     r_g <- runif(n = n, min = 0, max = 1)
##     velocity_star <- w * particle$velocity +
##       c_p * r_p * (particle$best_position - particle$position) +
##       c_g * r_g * (swarm$best_position - particle$position)
##     s <- 1.0 / (1.0 + exp(-velocity_star))
##     r <- runif(n = n, min = 0, max = 1)
##     logsub <- r < s

##     ## Update position and value
##     particle$position[logsub] <- TRUE
##     particle$position[!logsub] <- FALSE
##     particle$value <- fun(particle)

##     ## Update best_position and best_value
##     if (particle$value > particle$best_value) {
##       particle$best_position <- particle$position
##       particle$best_value <- particle$value
##     }

##     particle
##   }

##   for (i in seq_along(swarm$particles)) {
##     swarm$particles[[i]] <- update_particle(particle = swarm$particles[[i]])
##   }

##   values <- vapply(X = swarm$particles, FUN = function(x) x$value, FUN.VALUE = numeric(1L))
##   max_ix <- which.max(values)
##   if (values[max_ix] > swarm$best_value) {
##     swarm$best_value <- values[max_ix]
##     swarm$best_position <- swarm$particles[[max_ix]]$position
##   }

##   swarm
## }


## BPSO <- function(fun, n, n_iter, n_particles = floor(10.0 + 2.0 * sqrt(n)),
##                  w = 1.0 / (2.0 * log(2.0)), c_p = 1.0 + log(2.0),
##                  c_g = 1.0 + log(2.0), verbose = TRUE) {

##   ## Initialization
##   swarm <- initialize_swarm(n = n, n_particles = n_particles, prob = 0.0, fun = fun)

##   ## Evolution
##   for (iter in seq_len(n_iter)) {

##     if (verbose && (iter == 1L || iter %% 5L == 0L)) {
##       cat(sprintf("iteration: %5d | best value = %.3f\n", iter, swarm$best_value))
##     }

##     swarm <- update_swarm(swarm = swarm, w = w, c_p = c_p,
##                           c_g = c_g, fun = fun)
##   }
##   list('value' = swarm$best_value,
##        'position' = swarm$best_position)
## }
