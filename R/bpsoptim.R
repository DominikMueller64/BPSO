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
#' Khanesar, M.A et al. (2007)
#'
#'  Default parameters are taken from:
#'  Cortez, P. (2014) (\url{http://www.springer.com/de/book/9783319082622})
#' }
#'
#' @examples
#' library('BPSO')
#'
#' set.seed(123L)
#' n <- 30L
#' m <- 100L
#' x <- replicate(n, sample(x = c(0L, 1L), size = m, replace = TRUE), simplify = FALSE)
#'
#' fn <- function (position) {
#'   if(!any(position))
#'     return(-Inf)
#'
#'   tmp <- do.call(what = rbind, args = x[position])
#'   tmp <- apply(X = tmp, MARGIN = 2L, FUN = var)
#'   val <- sum(tmp, na.rm = TRUE)
#'   if (is.na(val)) return(-Inf)
#'   val
#' }
#'
#' fm <- bpsoptim(par = rep(FALSE, m), fn = fn,
#'                control = list(maxit = 100L, REPORT = 5L))
#'
#' str(fm)
#'
#' @importFrom stats runif
#' @import futile.logger
#' @export
bpsoptim <- function(par, fn, ..., control = list(), debug = FALSE) {

  ## Get loggers
  flog.logger('console', INFO) ## logger for informative messages

  logfile <- 'bpso.log'
  if (file.exists(logfile)) file.remove(logfile) ## Remove any previous logfile
  flog.logger('logfile', DEBUG, appender = appender.file(logfile)) ## logger for debugging

  if (!debug) {
    flog.threshold(0L, name = 'debug')
  }

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

        flog.info('iteration = %4d | value = %.3f | P(change) = %1.3f',
                  it, swarm$best_value, prob_bit_change,
                  name = 'info')

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
