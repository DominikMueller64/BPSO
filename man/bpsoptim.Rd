% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bpsoptim.R
\name{bpsoptim}
\alias{bpsoptim}
\title{Binary Particle Swarm Optimization}
\usage{
bpsoptim(par, fn, ..., control = list(), debug = FALSE)
}
\arguments{
\item{par}{Vector with length defining the dimensionality of the
optimization problem. Providing actual values of \code{par} are not
necessary (\code{NA} is just fine). If values are provided, the first
particle will be initialized to the position provided by \code{par}.}

\item{fn}{A function to be minimized (or maximized), with first argument the vector
of parameters over which minimization is to take place.  It should return
a scalar result.}

\item{\dots}{Further arguments to be passed to \code{fn}.}

\item{control}{A list of control parameters. See \dQuote{Details}.}
}
\description{
This is an implementation of the binary particle swarm optimization as
described in Khanesar et al. (2007).
}
\details{
By default this function performs minimization using a binary particle swarm
algorithm, but it will maximize if \code{control$fnscale} is negative.

The \code{control} argument is a list that can supply any of the following components:
\describe{
  \item{trace:}{Logical; if \code{TRUE}, tracing information on
    the progress of the optimization is produced. Defaults to \code{TRUE}.
  }
 \item{fnscale:}{An overall scaling to be applied to the value of \code{fn}
   during optimization. If negative, turns the problem into a maximization problem.
   Optimization is performed on \code{fn(par)/fnscale}. Defaults to \code{1}.
 }
 \item{maxit:}{
   The maximum number of iterations. Defaults to \code{1000}.
 }
 \item{REPORT:}{
   The frequency for reports if \code{control$trace} is
   positive. Defaults to \code{10}.
 }
 \item{trace.stats:}{Logical; if \code{TRUE} statistics at every
   reporting step are collected and returned. Defaults to \code{FALSE}.
 }
 \item{s:}{
   The swarm size. Defaults to \code{floor(10+2*sqrt(length(par)))}.
 }
 \item{w:}{
   The inertia weight. Defaults to \code{1/(2*log(2))}.
 }
 \item{c.p:}{
   The local exploration constant. Defaults to \code{1+log(2)}.
 }
 \item{c.g:}{
   The global exploration constant. Defaults to \code{1+log(2)}.
 }
 \item{v.max:}{
   The maximal (euclidean) length of the velocity vector. Defaults to \code{NA}
   which disables clamping of the velocity. Currently not implemented!
 }
 \item{maxit.stagnate:}{
   The maximum number of iterations without improvement.
   Defaults to \code{Inf}.
 }
}
\value{
 A list with components:
 \item{par}{
    The best set of parameters found.
 }
 \item{value}{
   The value of \code{fn} corresponding to \code{par}.
 }
}
If \code{trace.stats} is \code{TRUE} additionally the component:
\item{stats}{
  A list of statistics collected at every reporting step with the following components:
  \describe{
    \item{\code{it}}{The iteration number}
    \item{\code{par}}{The currently best set of parameters.}
    \item{\code{value}}{The corresponding value of \code{fn}.}
    \item{\code{f}}{The current values of \code{fn} of each particle.}
    \item{\code{x}}{The current particle positions as a matrix.}
  }
}

\references{
Algorithm is taken from:
Khanesar, M.A et al. (2007)

 Default parameters are taken from:
 Cortez, P. (2014) (\url{http://www.springer.com/de/book/9783319082622})
}
}
\examples{
library('BPSO')

set.seed(123L)
n <- 30L
m <- 100L
x <- replicate(n, sample(x = c(0L, 1L), size = m, replace = TRUE), simplify = FALSE)

fn <- function (position) {
  if(!any(position))
    return(-Inf)

 tmp <- do.call(what = rbind, args = x[position])
 mean(apply(X = tmp, MARGIN = 2L, FUN = var))
}

fm <- bpsoptim(par = rep(FALSE, m), fn = fn,
               control = list(maxit = 100L, REPORT = 5L),
               debug = TRUE)

str(fm)

}
