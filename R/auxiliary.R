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

    if (is.na(particle$value))
      flog.error('The function value is NA!', name = 'console')

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
