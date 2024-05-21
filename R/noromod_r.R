#' @title Seasonal forcing coefficient for the norovirus model
#' @param t Double value for the time, taken as the ordinal day
#' @param w1 Double value for the first weighting factor
#' @param w2 Double value for the second weighting factor
#' @return Double value for the coefficient of seasonal forcing
#' @keywords internal
seasonal_forcing <- function(t, w1, w2) {
  z <- 1 + w1 * cos((2 * pi * t / 364) + w2)
  return(z)
}

#' @title Compartmental model function for a norovirus epidemic
#' @description A compartmental norovirus model with four age groups. This
#' function is intended to be passed to the \pkg{deSolve} function `lsoda`.
#' @param t Double value for the time.
#' @param state Vector of the initial state. The form of the vector should be
#' \eqn{X_i, for X \in S, E, I_s, I_a, R, I_new, I_re}, where values for each
#' of the compartments for each age group \eqn{i} are given consecutively.
#' @param parameters List for the parameters used in the simulation. See
#' [default_parameters()] for a list of parameters used in development.
#' @return A list with a single numeric vector of the same size as `state`,
#' suitable as output for [deSolve::lsoda()].
#' @export
norovirus_model_r <- function(t, state, parameters) {
  # prepare initial conditions
  n_age_groups <- parameters[["n_age_groups"]]

  # compartmental structure is repeated over array slices
  # nolint start
  # 1|2| 3| 4|5|      6|  7
  # S|E|Ia|Is|R|new_inf|re_inf
  # nolint end

  # hardcoding the dimensions although this could be specified in parameters
  state <- array(
    state,
    dim = c(n_age_groups, 7L, 3L) # 4 age groups, 7 states, 3 levels of vax
  )

  # count the total population as the sum of the first five columns
  # hardcoding the number of epidemiological compartments - S,E,Is,Ia,R
  total_pop <- rowSums(state[, seq(5L), ])

  # NOTE: parameters are handled here for clarity, some performance loss
  # is expected due to repeated calculations
  delta <- 1 / (parameters[["D_immun"]] * 365)
  w1 <- (parameters[["season_amp"]] / 100)
  w2_values <- (parameters[["season_offset"]] / 100)
  q1 <- exp(parameters[["probT_under5"]])
  q2 <- exp(parameters[["probT_over5"]])
  q <- c(q1, q2, q2, q2)

  # more parameters
  rho <- parameters[["rho"]]
  b <- parameters[["b"]] # 60286751
  d <- parameters[["d"]]
  sigma <- parameters[["sigma"]]
  epsilon <- parameters[["epsilon"]]
  psi <- parameters[["psi"]]
  gamma <- parameters[["gamma"]]
  aging <- parameters[["aging"]]

  # vaccination rates
  phi <- parameters[["phi"]]

  # waning rates
  upsilon <- 1 / (parameters[["upsilon"]] * 365)
  upsilon[is.infinite(upsilon)] <- 0

  # contact matrix
  cm <- parameters[["contacts"]]

  # calculate the current w2 using `season_offset_intervals` and `season_offset`
  # as a lookup table
  # use Position to search the list of intervals for the index where t is less
  # than the change point. Note default behaviour of Position() is to return the
  # first index that satisfies the condition
  # NOTE: t taken from function scope, this is the simulation time
  w2_current <- w2_values[
    Position(
      f = function(x) t <= x, x = parameters[["season_change_points"]]
    )
  ]

  seasonal_term <- seasonal_forcing(t = t, w1 = w1, w2 = w2_current)
  infection_potential <- q * seasonal_term * (
    cm %*% (rowSums(state[, 3L, ]) + rowSums(state[, 4L, ]) * rho)
  )

  # calculate new infections by vax status
  new_infections <- state[, 1, ] * c(infection_potential)

  # calculate re-infections by vax status
  re_infections <- state[, 5, ] * c(infection_potential)

  # Calculate births for addition only to the first age group
  # TODO: reconsider how births are calculated. Consider `b * [14 - 65]` only?
  # TODO: consider adding `b` to ageing matrix as inflow rate?
  births <- c(b * sum(total_pop), rep(0, n_age_groups - 1))

  # create empty array of the dimensions of state
  dState <- array(0, dim = dim(state))

  # waning immunity: δ*R + aging out: aging * S
  # change in susceptibles
  dState[, 1, ] <- (delta * state[, 5L, ]) + (aging %*% state[, 1L, ]) -
    new_infections
  dState[, 1, 1] <- dState[, 1, 1] + births

  # account for vaccinations and vaccine immunity waning in and out
  dS_vax_out <- state[, 1, ] * phi
  dS_waning_out <- state[, 1, ] * upsilon

  dState[, 1, ] <- dState[, 1, ] - (dS_vax_out + dS_waning_out) +
    dS_vax_out[, c(3, 1, 2)] + dS_waning_out[, c(2, 3, 1)]

  # direct waning from RV1 and RV2
  dState[, 1, c(1, 2)] <- dState[, 1, c(1, 2)] +
    (state[, 5, c(2, 3)] * delta * gamma)

  # exposed to infectious: ɛ*E + aging out: aging * E
  # change in exposed
  dState[, 2, ] <- -(epsilon * state[, 2L, ]) + (aging %*% state[, 2L, ]) +
    new_infections

  # symptomatic to asymptomatic: ψ*Is + # aging out: aging * Is +
  # exposed to infectious symptomatic: ɛ*σ*E
  # change in infectious symptomatic
  dState[, 3, ] <- -(psi * state[, 3L, ]) + (aging %*% state[, 3L, ]) +
    (epsilon * state[, 2, ] %*% diag(sigma))

  # symptomatic to asymptomatic: ψ*Is - recovery: γ*Ia + aging out: aging * Ia
  # exposed to infectious asymptomatic: ɛ*(1-σ)*E
  # change in infectious asymptomatic
  dState[, 4, ] <- (psi * state[, 3L, ]) - (gamma * state[, 4L, ]) +
    (aging %*% state[, 4L, ]) + re_infections +
    (epsilon * state[, 2, ] %*% diag((1 - sigma)))

  # recovery: γ*Ia + aging out: aging * R
  # change in recovered
  dState[, 5, ] <- (gamma * state[, 4L, ]) + (aging %*% state[, 5L, ]) -
    re_infections

  # account for vaccinations and vaccine immunity waning in and out
  dR_vax_out <- state[, 5, ] * phi
  dR_waning_out <- state[, 5, ] * upsilon

  dState[, 5, ] <- dState[, 5, ] - (dR_vax_out + dR_waning_out) +
    dR_vax_out[, c(3, 1, 2)] + dR_waning_out[, c(2, 3, 1)]

  # direct waning to S and SV1
  dState[, 5, c(2, 3)] <- dState[, 5, c(2, 3)] -
    (state[, 5, c(2, 3)] * delta * gamma)

  # store new infections and re-infections
  dState[, 6, ] <- new_infections
  dState[, 7, ] <- re_infections

  # background mortality in all epi compartments in all vax strata
  dState[, 1:5, ] <- dState[, 1:5, ] - (state[, 1:5, ] * d)

  # return in the same order as state
  return(list(c(dState)))
}
