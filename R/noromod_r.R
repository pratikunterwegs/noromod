
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

  # hardcoding the number of columns
  state <- matrix(
    state,
    nrow = n_age_groups, ncol = 7L
  )
  susceptible <- state[, 1]
  exposed <- state[, 2]
  infect_symp <- state[, 3]
  infect_asymp <- state[, 4]
  recovered <- state[, 5]

  # count the total population as the sum of the first five columns
  # hardcoding the number of epidemiological compartments - S,E,Is,Ia,R
  total_pop <- rowSums(state[, seq(5)])

  # some parameters
  delta <- 1 / (parameters[["D_immun"]] * 365)
  w1 <- (parameters[["season_amp"]] / 100)
  w2 <- (parameters[["season_offset"]] / 10)
  q1 <- (parameters[["probT_under5"]] / 10)
  q2 <- (parameters[["probT_over5"]] / 100)
  q <- c(q1, q2, q2, q2)

  # more parameters
  rho <- parameters[["rho"]]
  b <- parameters[["b"]]*60286751
  d <- parameters[["d"]]
  sigma <- parameters[["sigma"]]
  epsilon <- parameters[["epsilon"]]
  psi <- parameters[["psi"]]
  gamma <- parameters[["gamma"]]
  aging <- parameters[["aging"]]
  
  # contact matrix
  cm <- parameters[["contacts"]]

  seasonal_term <- seasonal_forcing(t = t, w1 = w1, w2 = w2)
  infection_potential <- q * seasonal_term * (
    cm %*% (infect_symp + infect_asymp * rho)
  )
  new_infections <- susceptible * infection_potential
  re_infections <- recovered * infection_potential

  dS[1] <- dS[1] + b # Adding births only to the first age group
  dS <- delta * recovered - new_infections - d * susceptible + aging %*% susceptible
  dE <- new_infections - epsilon * (1 - sigma) * exposed - epsilon *
    sigma * exposed - d * exposed + aging %*% exposed
  dIs <- epsilon * sigma * exposed - psi * infect_symp - d * infect_symp + aging %*% infect_symp
  dIa <- epsilon * (1 - sigma) * exposed + psi * infect_symp - gamma *
    infect_asymp - d * infect_asymp + re_infections + aging %*% infect_asymp
  dR <- gamma * infect_asymp - delta * recovered - d * recovered - re_infections + aging %*% recovered

  return(list(c(dS, dE, dIs, dIa, dR, re_infections, new_infections)))
}
