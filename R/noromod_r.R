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

  # compartmental structure
  # S|E|Ia|Is|R|V1|Ev1|Iav1|Isv1|Rv1|V2|Ev2|Iav2|Isv2|Rv2

  # hardcoding the number of columns
  state <- matrix(
    state,
    nrow = n_age_groups, ncol = 21L
  )
  susceptible <- state[, 1]
  exposed <- state[, 2]
  infect_symp <- state[, 3]
  infect_asymp <- state[, 4]
  recovered <- state[, 5]

  vax_1 <- state[, 6]
  exposed_v1 <- state[, 7]
  infect_symp_v1 <- state[, 8]
  infect_asymp_v1 <- state[, 9]
  recovered_v1 <- state[, 10]

  vax_2 <- state[, 11]
  exposed_v2 <- state[, 12]
  infect_symp_v2 <- state[, 13]
  infect_asymp_v2 <- state[, 14]
  recovered_v2 <- state[, 15]

  # count the total population as the sum of the first five columns
  # hardcoding the number of epidemiological compartments - S,E,Is,Ia,R
  total_pop <- rowSums(state[, seq(15)])

  # some parameters
  delta <- 1 / (parameters[["D_immun"]] * 365)
  w1 <- (parameters[["season_amp"]] / 100)
  w2_values <- (parameters[["season_offset"]] / 100)
  q1 <- exp(parameters[["probT_under5"]])
  q2 <- exp(parameters[["probT_over5"]])
  q <- c(q1, q2, q2, q2)

  # more parameters
  rho <- parameters[["rho"]]
  b <- parameters[["b"]] #* 60286751
  d <- parameters[["d"]]
  sigma <- parameters[["sigma"]][1]
  sigma_v1 <- parameters[["sigma"]][2]
  sigma_v2 <- parameters[["sigma"]][3]
  epsilon <- parameters[["epsilon"]]
  psi <- parameters[["psi"]]
  gamma <- parameters[["gamma"]]
  aging <- parameters[["aging"]]

  # vaccination rates
  phi_1 <- parameters[["phi_1"]]
  phi_2 <- parameters[["phi_2"]]

  # waning rates
  upsilon_1 <- 1 / (parameters[["upsilon"]][1] * 365)
  upsilon_2 <- 1 / (parameters[["upsilon"]][2] * 365)

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
    cm %*% (
      (infect_symp + infect_symp_v1 + infect_symp_v2) +
        (infect_asymp + infect_asymp_v1 + infect_asymp_v2) * rho
    )
  )
  # calculate new infections by vax status
  new_infections <- susceptible * infection_potential
  new_infections_v1 <- vax_1 * infection_potential
  new_infections_v2 <- vax_2 * infection_potential

  # calculate re-infections by vax status
  re_infections <- recovered * infection_potential
  re_infections_v1 <- recovered_v1 * infection_potential
  re_infections_v2 <- recovered_v2 * infection_potential

  # Calculate births for addition only to the first age group
  # TODO: reconsider how births are calculated. Consider `b * [14 - 65]` only?
  # TODO: consider adding `b` to ageing matrix as inflow rate?
  births <- c(b * sum(total_pop), rep(0, n_age_groups - 1))

  # TODO: annotate transitions
  # compartmental transitions - non vaccinated
  dS <- births + (delta * recovered) - new_infections - (d * susceptible) +
    (aging %*% susceptible) - (phi_1 * susceptible) + (upsilon_1 * vax_1) + (recovered_v1 * (delta * upsilon_1))
  dE <- new_infections - epsilon * (1 - sigma) * exposed - epsilon *
    sigma * exposed - d * exposed + (aging %*% exposed)
  dIs <- epsilon * sigma * exposed - psi * infect_symp - d * infect_symp +
    (aging %*% infect_symp)
  dIa <- epsilon * (1 - sigma) * exposed + psi * infect_symp - gamma *
    infect_asymp - d * infect_asymp + re_infections + (aging %*% infect_asymp)
  dR <- gamma * infect_asymp - delta * recovered - d * recovered -
    re_infections + (aging %*% recovered) + (upsilon_1 * recovered_v1) -
    (phi_1 * recovered)

  # compartmental transitions - vax 1
  dV1 <- (delta * recovered_v1) - new_infections_v1 - (d * vax_1) +
    (aging %*% vax_1) - (phi_2 * vax_1) + (upsilon_2 * vax_2) +
    (phi_1 * susceptible) - (upsilon_1 * vax_1)
  dEv1 <- new_infections_v1 - epsilon * (1 - sigma_v1) * exposed_v1 - epsilon *
    sigma_v1 * exposed_v1 - d * exposed_v1 + (aging %*% exposed_v1)
  dIsv1 <- epsilon * sigma_v1 * exposed_v1 - psi * infect_symp_v1 -
    d * infect_symp_v1 + (aging %*% infect_symp_v1)
  dIav1 <- epsilon * (1 - sigma_v1) * exposed_v1 + psi * infect_symp_v1 -
    gamma * infect_asymp_v1 - d * infect_asymp_v1 + re_infections_v1 +
    (aging %*% infect_asymp_v1)
  dRv1 <- gamma * infect_asymp_v1 - delta * recovered_v1 - d * recovered_v1 -
    re_infections_v1 + (aging %*% recovered_v1) + (upsilon_2 * recovered_v2) -
    (phi_2 * recovered_v2) - (upsilon_1 * recovered_v1) + (phi_1 * recovered) - (recovered_v1 * (delta * upsilon_1))

  # compartmental transitions - vax 2
  dV2 <- (delta * recovered_v2) - new_infections_v2 - (d * vax_2) +
    (aging %*% vax_2) + (phi_2 * vax_1) - (upsilon_2 * vax_2)
  dEv2 <- new_infections_v2 - epsilon * (1 - sigma_v2) * exposed_v2 - epsilon *
    sigma_v2 * exposed_v2 - d * exposed_v2 + (aging %*% exposed_v2)
  dIsv2 <- epsilon * sigma_v2 * exposed_v2 - psi * infect_symp_v2 -
    d * infect_symp_v2 + (aging %*% infect_symp_v2)
  dIav2 <- epsilon * (1 - sigma_v2) * exposed_v2 + psi * infect_symp_v2 -
    gamma * infect_asymp_v2 - d * infect_asymp_v2 + re_infections_v2 +
    (aging %*% infect_asymp_v2)
  dRv2 <- gamma * infect_asymp_v2 - delta * recovered_v2 - d * recovered_v2 -
    re_infections_v2 + (aging %*% recovered_v2) - (upsilon_2 * recovered_v2) +
    (phi_2 * recovered_v1)

  return(list(c(
    dS, dE, dIs, dIa, dR,
    dV1, dEv1, dIsv1, dIav1, dRv1,
    dV2, dEv2, dIsv2, dIav2, dRv2,
    re_infections, new_infections,
    re_infections_v1, new_infections_v1,
    re_infections_v2, new_infections_v2
  )))
}
