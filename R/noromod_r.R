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
  
  # New structure (7 columns per section):
  # Non-vaccinated: S|E|Is|Ia|R|re_infections|new_infections
  # Vax1: V1|Ev1|Isv1|Iav1|Rv1|re_infections_v1|new_infections_v1
  # Vax2: V2|Ev2|Isv2|Iav2|Rv2|re_infections_v2|new_infections_v2
  # Plus 2 vaccination tracking columns at the end: new_vaccinations_v1|new_vaccinations_v2
  
  state <- matrix(state, nrow = n_age_groups, ncol = 23L)  # 7 + 7 + 7 + 2 = 23 columns
  
  # Non-vaccinated compartments (columns 1-7)
  susceptible <- state[, 1]
  exposed <- state[, 2]
  infect_symp <- state[, 3]
  infect_asymp <- state[, 4]
  recovered <- state[, 5]
  re_infections_track <- state[, 6]
  new_infections_track <- state[, 7]
  
  # Vaccine 1 compartments (columns 8-14)
  vax_1 <- state[, 8]
  exposed_v1 <- state[, 9]
  infect_symp_v1 <- state[, 10]
  infect_asymp_v1 <- state[, 11]
  recovered_v1 <- state[, 12]
  re_infections_v1_track <- state[, 13]
  new_infections_v1_track <- state[, 14]
  
  # Vaccine 2 compartments (columns 15-21)
  vax_2 <- state[, 15]
  exposed_v2 <- state[, 16]
  infect_symp_v2 <- state[, 17]
  infect_asymp_v2 <- state[, 18]
  recovered_v2 <- state[, 19]
  re_infections_v2_track <- state[, 20]
  new_infections_v2_track <- state[, 21]
  
  # Vaccination tracking (columns 22-23)
  new_vaccinations_v1_track <- state[, 22]
  new_vaccinations_v2_track <- state[, 23]
  
  # count the total population
  total_pop <- rowSums(state[, c(1:5, 8:12, 15:19)])
  
  # Parameters
  delta <- 1 / (parameters[["D_immun"]] * 365)
  w1 <- (parameters[["season_amp"]] / 100)
  w2_values <- (parameters[["season_offset"]] / 100)
  w3 <- (parameters[["season_amp_over65"]])
  q1 <- exp(parameters[["probT_under5"]])
  q2 <- exp(parameters[["probT_over5"]])
  q <- c(q1, q2, q2, q2)
  
  rho <- parameters[["rho"]]
  b <- parameters[["b"]]
  d <- parameters[["d"]]
  sigma <- parameters[["sigma"]][1]
  sigma_v1 <- parameters[["sigma"]][2]
  sigma_v2 <- parameters[["sigma"]][3]
  epsilon <- parameters[["epsilon"]]
  psi <- parameters[["psi"]]
  gamma <- parameters[["gamma"]]
  aging <- parameters[["aging"]]
  vacc_start <- parameters[["vacc_start"]]
  
  # Vaccination parameters
  if (t > vacc_start) {
    phi_1 <- parameters[["phi_1"]]
    phi_2 <- parameters[["phi_2"]]
    upsilon_1 <- parameters[["upsilon"]][1]
    upsilon_2 <- parameters[["upsilon"]][2]
  } else {
    phi_1 <- 0
    phi_2 <- 0
    upsilon_1 <- 0
    upsilon_2 <- 0
  }
  
  # Contact matrix
  cm <- parameters[["contacts"]]
  
  # Calculate seasonal term
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
  
  # Calculate infection potentials
  infection_potential <- q * seasonal_term * (
    cm %*% (
      (infect_symp + infect_symp_v1 + infect_symp_v2) +
        (infect_asymp + infect_asymp_v1 + infect_asymp_v2) * rho
    )
  )
  
  # apply seperate infection potential to the 4th age group
  infection_potential[4] <- infection_potential[4] * w3
  
  # Calculate all infections
  new_infections <- susceptible * infection_potential
  new_infections_v1 <- vax_1 * infection_potential
  new_infections_v2 <- vax_2 * infection_potential
  
  re_infections <- recovered * infection_potential
  re_infections_v1 <- recovered_v1 * infection_potential
  re_infections_v2 <- recovered_v2 * infection_potential
  
  # Calculate vaccinations
  new_vaccinations_v1_S <- phi_1 * susceptible
  new_vaccinations_v1_R <- phi_1 * recovered
  
  # Calculate births
  births <- c(b * sum(total_pop), rep(0, n_age_groups - 1))
  
  # Non-vaccinated compartment changes
  dS <- births + (delta * recovered) - new_infections - (d * susceptible) + 
    (aging %*% susceptible) - (phi_1 * susceptible) + (upsilon_1 * vax_1)
  dE <- new_infections - (epsilon * exposed) - (d * exposed) +
    (aging %*% exposed)
  dIs <- (epsilon * sigma * exposed) - (psi * infect_symp) - (d * infect_symp) +
    (aging %*% infect_symp)
  dIa <- (epsilon * (1 - sigma) * exposed) + (psi * infect_symp) -
    (gamma * infect_asymp) - (d * infect_asymp) + re_infections +
    (aging %*% infect_asymp)
  dR <- (gamma * infect_asymp) - (delta * recovered) - (d * recovered) -
    re_infections + (aging %*% recovered) - (phi_1 * recovered) + 
    (upsilon_1 * recovered_v1)
  
  # Vaccine 1 compartment changes
  dV1 <- (delta * recovered_v1) - new_infections_v1 - (d * vax_1) + 
    (aging %*% vax_1) - (phi_2 * vax_1) + (upsilon_2 * vax_2) + 
    (phi_1 * susceptible) - (upsilon_1 * vax_1) + 
    (recovered_v2 * (delta * upsilon_2))
  dEv1 <- new_infections_v1 - (epsilon * exposed_v1) -
    (d * exposed_v1) + (aging %*% exposed_v1)
  dIsv1 <- (epsilon * sigma_v1 * exposed_v1) - (psi * infect_symp_v1) -
    (d * infect_symp_v1) + (aging %*% infect_symp_v1)
  dIav1 <- (epsilon * (1 - sigma_v1) * exposed_v1) + (psi * infect_symp_v1) -
    (gamma * infect_asymp_v1) - (d * infect_asymp_v1) + re_infections_v1 +
    (aging %*% infect_asymp_v1)
  dRv1 <- (gamma * infect_asymp_v1) - (delta * recovered_v1) - 
    (d * recovered_v1) - re_infections_v1 + (aging %*% recovered_v1) + 
    (upsilon_2 * recovered_v2) - (phi_2 * recovered_v1) - 
    (upsilon_1 * recovered_v1) + (phi_1 * recovered)
  
  
  # Vaccine 2 compartment changes
  dV2 <- (delta * recovered_v2) - new_infections_v2 - (d * vax_2) +
    (aging %*% vax_2) + (phi_2 * vax_1) - (upsilon_2 * vax_2)
  dEv2 <- new_infections_v2 - (epsilon * exposed_v2) -
    (d * exposed_v2) + (aging %*% exposed_v2)
  dIsv2 <- (epsilon * sigma_v2 * exposed_v2) - (psi * infect_symp_v2) -
    (d * infect_symp_v2) + (aging %*% infect_symp_v2)
  dIav2 <- (epsilon * (1 - sigma_v2) * exposed_v2) + (psi * infect_symp_v2) -
    (gamma * infect_asymp_v2) - (d * infect_asymp_v2) + re_infections_v2 +
    (aging %*% infect_asymp_v2)
  dRv2 <- (gamma * infect_asymp_v2) - (delta * recovered_v2) -
    (d * recovered_v2) - re_infections_v2 + (aging %*% recovered_v2) -
    (upsilon_2 * recovered_v2) + (phi_2 * recovered_v1) -
    (recovered_v2 * (delta * upsilon_2))
  
  
  # Return exactly 92 values (23 columns × 4 rows)
  return(list(c(
    # Non-vaccinated compartments and tracking (7 values)
    dS, dE, dIs, dIa, dR,
    re_infections, new_infections,
    
    # Vaccine 1 compartments and tracking (7 values)
    dV1, dEv1, dIsv1, dIav1, dRv1,
    re_infections_v1, new_infections_v1,
    
    # Vaccine 2 compartments and tracking (7 values)
    dV2, dEv2, dIsv2, dIav2, dRv2,
    re_infections_v2, new_infections_v2,
    
    # Vaccination tracking (2 values)
    new_vaccinations_v1_S, new_vaccinations_v1_R
  )))
}
