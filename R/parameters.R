#' @title Useful default parameters
#' @export
default_parameters <- function() {
  # prepare aging matrix
  ages <- c(4, 14, 64, 80)
  da <- diff(c(0, ages))
  length(ages)
  aging <- diag(-1 / da)
  aging[row(aging) - col(aging) == 1] <- 1 / utils::head(da, -1)
  # No ageing in last group - flow out via mortality rate

  upsilon <- c(4.5, 4.5)
  upsilon <- 1 / (upsilon * 365)
  upsilon[is.infinite(upsilon)] <- 0.0

  params <- list(
    target_coverage_under5 = 0.9,
    target_coverage_over65 = 0.75,
    days_to_target = 365,
    contacts = matrix(1),
    sigma = c(0.82, 0.82/2, 0.82/2),
    phi_1 = c(-log(1 - target_coverage_under5) / days_to_target,0,0,-log(1 - target_coverage_over65) / days_to_target), # vector, one value per age group
    phi_2 = c(0, 0, 0, 0),
    upsilon = upsilon,
    rho = 0.05,
    season_amp = 3.9,
    season_amp_over65 = 1.6,
    season_offset = c(8.4),
    # NOTE: only need change points, i.e., final values of each season
    # season_change_points = c(8580, 8944, 9315, 9679, 10043, 10407, 10771),
    season_change_points = c(11000, 0, 0, 0, 0, 0, 0),
    D_immun = 9.01,
    probT_under5 = log(0.195),
    probT_over5 = log(0.039),
    b = (11.4 / 1000) / 365,
    # background mortality must be a vector for C++ implementations
    # NOTE: this is not age-related mortality
    d = rep(0, length(ages)),
    epsilon = 1,
    psi = 1 / 2,
    gamma = 1 / 10,
    n_age_groups = 4,
    aging = aging / 365,
    vacc_start = 3535
  )

  params
}
