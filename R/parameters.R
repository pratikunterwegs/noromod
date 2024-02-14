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

  params <- list(
    contacts = matrix(1),
    sigma = 0.72,
    rho = 0.070,
    season_amp = 7.45,
    season_offset = c(0.8, 0.1, 10, 0.8, 0.1, 0.1, 10),
    # NOTE: only need change points, i.e., final values of each season
    season_change_points = c(8580, 8944, 9315, 9679, 10043, 10407, 10771),
    D_immun = 6.8,
    probT_under5 = log(0.18),
    probT_over5 = log(0.036),
    b = (11.4 / 1000) / 365,
    # background mortality must be a vector for C++ implementations
    # NOTE: this is not age-related mortality
    d = rep(0, length(ages)),
    epsilon = 1,
    psi = 1 / 2,
    gamma = 1 / 10,
    n_age_groups = 4,
    aging = aging / 365
  )

  params
}
