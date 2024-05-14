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
    sigma = c(0.82, 0.41, 0.41),
    phi_1 = c(1e-4, 0, 0, 1e-4), # vector, one value per age group
    phi_2 = c(0, 0, 0, 0),
    upsilon = c(4.4, 0),
    rho = 0.05,
    season_amp = 3.6,
    season_offset = c(5.76, 0, 0, 0, 0, 0, 0),
    # NOTE: only need change points, i.e., final values of each season
    # alternative season_change_points =
    # c(8580, 8944, 9315, 9679, 10043, 10407, 10771) # nolint
    season_change_points = c(11000, 0, 0, 0, 0, 0, 0),
    D_immun = 4.4,
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
    aging = aging / 365
  )

  params
}
