#' @title Useful default parameters
#' @param population An optional argument taking an object of class
#' `<population>` from the \pkg{epidemics}. Used to obtain a contact matrix.
#' package.
#' @export
default_parameters <- function(population = NULL) {
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
    season_offset = 0.08,
    D_immun = 6.8,
    probT_under5 = 1.8,
    probT_over5 = 3.6,
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
  if (!is.null(population)) {
    checkmate::assert_class(population, "population")
    contacts <- population$contact_matrix / population$demography_vector
    params[["contacts"]] <- contacts
  }
  params
}
