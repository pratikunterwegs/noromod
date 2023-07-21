
#' @title Useful default parameters
#' @param population An object of class `<population>` from the \pkg{epidemics}
#' package.
#' @export
default_parameters <- function(population) {
  checkmate::assert_class(population, "population")
  params <- list(
    contacts = matrix(1),
    sigma = 0.72,
    rho = 0.070,
    season_amp = 7.45,
    season_offset = 0.08,
    D_immun = 6.8,
    probT_under5 = 1.8,
    probT_over5 = 3.6,
    rho = 0.05,
    b = (11.4 / 1000) / 365,
    d = (11.4 / 1000) / 365,
    sigma = 0.5,
    epsilon = 1,
    psi = 1 / 2,
    gamma = 1 / 10,
    n_age_groups = 4
  )
  if (!missing(population)) {
    contacts <- population$contact_matrix / population$demography_vector
    params[["contacts"]] <- contacts
  }
  params
}
