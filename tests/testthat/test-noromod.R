#### Basic tests for equivalence and correctness ####

test_that("Basic expectations", {
  # seven compartments
  uk_pop <- epidemics::population(
    name = "UK population",
    contact_matrix = matrix(1, 4, 4),
    demography_vector = 67e6 * rep(0.25, 4),
    initial_conditions = matrix(
      c(1 - 1e-6, 0, 1e-6, 0, 0, 0, 0),
      nrow = 4, ncol = 7, byrow = TRUE
    )
  )

  noromod_r <- norovirus_model_r(
    t = 1,
    state = uk_pop$initial_conditions * uk_pop$demography_vector,
    default_parameters(population = uk_pop)
  )

  noromod_cpp <- norovirus_model_cpp(
    t = 1,
    state = uk_pop$initial_conditions * uk_pop$demography_vector,
    default_parameters(population = uk_pop)
  )

  # expect integration with lsoda
  times <- seq(1000)
  expect_no_condition(
    deSolve::lsoda(
      y = uk_pop$initial_conditions * uk_pop$demography_vector,
      times = times,
      func = norovirus_model_r,
      parms = default_parameters(population = uk_pop)
    )
  )
  expect_no_condition(
    deSolve::lsoda(
      y = uk_pop$initial_conditions * uk_pop$demography_vector,
      times = times,
      func = norovirus_model_cpp,
      parms = default_parameters(population = uk_pop)
    )
  )
  # expect identical
  expect_identical(
    noromod_cpp, noromod_r
  )
  expect_vector(
    noromod_cpp, list()
  )
})
