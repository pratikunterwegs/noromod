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

  tstep <- sample(seq(100), 1)
  noromod_r <- norovirus_model_r(
    t = tstep,
    state = uk_pop$initial_conditions * uk_pop$demography_vector,
    default_parameters(population = uk_pop)
  )

  noromod_cpp <- norovirus_model_cpp(
    t = tstep,
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
    noromod_cpp, noromod_r,
    tolerance = 1e-6
  )

  # expect values over time are identical (allow greater tolerance)
  expect_identical(
    deSolve::lsoda(
      y = uk_pop$initial_conditions * uk_pop$demography_vector,
      times = times,
      func = norovirus_model_cpp,
      parms = default_parameters(population = uk_pop)
    ),
    deSolve::lsoda(
      y = uk_pop$initial_conditions * uk_pop$demography_vector,
      times = times,
      func = norovirus_model_r,
      parms = default_parameters(population = uk_pop)
    ),
    tolerance = 1
  )

  expect_vector(
    noromod_cpp, list()
  )

  # expect snapshots
  expect_snapshot(
    norovirus_model_r(
      t = 100,
      state = uk_pop$initial_conditions * uk_pop$demography_vector,
      default_parameters(population = uk_pop)
    )
  )
  expect_snapshot(
    norovirus_model_cpp(
      t = 100,
      state = uk_pop$initial_conditions * uk_pop$demography_vector,
      default_parameters(population = uk_pop)
    )
  )
})

test_that("Model with Boost C++ solvers", {
  init <- c(
    3857263, 8103718, 42460865, 12374961,
    100, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0
  )
  init_mat <- matrix(
    init,
    nrow = 4, ncol = 7
  )

  age_groups <- c(0, 5, 15, 65)
  polymod <- socialmixr::polymod
  UK_structure <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(age_groups),
    symmetric = TRUE
  )

  # Symmetrical contact matrix
  uk_contact_rate_matrix <- as.matrix(UK_structure$matrix)
  demography <- UK_structure$demography$population

  uk_contact_rate_matrix <- t(t(uk_contact_rate_matrix) / demography)

  # seven compartments
  uk_pop <- epidemics::population(
    name = "UK population",
    contact_matrix = uk_contact_rate_matrix,
    demography_vector = demography,
    initial_conditions = matrix(
      c(1 - 1e-6, 0, 1e-6, 0, 0, 0, 0),
      nrow = 4, ncol = 7, byrow = TRUE
    )
  )

  # add contact matrix to pop
  params <- default_parameters(uk_pop)
  params[["contacts"]] <- uk_contact_rate_matrix

  # run model
  data <- noromod_cpp_boost(
    initial_conditions = matrix(init, nrow = 4, ncol = 7),
    params = params, time_end = 1100, increment = 1
  )

  # expect output is a list
  expect_vector(data, list())
  # expect dataframes
  data <- output_to_df(data)
  expect_s3_class(data, "data.frame")

  # expect all numeric and greater than 0
  expect_true(
    all(
      apply(data, 2, function(x) {
        all(is.numeric(x)) && all(x >= 0)
      })
    )
  )
})
