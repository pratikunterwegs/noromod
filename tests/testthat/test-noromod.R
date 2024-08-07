#### Basic tests for equivalence and correctness ####

# prepare initial conditions from example in Readme
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

# add contact matrix to pop
params <- default_parameters()
params[["contacts"]] <- uk_contact_rate_matrix

test_that("Basic expectations for R and Rcpp ODE systems for deSolve", {
  # pick a random timestep
  tstep <- sample(seq(100), 1)
  noromod_r <- norovirus_model_r(
    t = tstep,
    state = init_mat,
    parameters = params
  )

  noromod_cpp <- norovirus_model_cpp(
    t = tstep,
    state = init_mat,
    parameters = params
  )

  # expect equivalence for a single timestep
  expect_identical(
    noromod_cpp, noromod_r,
    tolerance = 1e-6
  )

  # expect identical integration with deSolve::lsoda
  times <- seq(1000)
  expect_no_condition(
    deSolve::lsoda(
      y = init_mat,
      times = times,
      func = norovirus_model_r,
      parms = params
    )
  )
  expect_no_condition(
    deSolve::lsoda(
      y = init_mat,
      times = times,
      func = norovirus_model_cpp,
      parms = params
    )
  )

  # expect values over time are identical (allow greater tolerance)
  expect_identical(
    deSolve::lsoda(
      y = init_mat,
      times = times,
      func = norovirus_model_cpp,
      parms = params
    ),
    deSolve::lsoda(
      y = init_mat,
      times = times,
      func = norovirus_model_r,
      parms = params
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
      state = init_mat,
      parameters = params
    )
  )
  expect_snapshot(
    norovirus_model_cpp(
      t = 100,
      state = init_mat,
      parameters = params
    )
  )
})

test_that("Expectations and equivalence for Boost C++ ODE system", {
  # run model and expect no conditions for a short run
  time_end <- 1100
  expect_no_condition(
    noromod_cpp_boost(
      initial_conditions = init_mat,
      params = params, time_end = time_end, increment = 1
    )
  )
  data <- noromod_cpp_boost(
    initial_conditions = init_mat,
    params = params, time_end = 1100, increment = 1
  )

  # expect output is a list
  expect_vector(data, list())

  # expect all values numeric and greater than 0
  expect_true(
    all(
      vapply(data[["x"]], FUN = function(x) {
        all(is.numeric(x)) && all(x >= 0)
      }, FUN.VALUE = TRUE)
    )
  )

  # expect equivalence with R model
  data_r <- deSolve::lsoda(
    y = init_mat,
    times = seq(0, time_end), # add zero to times
    func = norovirus_model_r,
    parms = params
  )

  # expect times are identical
  expect_identical(
    data[["time"]],
    data_r[, "time"]
  )

  # expect values are identical
  # convert output to comparable form
  data_as_matrix <- matrix(data[["time"]], nrow = length(data[["time"]]))
  data_as_matrix <- cbind(
    data_as_matrix, do.call(rbind, lapply(data[["x"]], as.vector))
  )

  # expect equal with high tolerance for solver differences
  expect_equal(
    data_as_matrix[1:10, 1:10],
    data_r[1:10, 1:10],
    ignore_attr = TRUE,
    tolerance = 1
  )
})
