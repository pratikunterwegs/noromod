// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <boost/numeric/odeint.hpp>
#include <cmath>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(epidemics)]]

/// @brief Return a coefficient of seasonal forcing
/// @param t Double value for the time, taken as the ordinal day
/// @param w1 Double value for the first weighting factor
/// @param w2 Double value for the second weighting factor
/// @return Double value for the coefficient of seasonal forcing
const double seasonal_forcing(const double &t, const double &w1,
                              const double &w2) {
  const double z = 1.0 + (w1 * std::cos((2.0 * M_PI * t / 364.0) + w2));
  return z;
}

//' @title Compartmental model function for a norovirus epidemic
//' @description A compartmental norovirus model with four age groups. This
//' function is intended to be passed to the \pkg{deSolve} function `lsoda`.
//' @param t Double value for the time.
//' @param state Vector of the initial state. The form of the vector should be
//' \eqn{X_i, for X \in S, E, I_s, I_a, R, I_new, I_re}, where values for each
//' of the compartments for each age group \eqn{i} are given consecutively.
//' @param parameters List for the parameters used in the simulation. See
//' [default_parameters()] for a list of parameters used in development.
//' @return A list with a single numeric vector of the same size as `state`,
//' suitable as output for [deSolve::lsoda()].
//' @export
// [[Rcpp::export]]
Rcpp::List norovirus_model_cpp(const double &t,
                               Eigen::VectorXd &state,  // NOLINT
                               const Rcpp::List &parameters) {
  // handle the initial conditions, which must be rehsaped into a matrix
  // matrix shape hardcoded to be 4 rows and 7 columns
  Eigen::MatrixXd state_matrix =
      Eigen::Map<Eigen::MatrixXd>(state.data(), 4L, 7L);

  // get current population size
  const Eigen::ArrayXd population_size =
      state_matrix.block(0, 0, 4, 5).rowwise().sum();

  // modify parameters
  const double delta = 1.0 / (Rcpp::as<double>(parameters["D_immun"]) * 365.0);
  const double w1 = Rcpp::as<double>(parameters["season_amp"]) / 100.0;
  const double w2 = Rcpp::as<double>(parameters["season_offset"]) / 10.0;
  const double q1 = Rcpp::as<double>(parameters["probT_under5"]) / 10.0;
  const double q2 = Rcpp::as<double>(parameters["probT_over5"]) / 100.0;
  // create contact matrix
  const Eigen::MatrixXd contacts = parameters["contacts"];

  // epi parameters
  const double rho = parameters["rho"];
  const double b = parameters["b"];
  const double d = parameters["d"];
  const double sigma = parameters["sigma"];
  const double epsilon = parameters["epsilon"];
  const double psi = parameters["psi"];
  const double gamma = parameters["gamma"];

  // prepare array for multiplication
  Eigen::ArrayXd param_(4L);
  param_ << q1, q2, q2, q2;

  // prepare seasonal forcing
  const double seasonal_term = seasonal_forcing(t, w1, w2);

  // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R
  // calculate new infections
  Eigen::ArrayXd sToE =
      param_ * seasonal_term * state_matrix.col(0).array() *
      (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho))).array();

  // calculate re-infections
  // recovered are column index 4 of the initial conditions
  Eigen::ArrayXd rToIa =
      param_ * seasonal_term * state_matrix.col(4).array() *
      (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho))).array();

  // compartmental transitions
  Eigen::ArrayXd births = (b * population_size);
  Eigen::ArrayXd rToS = (delta * state_matrix.col(4));
  Eigen::ArrayXd eToIa =
      ((1.0 - sigma) * epsilon) * state_matrix.col(1).array();
  Eigen::ArrayXd eToIs = (sigma * epsilon) * state_matrix.col(1).array();
  Eigen::ArrayXd isToIa = psi * state_matrix.col(2).array();
  Eigen::ArrayXd iaToR = gamma * state_matrix.col(3).array();

  // compartmental changes
  Eigen::ArrayXd dS = births + rToS - sToE - (state_matrix.col(0).array() * d);
  Eigen::ArrayXd dE = sToE - eToIa - eToIs - (state_matrix.col(1).array() * d);
  Eigen::ArrayXd dIs = eToIs - isToIa - (state_matrix.col(2).array() * d);
  Eigen::ArrayXd dIa =
      eToIa + isToIa + rToIa - iaToR - (state_matrix.col(3).array() * d);
  Eigen::ArrayXd dR = iaToR - rToS - rToIa - (state_matrix.col(4).array() * d);
  // must also return re-infections as rToIa, and new infections as sToE

  // create return array
  Eigen::VectorXd vec_changes(state_matrix.size());
  vec_changes << dS, dE, dIs, dIa, dR, rToIa, sToE;

  return Rcpp::List::create(vec_changes);
}

struct norovirus_model {
  const double rho, b, d, sigma, epsilon, psi, gamma;
  double delta, w1, w2, q1, q2;
  Eigen::MatrixXd contacts;
  Eigen::ArrayXd param_;
  // npi, interv, pop
  explicit norovirus_model(const Rcpp::List params)
      : rho(params["rho"]),
        b(params["b"]),
        d(params["d"]),
        sigma(params["sigma"]),
        epsilon(params["epsilon"]),
        psi(params["psi"]),
        gamma(params["gamma"]),
        delta(params["D_immun"]),
        w1(params["season_amp"]),
        w2(params["season_offset"]),
        q1(params["probT_under5"]),
        q2(params["probT_over5"]),
        contacts(Rcpp::as<Eigen::MatrixXd>(params["contacts"])) {}

  void init_model() {
    // parameters
    delta = 1.0 / (delta * 365.0);
    w1 = w1 / 100.0;
    w2 = w2 / 10.0;
    q1 = q1 / 10.0;
    q2 = q2 / 100.0;

    // param vector
    param_ = Eigen::ArrayXd (4);
    param_ << q1, q2, q2, q2;

  }

  void operator()(const odetools::state_type &state_matrix,
                  odetools::state_type &dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of state_matrix
    dxdt.resize(state_matrix.rows(), state_matrix.cols());

    // prepare seasonal forcing
    const double seasonal_term = seasonal_forcing(t, w1, w2);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations

    // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R
    // calculate new infections
    Eigen::ArrayXd sToE =
      param_ * seasonal_term * state_matrix.col(0).array() *
      (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho))).array();

    // calculate re-infections
    // recovered are column index 4 of the initial conditions
    Eigen::ArrayXd rToIa =
      param_ * seasonal_term * state_matrix.col(4).array() *
      (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho))).array();

    // get current population size
    const Eigen::ArrayXd population_size =
        state_matrix.block(0, 0, 4, 5).rowwise().sum();

    // compartmental transitions
    Eigen::ArrayXd births = (b * population_size);
    Eigen::ArrayXd rToS = (delta * state_matrix.col(4));
    Eigen::ArrayXd eToIa =
        ((1.0 - sigma) * epsilon) * state_matrix.col(1).array();
    Eigen::ArrayXd eToIs = (sigma * epsilon) * state_matrix.col(1).array();
    Eigen::ArrayXd isToIa = psi * state_matrix.col(2).array();
    Eigen::ArrayXd iaToR = gamma * state_matrix.col(3).array();

    // compartmental changes accounting for contacts (for dS and dE)
    dxdt.col(0) = births + rToS - sToE - (state_matrix.col(0).array() * d);
    dxdt.col(1) = sToE - eToIa - eToIs - (state_matrix.col(1).array() * d);
    dxdt.col(2) = eToIs - isToIa - (state_matrix.col(2).array() * d);
    dxdt.col(3) =
        eToIa + isToIa + rToIa - iaToR - (state_matrix.col(3).array() * d);
    dxdt.col(4) = iaToR - rToS - rToIa - (state_matrix.col(4).array() * d);
    dxdt.col(5) = rToIa;
    dxdt.col(6) = sToE;
  }
};

//' @title Run an age-structured SEIRS for norovirus
//'
//' @param initial_conditions An integer matrix holding the initial conditions
//' of the simulation. Each column should represent one compartment, in the
//' order: susceptible, exposed, infectious and symptomatic, infectious and
//' asymptomatic, and recovered. Two extra columns for re-infections and
//' new infections are are required.
//' Rows must represent age groups.
//' @param params A `list` object with infection parameters.
//' @param time_end The maximum time, defaults to 200.0.
//' @param increment The increment time, defaults to 0.1.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix.
//' The second list element is a vector of timesteps.
//' Pass this to the function `epidemics:::output_to_df()` to get a `data.table`
//' of simulation results.
//' @export
// [[Rcpp::export(name="noromod_cpp_boost")]]
Rcpp::List noromod_cpp_boost(
    const Eigen::MatrixXd &initial_conditions, const Rcpp::List &params,
    const double &time_end = 200.0,  // double required by boost solver
    const double &increment = 1.0) {
  // initial conditions from input
  odetools::state_type x = initial_conditions;

  // create a default epidemic with parameters
  norovirus_model this_model(params);
  this_model.init_model();

  // prepare storage containers for the observer
  std::vector<odetools::state_type> x_vec;  // is a vector of MatrixXd
  std::vector<double> times;

  // a controlled stepper for constant step sizes
  boost::numeric::odeint::runge_kutta4<
      odetools::state_type, double, odetools::state_type, double,
      boost::numeric::odeint::vector_space_algebra>
      stepper;

  // run the function without assignment
  boost::numeric::odeint::integrate_const(stepper, this_model, x, 0.0, time_end,
                                          increment,
                                          odetools::observer(x_vec, times));

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x_vec),
                            Rcpp::Named("time") = Rcpp::wrap(times));
}
