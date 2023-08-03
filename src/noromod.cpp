// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

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
