// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

/// @brief
/// @param t
/// @param w1
/// @param w2
/// @return
const double seasonal_forcing(const double &t, const double &w1,
                              const double &w2) {
  const double z = 1.0 + (w1 * std::cos((2.0 * M_PI * t / 364.0) + w2));
  return z;
}

//' @title
//' @description
//' @param
//' @return
//' @export
// [[Rcpp::export]]
Rcpp::List norovirus_model(const double &t, const Eigen::MatrixXd &state,
                           Rcpp::NumericVector parameters) {
  // modify parameters
  const double delta = 1.0 / (parameters["D_immun"] * 365.0);
  const double w1 = parameters["season_amp"] / 100.0;
  const double w2 = parameters["season_offset"] / 10.0;
  const double q1 = parameters["probT_under5"] / 10.0;
  const double q2 = parameters["probT_over5"] / 100.0;
  const double contacts = parameters["contacts"];
  
  // epi parameters
  const double rho = parameters["rho"];
  const double beta = parameters["beta"];
  const double d = parameters["d"];

  // prepare array for multiplication
  Eigen::ArrayXd param_ << q1, q2, q2, q2;

  // prepare seasonal forcing
  const double seasonal_term = seasonal_forcing(t, parameters["season_amp"],
                                                parameters["season_offset"]);

  // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R
  // calculate new infections
  Eigen::ArrayXd sToE =
      state.col(0).array() *
      (contacts * (state.col(2).array() + state.col(3).array() * rho)).array();
  // multiply parameters and seasonal effect
  sToE = param_ * seasonal_term * new_infections();

  // calculate re-infections
  // recovered are column index 4 of the initial conditions
  Eigen::ArrayXd re_infections =
      state.col(4).array() *
      (contacts * (state.col(2).array() + state.col(3).array() * rho)).array();
  // multiply parameters and seasonal effect
  re_infections = param_ * seasonal_term * re_infections();

  // compartmental transitions
  Eigen::ArrayXd rToS = (b * N) + (delta * state.col(4));
  Eigen::ArrayXd sToE = 

  // compartmental changes
  Eigen::ArrayXd dS = 
  dxdt.col(0) = -sToE - sToV;  // -β*S*contacts*I - ν*S
  dxdt.col(1) = sToE - eToI;   // β*S*contacts*I - α*E
  dxdt.col(2) = eToI - iToR;   // α*E - γ*I
  dxdt.col(3) = iToR;          // γ*I
  dxdt.col(4) = sToV;          // ν*S
}
