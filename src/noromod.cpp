// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

// clang-format on

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

typedef std::vector<double> state_type;

/// @brief A simple observer for the integrator
struct observer {
  std::vector<state_type> &m_states;
  std::vector<double> &m_times;

  /// @brief Constructor for the observer
  /// @param states A vector of `state_type`, the model compartments
  /// @param times A vector of doubles, the model times logged
  observer(std::vector<state_type> &states,  // NOLINT
           std::vector<double> &times)       // NOLINT
      : m_states(states), m_times(times) {}

  void operator()(const state_type &x, double t) {
    m_states.push_back(x);
    m_times.push_back(t);
  }
};

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
      Eigen::Map<Eigen::MatrixXd>(state.data(), 4L, 21L);

  // get ageing matrix
  Eigen::MatrixXd aging = Rcpp::as<Eigen::MatrixXd>(parameters["aging"]);

  // get current population size
  const Eigen::ArrayXd population_size =
      state_matrix.block(0, 0, 4, 15).rowwise().sum();

  // modify parameters
  const double delta = 1.0 / (Rcpp::as<double>(parameters["D_immun"]) * 365.0);
  const double w1 = Rcpp::as<double>(parameters["season_amp"]) / 100.0;
  const Rcpp::NumericVector w2_values = parameters["season_offset"];
  // std::exp operates after conversion to STL double
  const double q1 = std::exp(Rcpp::as<double>(parameters["probT_under5"]));
  const double q2 = std::exp(Rcpp::as<double>(parameters["probT_over5"]));
  // create contact matrix
  const Eigen::MatrixXd contacts = parameters["contacts"];

  // epi parameters
  const double rho = parameters["rho"];
  const double b = parameters["b"];
  const Eigen::ArrayXd d = Rcpp::as<Eigen::ArrayXd>(parameters["d"]);
  const double sigma = parameters["sigma"][0];
  const double sigma_v1 = parameters["sigma"][1];
  const double sigma_v2 = parameters["sigma"][2];
  const double epsilon = parameters["epsilon"];
  const double psi = parameters["psi"];
  const double gamma = parameters["gamma"];

  // vaccination rates
  const double phi_1 = parameters["phi_1"];
  const double phi_2 = parameters["phi_2"];
  
  // waning rates
  const double upsilon_1 = 1.0 / (Rcpp::NumericVector(parameters["upsilon"])[0] * 365.0);
  const double upsilon_2 = 1.0 / (Rcpp::NumericVector(parameters["upsilon"])[1] * 365.0);
  
  // get seasonal change points
  const Rcpp::NumericVector season_change_points =
      parameters["season_change_points"];

  // prepare current w2 value
  double w2_current = 0.0;
  // expect that intervals sequences are in order
  for (size_t i = 0; i < season_change_points.size(); i++) {
    if (t <= season_change_points[i]) {
      w2_current = w2_values[i] / 100.0;  // division by 100.0 moved here
      break;                              // exit the loop
    }
  }

  // prepare array for multiplication
  Eigen::ArrayXd param_(4L);
  param_ << q1, q2, q2, q2;

  // prepare seasonal forcing
  const double seasonal_term = seasonal_forcing(t, w1, w2_current);

  // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R, 5:rToIa, 6:sToE, 7:Sv1, 8:Ev1, 9:Isv1, 
  // 10:Iav1, 11:Rv1, 12:rv1ToIav1, 13:sv1ToEv1, 14:Sv2, 15:Ev2, 16:Isv2, 17:Iav2, 
  // 18:Rv2, 19:rv2ToIav2, 20:sv2ToEv2
  // calculate new infections
  Eigen::ArrayXd sToE =
      param_ * seasonal_term * state_matrix.col(0).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

  Eigen::ArrayXd s_v1ToE_v1 =
      param_ * seasonal_term * state_matrix.col(7).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

  Eigen::ArrayXd s_v2ToE_v2 =
      param_ * seasonal_term * state_matrix.col(14).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

  // calculate re-infections
  // recovered are column index 4 of the initial conditions
  Eigen::ArrayXd rToIa =
      param_ * seasonal_term * state_matrix.col(4).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

  Eigen::ArrayXd r_v1ToIa_v1 =
      param_ * seasonal_term * state_matrix.col(11).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();
  
  Eigen::ArrayXd r_v2ToIa_v2 =
      param_ * seasonal_term * state_matrix.col(18).array() *
      (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
      ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();
         
  // compartmental transitions
  double births_ = (b * population_size.sum());
  Eigen::ArrayXd births(4);  // hardcoded for four age groups
  births << births_, 0.0, 0.0, 0.0;

  Eigen::ArrayXd rToS = (delta * state_matrix.col(4));
  Eigen::ArrayXd eToIa =
      ((1.0 - sigma) * epsilon) * state_matrix.col(1).array();
  Eigen::ArrayXd eToIs = (sigma * epsilon) * state_matrix.col(1).array();
  Eigen::ArrayXd isToIa = psi * state_matrix.col(2).array();
  Eigen::ArrayXd iaToR = gamma * state_matrix.col(3).array();
  
  // v1 transitions
  Eigen::ArrayXd r_v1ToS_v1 = (delta * state_matrix.col(11));
  Eigen::ArrayXd e_v1ToIa_v1 =
    ((1.0 - sigma_v1) * epsilon) * state_matrix.col(8).array();
  Eigen::ArrayXd e_v1ToIs_v1 = (sigma_v1 * epsilon) * state_matrix.col(8).array();
  Eigen::ArrayXd is_v1ToIa_v1 = psi * state_matrix.col(9).array();
  Eigen::ArrayXd ia_v1ToR_v1 = gamma * state_matrix.col(10).array();
  
  Eigen::ArrayXd sToS_v1 = phi_1 * state_matrix.col(0).array();
  Eigen::ArrayXd s_v1ToS = upsilon_1 * state_matrix.col(7).array();
  Eigen::ArrayXd rToRS_v1 = phi_1 * state_matrix.col(4).array();
  Eigen::ArrayXd r_v1ToR = upsilon_1 * state_matrix.col(11).array();
  Eigen::ArrayXd r_v1ToS = (delta * upsilon_1) * state_matrix.col(11).array();
  
  //v2 transitions
  Eigen::ArrayXd r_v2ToS_v2 = (delta * state_matrix.col(18));
  Eigen::ArrayXd e_v2ToIa_v2 =
    ((1.0 - sigma_v1) * epsilon) * state_matrix.col(15).array();
  Eigen::ArrayXd e_v2ToIs_v2 = (sigma_v2 * epsilon) * state_matrix.col(15).array();
  Eigen::ArrayXd is_v2ToIa_v2 = psi * state_matrix.col(16).array();
  Eigen::ArrayXd ia_v2ToR_v2 = gamma * state_matrix.col(17).array();
  
  Eigen::ArrayXd r_v2ToS_v2 = (delta * state_matrix.col(18));
  Eigen::ArrayXd s_v1ToS_v2 = phi_2 * state_matrix.col(7).array();
  Eigen::ArrayXd s_v2ToS_v1 = upsilon_2 * state_matrix.col(14).array();
  Eigen::ArrayXd r_v1ToR_v2 = phi_2 * state_matrix.col(11).array();
  Eigen::ArrayXd r_v2ToR_v1 = upsilon_2 * state_matrix.col(18).array();
  Eigen::ArrayXd r_v2ToS_v1 = (delta * upsilon_2) * state_matrix.col(18).array();

  // compartmental changes
  Eigen::ArrayXd dS = births + rToS - sToE - sToSv1 + sv1ToS + rv1ToS - (state_matrix.col(0).array() * d) +
                      (aging * state_matrix.col(0)).array();
  Eigen::ArrayXd dE = sToE - eToIa - eToIs - (state_matrix.col(1).array() * d) +
                      (aging * state_matrix.col(1)).array();
  Eigen::ArrayXd dIs = eToIs - isToIa - (state_matrix.col(2).array() * d) +
                       (aging * state_matrix.col(2)).array();
  Eigen::ArrayXd dIa = eToIa + isToIa + rToIa - iaToR -
                       (state_matrix.col(3).array() * d) +
                       (aging * state_matrix.col(3)).array();
  Eigen::ArrayXd dR = iaToR - rToS - rToIa + rv1ToR - (state_matrix.col(4).array() * d) +
                      (aging * state_matrix.col(4)).array();
  
  Eigen::ArrayXd dS_v1 = r_v1ToSv1 - s_v1ToE_v1 - sToS_v1 + s_v1ToS - s_v1ToS_v2 + s_v2ToS_v1 -
                      (state_matrix.col(7).array() * d) + (aging * state_matrix.col(7)).array();
  Eigen::ArrayXd dE_v1 = s_v1ToEv1 - e_v1ToIa_v1 - e_v1ToIs_v1 - (state_matrix.col(8).array() * d);
  Eigen::ArrayXd dIs_v1 = e_v1ToIs_v1 - is_v1ToIa_v1 - (state_matrix.col(9).array() * d) +
                       (aging * state_matrix.col(9)).array();
  Eigen::ArrayXd dIa_v1 = e_v1ToIa_v1 + is_v1ToIa_v1 + r_v1ToIa_v1 - ia_v1ToR_v1 - 
                       (state_matrix.col(10).array() * d) + (aging * state_matrix.col(10)).array();
  Eigen::ArrayXd dR_v1 = ia_v1ToR_v1 - r_v1ToS_v1 - r_v1ToIa_v1 - r_v1ToR - r_v1ToS  + rToR_v1 + r_v2ToRv1 - r_v1ToR_v2
                      - (state_matrix.col(11).array() * d) + (aging * state_matrix.col(11)).array();)
    
  Eigen::ArrayXd dS_v2 = r_v2ToSv2 - s_v2ToE_v2 + s_v1ToS_v2 - s_v2ToS_v1 - (state_matrix.col(14).array() * d) +
                      (aging * state_matrix.col(14)).array();
  Eigen::ArrayXd dE_v2 = s_v2ToEv2 - e_v2ToIa_v2 - e_v2ToIs_v2 - (state_matrix.col(15).array() * d);
  Eigen::ArrayXd dIs_v2 = e_v2ToIs_v2 - is_v2ToIa_v2 - (state_matrix.col(16).array() * d) +
                       (aging * state_matrix.col(16)).array();
  Eigen::ArrayXd dIa_v2 = e_v2ToIa_v2 + is_v2ToIa_v2 + r_v2ToIa_v2 - ia_v2ToR_v2 - 
                       (state_matrix.col(17).array() * d) + (aging * state_matrix.col(17)).array();
  Eigen::ArrayXd dR_v2 = ia_v2ToR_v2 - r_v2ToS_v2 - r_v2ToIa_v2 - r_v2ToS_v1  + r_v1ToRv2 - r_v2ToR_v1
                      - (state_matrix.col(18).array() * d) + (aging * state_matrix.col(18)).array();)
  
  // must also return re-infections as rToIa, and new infections as sToE

  // create return array
  Eigen::VectorXd vec_changes(state_matrix.size());
  vec_changes << dS, dE, dIs, dIa, dR, rToIa, sToE, dS_v1, dE_v1, dIs_v1, dIa_v1, dR_v1, r_v1ToIa_v1, s_v1ToE_v1, dS_v2, dE_v2, dIs_v2, dIa_v2, dR_v2, r_v2ToIa_v2, s_v2ToE_v2;;

  return Rcpp::List::create(vec_changes);
}

struct norovirus_model {
  const double rho, b;
  Eigen::ArrayXd d;
  const double sigma, epsilon, psi, gamma, phi, upsilon;
  double delta, w1, q1, q2;
  const std::vector<double> w2_values, season_change_points;
  Eigen::MatrixXd contacts, aging;
  Eigen::ArrayXd param_;
  // npi, interv, pop
  explicit norovirus_model(const Rcpp::List params)
      : rho(params["rho"]),
        b(params["b"]),
        d(Rcpp::as<Eigen::ArrayXd>(params["d"])),
        sigma(params["sigma"]),
        epsilon(params["epsilon"]),
        psi(params["psi"]),
        gamma(params["gamma"]),
        delta(params["D_immun"]),
        w1(params["season_amp"]),
        q1(params["probT_under5"]),
        q2(params["probT_over5"]),
        phi(params["phi"]),
        upsilon(params["upsilon"]),
        q2(params["probT_over5"]),
        w2_values(Rcpp::as<std::vector<double> >(params["season_offset"])),
        season_change_points(
            Rcpp::as<std::vector<double> >(params["season_change_points"])),
        contacts(Rcpp::as<Eigen::MatrixXd>(params["contacts"])),
        aging(Rcpp::as<Eigen::MatrixXd>(params["aging"])) {}

  void init_model() {
    // parameters
    // w2 now processed in operator
    delta = 1.0 / (delta * 365.0);
    w1 = w1 / 100.0;
    q1 = std::exp(q1);
    q2 = std::exp(q2);

    // param vector
    param_ = Eigen::ArrayXd(4);
    param_ << q1, q2, q2, q2;
  }

  void operator()(const odetools::state_type &state_matrix,
                  odetools::state_type &dxdt,  // NOLINT
                  const double t) {
    // resize the dxdt vector to the dimensions of state_matrix
    dxdt.resize(state_matrix.rows(), state_matrix.cols());

    // prepare w2_current, initially first value
    double w2_current = w2_values[0] / 100.0;
    // expect that intervals sequences are in order
    for (size_t i = 0; i < season_change_points.size(); i++) {
      if (t <= season_change_points[i]) {
        w2_current = w2_values[i] / 100.0;  // division by 100.0 now here
        break;                              // exit the loop
      }
    }

    // prepare seasonal forcing
    const double seasonal_term = seasonal_forcing(t, w1, w2_current);

    // NB: Casting initial conditions matrix columns to arrays is necessary
    // for vectorised operations

    // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R
    // calculate new infections
    Eigen::ArrayXd sToE =
        param_ * seasonal_term * state_matrix.col(0).array() *
        (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho)))
            .array();
        
    Eigen::ArrayXd s_v1ToE_v1 =
        param_ * seasonal_term * state_matrix.col(7).array() *
        (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
        ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();
           
    Eigen::ArrayXd s_v2ToE_v2 =
        param_ * seasonal_term * state_matrix.col(14).array() *
        (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
        ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

    // calculate re-infections
    // recovered are column index 4 of the initial conditions
    Eigen::ArrayXd rToIa =
        param_ * seasonal_term * state_matrix.col(4).array() *
        (contacts * (state_matrix.col(2) + (state_matrix.col(3) * rho)))
            .array();
        
    Eigen::ArrayXd r_v1ToIa_v1 =
        param_ * seasonal_term * state_matrix.col(11).array() *
        (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
        ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();
           
    Eigen::ArrayXd r_v2ToIa_v2 =
         param_ * seasonal_term * state_matrix.col(18).array() *
         (contacts * ((state_matrix.col(2) + state_matrix.col(9) + state_matrix.col(16)) +
         ((state_matrix.col(3) + state_matrix.col(10) + state_matrix.col(17) * rho))).array();

    // get current population size
    const Eigen::ArrayXd population_size =
        state_matrix.block(0, 0, 4, 5).rowwise().sum();

    // compartmental transitions
    double births_ = (b * population_size.sum());
    Eigen::ArrayXd births(4);  // hardcoded for four age groups
    births << births_, 0.0, 0.0, 0.0;

    Eigen::ArrayXd dS = births + rToS - sToE - sToSv1 + sv1ToS + rv1ToS - (state_matrix.col(0).array() * d) +
      (aging * state_matrix.col(0)).array();
    Eigen::ArrayXd dE = sToE - eToIa - eToIs - (state_matrix.col(1).array() * d) +
      (aging * state_matrix.col(1)).array();
    Eigen::ArrayXd dIs = eToIs - isToIa - (state_matrix.col(2).array() * d) +
      (aging * state_matrix.col(2)).array();
    Eigen::ArrayXd dIa = eToIa + isToIa + rToIa - iaToR -
      (state_matrix.col(3).array() * d) +
      (aging * state_matrix.col(3)).array();
    Eigen::ArrayXd dR = iaToR - rToS - rToIa + rv1ToR - (state_matrix.col(4).array() * d) +
      (aging * state_matrix.col(4)).array();
    
    Eigen::ArrayXd dS_v1 = r_v1ToSv1 - s_v1ToE_v1 - sToS_v1 + s_v1ToS - s_v1ToS_v2 + s_v2ToS_v1 -
      (state_matrix.col(7).array() * d) + (aging * state_matrix.col(7)).array();
    Eigen::ArrayXd dE_v1 = s_v1ToEv1 - e_v1ToIa_v1 - e_v1ToIs_v1 - (state_matrix.col(8).array() * d);
    Eigen::ArrayXd dIs_v1 = e_v1ToIs_v1 - is_v1ToIa_v1 - (state_matrix.col(9).array() * d) +
      (aging * state_matrix.col(9)).array();
    Eigen::ArrayXd dIa_v1 = e_v1ToIa_v1 + is_v1ToIa_v1 + r_v1ToIa_v1 - ia_v1ToR_v1 - 
      (state_matrix.col(10).array() * d) + (aging * state_matrix.col(10)).array();
    Eigen::ArrayXd dR_v1 = ia_v1ToR_v1 - r_v1ToS_v1 - r_v1ToIa_v1 - r_v1ToR - r_v1ToS  + rToR_v1 + r_v2ToRv1 - r_v1ToR_v2
    - (state_matrix.col(11).array() * d) + (aging * state_matrix.col(11)).array();)
      
      Eigen::ArrayXd dS_v2 = r_v2ToSv2 - s_v2ToE_v2 + s_v1ToS_v2 - s_v2ToS_v1 - (state_matrix.col(14).array() * d) +
      (aging * state_matrix.col(14)).array();
    Eigen::ArrayXd dE_v2 = s_v2ToEv2 - e_v2ToIa_v2 - e_v2ToIs_v2 - (state_matrix.col(15).array() * d);
    Eigen::ArrayXd dIs_v2 = e_v2ToIs_v2 - is_v2ToIa_v2 - (state_matrix.col(16).array() * d) +
      (aging * state_matrix.col(16)).array();
    Eigen::ArrayXd dIa_v2 = e_v2ToIa_v2 + is_v2ToIa_v2 + r_v2ToIa_v2 - ia_v2ToR_v2 - 
      (state_matrix.col(17).array() * d) + (aging * state_matrix.col(17)).array();
    Eigen::ArrayXd dR_v2 = ia_v2ToR_v2 - r_v2ToS_v2 - r_v2ToIa_v2 - r_v2ToS_v1  + r_v1ToRv2 - r_v2ToR_v1
    - (state_matrix.col(18).array() * d) + (aging * state_matrix.col(18)).array();)

    // compartmental changes accounting for contacts (for dS and dE)
    dxdt.col(0) = births + rToS - sToE - (state_matrix.col(0).array() * d) +
                  (aging * state_matrix.col(0)).array();
    dxdt.col(1) = sToE - eToIa - eToIs - (state_matrix.col(1).array() * d) +
                  (aging * state_matrix.col(1)).array();
    dxdt.col(2) = eToIs - isToIa - (state_matrix.col(2).array() * d) +
                  (aging * state_matrix.col(2)).array();
    dxdt.col(3) = eToIa + isToIa + rToIa - iaToR -
                  (state_matrix.col(3).array() * d) +
                  (aging * state_matrix.col(3)).array();
    dxdt.col(4) = iaToR - rToS - rToIa - (state_matrix.col(4).array() * d) +
                  (aging * state_matrix.col(4)).array();
    dxdt.col(5) = rToIa;
    dxdt.col(6) = sToE;
    dxdt.col(7) = r_v1ToSv1 - s_v1ToE_v1 - sToS_v1 + s_v1ToS - s_v1ToS_v2 + s_v2ToS_v1 -
      (state_matrix.col(7).array() * d) + (aging * state_matrix.col(7)).array();
    dxdt.col(8) = s_v1ToEv1 - e_v1ToIa_v1 - e_v1ToIs_v1 - (state_matrix.col(8).array() * d);
    dxdt.col(9) = e_v1ToIs_v1 - is_v1ToIa_v1 - (state_matrix.col(9).array() * d) +
      (aging * state_matrix.col(9)).array();
    dxdt.col(10) = e_v1ToIa_v1 + is_v1ToIa_v1 + r_v1ToIa_v1 - ia_v1ToR_v1 - 
      (state_matrix.col(10).array() * d) + (aging * state_matrix.col(10)).array();
    dx.dt(11) = ia_v1ToR_v1 - r_v1ToS_v1 - r_v1ToIa_v1 - r_v1ToR - r_v1ToS  + rToR_v1 + r_v2ToRv1 - r_v1ToR_v2
    - (state_matrix.col(11).array() * d) + (aging * state_matrix.col(11)).array();)
    dx.dt(12) = r_v1ToIa_v1;
    dx.dt(13) = s_v1ToE_v1;
    dx.dt(14) = r_v2ToSv2 - s_v2ToE_v2 + s_v1ToS_v2 - s_v2ToS_v1 - (state_matrix.col(14).array() * d) +
      (aging * state_matrix.col(14)).array();
    dx.dt(15) = s_v2ToEv2 - e_v2ToIa_v2 - e_v2ToIs_v2 - (state_matrix.col(15).array() * d);
    dx.dt(16) = e_v2ToIs_v2 - is_v2ToIa_v2 - (state_matrix.col(16).array() * d) +
      (aging * state_matrix.col(16)).array();
    dx.dt(17) = e_v2ToIa_v2 + is_v2ToIa_v2 + r_v2ToIa_v2 - ia_v2ToR_v2 - 
      (state_matrix.col(17).array() * d) + (aging * state_matrix.col(17)).array();
    dx.dt(18) = ia_v2ToR_v2 - r_v2ToS_v2 - r_v2ToIa_v2 - r_v2ToS_v1  + r_v1ToRv2 - r_v2ToR_v1 - (state_matrix.col(18).array() * d) + (aging * state_matrix.col(18)).array();)
    dx.dt(19) = r_v2ToIa_v2;
    dx.dt(20) = s_v2ToE_v2;
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
