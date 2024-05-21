// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>

#include <unsupported/Eigen/CXX11/Tensor>
#include <cmath>
#include <vector>
#include <boost/numeric/odeint.hpp>
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

  /// @brief Overloaded operator for the observer structure
  /// @param x The current system state x.
  /// @param t The current system time t.
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

/// @brief A structure for the RHS of the ODE system. Must have an overloaded
/// operator function to be a `FunctionObject`, allowing it to be passed to an
/// ODE integrator from `boost::odeint`.
struct norovirus_model {
  const double rho, b, d;  // single value for d (background mortality)
  double births = 0.0;
  Rcpp::NumericVector sigma_vec;  // has size 3L for 3 levels of protection
  Eigen::Tensor<double, 2> sigma =
      Eigen::Tensor<double, 2>(3, 3);  // for 3 levels of protection
  const double epsilon, psi, gamma;
  double delta, w1, q1, q2;
  double seasonal_term = 0.0;
  const std::vector<double> w2_values, season_change_points;

  std::vector<double> contacts_vec, aging_vec;
  // vaccination rate and vaccine-immunity waning rate
  std::vector<double> phi_vec, upsilon_vec;

  Eigen::Tensor<double, 2> contacts, aging = Eigen::Tensor<double, 2>(4, 4);
  Eigen::Tensor<double, 2> phi, upsilon = Eigen::Tensor<double, 2>(4, 3);
  Eigen::Tensor<double, 1> param_ = Eigen::Tensor<double, 1>(4);
  // npi, interv, pop

  /// @brief Constructor list for the norovirus model structure.
  /// @param params A named list of model parameters; this is taken directly
  /// from the parameter list passed from R.
  explicit norovirus_model(const Rcpp::List &params)
      : rho(params["rho"]),
        b(params["b"]),
        d(params["d"]),
        sigma_vec(params["sigma"]),
        epsilon(params["epsilon"]),
        psi(params["psi"]),
        gamma(params["gamma"]),
        delta(params["D_immun"]),
        w1(params["season_amp"]),
        q1(params["probT_under5"]),
        q2(params["probT_over5"]),
        phi_vec(Rcpp::as<std::vector<double>>(params["phi"])),
        upsilon_vec(Rcpp::as<std::vector<double>>(params["upsilon"])),
        w2_values(Rcpp::as<std::vector<double>>(params["season_offset"])),
        season_change_points(
            Rcpp::as<std::vector<double>>(params["season_change_points"])),
        contacts_vec(Rcpp::as<std::vector<double>>(params["contacts"])),
        aging_vec(Rcpp::as<std::vector<double>>(params["aging"])) {}

  /// @brief A helper function for some parameter value transformations that are
  /// not allowed in the constructor.
  void init_model() {
    // parameters
    // w2 now processed in operator
    delta = 1.0 / (delta * 365.0);
    w1 = w1 / 100.0;
    q1 = std::exp(q1);
    q2 = std::exp(q2);

    // param Tensor
    param_.setValues({q1, q2, q2, q2});

    // Sigma Tensor initialisation, create a diagonal matrix
    sigma.setZero();
    for (size_t i = 0; i < sigma_vec.size(); i++) {
      sigma(i, i) = sigma_vec[i];
    }

    // map to tensors; number of age groups is hardcoded to 4 giving 4x4
    // matrices
    contacts = Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
        &contacts_vec[0], 4, 4);
    aging = Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
        &aging_vec[0], 4, 4);

    // NOTE: 4 x 3 2D tensors for 4 ages and 3 vaccine strata
    phi = Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
        &phi_vec[0], 4, 3);
    upsilon = Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
        &upsilon_vec[0], 4, 3);
  }

  // add nolint flags to prevent linting of x and dx as non-const references
  /// @brief An overloaded operator function that represents the RHS of the ODE
  /// system.
  /// @param x The system state; `state_type` is an `Eigen::VectorXd`. Passed as
  /// a non-const reference so that it can be mapped to a 3D Eigen::Tensor.
  /// @param dx The system changes; `state_type` is an `Eigen::VectorXd`. Passed
  /// as a non-const reference so that it can be mapped to a 3D Eigen::Tensor.
  /// Eigen::Tensor.
  /// @param t The system time t; required by the `boost::odeint` integrator.
  void operator()(state_type &x,   // NOLINT
                  state_type &dx,  // NOLINT
                  const double &t) {
    // map a tensor to the state vector with required dims (4, 7, 3)
    // for age groups, epi compartments, vaccination strata
    auto x_tensor = Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(
        &x[0], 4, 7, 3);
    auto dx_tensor =
        Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(&dx[0], 4,
                                                                    7, 3);
    dx_tensor.setZero();

    // prepare w2_current, initially first value
    double w2_current = w2_values[0] / 100.0;
    // expect that intervals sequences are in order
    for (size_t i = 0; i < season_change_points.size(); i++) {
      if (t <= season_change_points[i]) {
        w2_current = w2_values[i] / 100.0;  // division by 100.0 now here
        break;                              // exit the loop
      }
    }

    // recalculate seasonal forcing term
    seasonal_term = seasonal_forcing(t, w1, w2_current);

    // prepare matrix product dimensions
    std::array<Eigen::IndexPair<int>, 1> product_dims = {
        Eigen::IndexPair<int>(1, 0)};

    // column indices: 0:S, 1:E, 2:Is, 3:Ia, 4:R, 5:new_infections, 6:re_infects
    // calculate infection potential
    // a 1D tensor of group-wise infections in all strata
    // reshape to two dimensions: 4 x 1
    auto groupwise_infected =
        (x_tensor.chip(2, 1) + (x_tensor.chip(3, 1) * rho))  // dims: 4 x 3
            .sum(std::array<int, 1>({1}))                    // dims: 4
            .reshape(std::array<int, 2>{4, 1});              // dims: 4 x 1

    // matrix multiply contacts x groupwise infected; result is 4 x 1
    // multiply with season coeff
    auto infection_potential =
        contacts.contract(groupwise_infected, product_dims) *
        seasonal_term;  // dims: 4 x 1

    // calculate new infections and reinfections in each vaccination stratum
    // dims are 4 * 3 for both
    auto new_infections = x_tensor.chip(0, 1) * infection_potential.broadcast(
                                                    std::array<int, 2>{1, 3});
    auto re_infections = x_tensor.chip(4, 1) * infection_potential.broadcast(
                                                   std::array<int, 2>{1, 3});

    // calculate births as b * total population; b = mean per-capita births
    // define starting offset, i.e., starting at 1st row, col, stratum
    std::array<Eigen::Index, 3> offsets = {0, 0, 0};
    // all four ages, first five epi compartments, all vax strata
    std::array<Eigen::Index, 3> extents = {4, 5, 3};
    // births is a double for the total number of births, while `b` is the rate
    Eigen::Tensor<double, 0> births_temp =
        x_tensor.slice(offsets, extents).sum() * b;
    births = births_temp.coeff();

    // changes in susceptibles from births and infections
    dx_tensor.chip(0, 1) = (delta * x_tensor.chip(4, 1)) - new_infections;
    // births enter the unvaccinated susceptible, lowest age group only
    dx_tensor(0, 0, 0) = dx_tensor(0, 0, 0) + births;

    // vaccination flows to and from susceptibles
    // outflows from strata due to vaccination or waning
    dx_tensor.chip(0, 1) =
        dx_tensor.chip(0, 1) - (x_tensor.chip(0, 1) * (phi + upsilon));
    // inflows due to waning and vaccination
    // Inflow into S from waning of V1 and RV1
    dx_tensor.chip(0, 1).chip(0, 1) =
        dx_tensor.chip(0, 1).chip(0, 1) +
        (x_tensor.chip(0, 1).chip(1, 1) * upsilon.chip(1, 1)) +
        (x_tensor.chip(4, 1).chip(1, 1) * delta * gamma);  // direct from RV1
    // Inflow into V1 from waning and vaccination
    dx_tensor.chip(0, 1).chip(1, 1) =
        dx_tensor.chip(0, 1).chip(1, 1) +
        // waning of V2
        (x_tensor.chip(0, 1).chip(2, 1) * upsilon.chip(2, 1)) +
        // vaccination of S
        (x_tensor.chip(0, 1).chip(0, 1) * phi.chip(0, 1)) +
        // direct waning of RV2
        (x_tensor.chip(4, 1).chip(2, 1) * delta * gamma);
    // Inflow into V2 from vaccination
    dx_tensor.chip(0, 1).chip(2, 1) =
        dx_tensor.chip(0, 1).chip(2, 1) +
        (x_tensor.chip(0, 1).chip(1, 1) * phi.chip(1, 1));

    // changes in exposed
    dx_tensor.chip(1, 1) = -(epsilon * x_tensor.chip(1, 1)) + new_infections;

    // changes in infectious asymptomatic
    dx_tensor.chip(2, 1) =
        (-psi * x_tensor.chip(2, 1)) +
        (epsilon * x_tensor.chip(2, 1).contract(sigma, product_dims));

    // changes in infectious symptomatic
    dx_tensor.chip(3, 1) =
        (psi * x_tensor.chip(2, 1)) - (gamma * x_tensor.chip(3, 1)) +
        (epsilon * x_tensor.chip(2, 1).contract(1.0 - sigma, product_dims));
    dx_tensor.chip(3, 1) = dx_tensor.chip(3, 1) + re_infections;

    // changes in recoveries
    dx_tensor.chip(4, 1) = (gamma * x_tensor.chip(3, 1)) - re_infections;

    // vaccination flows to and from recovered
    // outflows from strata due to vaccination or waning
    dx_tensor.chip(4, 1) =
        dx_tensor.chip(4, 1) - (x_tensor.chip(4, 1) * (phi + upsilon));
    // inflows due to waning and vaccination
    // Inflow into R from waning of RV1
    dx_tensor.chip(4, 1).chip(0, 1) =
        dx_tensor.chip(4, 1).chip(0, 1) +
        x_tensor.chip(4, 1).chip(1, 1) * upsilon.chip(1, 1);
    // Inflow into RV1 from waning and vaccination, direct outflow to S
    dx_tensor.chip(4, 1).chip(1, 1) =
        dx_tensor.chip(4, 1).chip(1, 1) +
        // waning of RV2
        (x_tensor.chip(4, 1).chip(2, 1) * upsilon.chip(2, 1)) +
        // vaccination of R
        (x_tensor.chip(4, 1).chip(0, 1) * phi.chip(0, 1)) -
        // direct out to S
        (x_tensor.chip(4, 1).chip(1, 1) * delta * gamma);
    // Inflow into RV2 from vaccination, direct outflow to V1
    dx_tensor.chip(4, 1).chip(2, 1) =
        dx_tensor.chip(4, 1).chip(2, 1) +
        (x_tensor.chip(4, 1).chip(1, 1) * phi.chip(1, 1)) -
        (x_tensor.chip(4, 1).chip(2, 1) * delta * gamma);

    // calculate background mortality in all epi compartments; uniform mortality
    // rate
    dx_tensor.slice(offsets, extents) = dx_tensor.slice(offsets, extents) -
                                        (x_tensor.slice(offsets, extents) * d);

    // apply aging-related flows to all compartments and all vaccination strata
    // NOTE: dx = x + aging; aging matrix coeffs have appropriate signs
    for (size_t i = 0; i < 3L; i++) {
      dx_tensor.slice(offsets, extents) =
          dx_tensor.slice(offsets, extents) +
          (aging.contract(x_tensor.slice(offsets, extents), product_dims));
    }

    // new infections
    dx_tensor.chip(5, 1) = new_infections;

    // reinfections
    dx_tensor.chip(6, 1) = re_infections;
  }
};

//' @title Run an age-structured SEIRSV model for norovirus
//'
//' @param initial_conditions A numeric vector holding the initial conditions
//' of the simulation. Each column should represent one compartment, in the
//' order: susceptible, exposed, infectious and symptomatic, infectious and
//' asymptomatic, and recovered. Two extra columns for re-infections and
//' new infections are also required.
//' Further columns must represent the same compartments for each vaccination
//' stratum.
//' Rows must represent age groups.
//' @param params A `list` object with infection parameters; see the convenience
//' function [default_parameters()].
//' @param time_end The maximum time, defaults to 200.0.
//' @param increment The increment time, defaults to 1.0.
//' @return A two element list, where the first element is a list of matrices
//' whose elements correspond to the numbers of individuals in each compartment
//' as specified in the initial conditions matrix.
//' The second list element is a vector of timesteps.
//' @export
// [[Rcpp::export(name="noromod_cpp_boost")]]
Rcpp::List noromod_cpp_boost(
    Eigen::VectorXd &initial_conditions,  // NOLINT
    const Rcpp::List &params,
    const double &time_end = 200.0,  // double required by boost solver
    const double &increment = 1.0) {
  // initial conditions from input
  state_type x = initial_conditions;

  // map to an Eigen Tensor
  auto state3d = Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::ColMajor>>(
      &initial_conditions[0], 2, 2, 3);

  // create a default epidemic with parameters
  norovirus_model this_model(params);
  this_model.init_model();

  // prepare storage containers for the observer
  std::vector<state_type> x_vec;  // is a vector of vectors
  std::vector<double> times;

  // a controlled stepper for constant step sizes
  boost::numeric::odeint::runge_kutta4<
      state_type, double, state_type, double,
      boost::numeric::odeint::vector_space_algebra>
      stepper;

  // run the function without assignment
  boost::numeric::odeint::integrate_const(stepper, this_model, x, 0.0, time_end,
                                          increment, observer(x_vec, times));

  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x_vec),
                            Rcpp::Named("time") = Rcpp::wrap(times));
}
