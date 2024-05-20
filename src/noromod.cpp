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

typedef Eigen::VectorXd state_type;

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

struct norovirus_model {
  const double rho, b, d;
  double births = 0.0;
  Eigen::Tensor<double, 1> d_rate;  // remove background mortality
  Rcpp::NumericVector sigma;        // has size 3L for 3 levels of protection
  Eigen::Tensor<double, 2> sigma_tensor =
      Eigen::Tensor<double, 2>(3, 3);  // for 3 levels of protection
  const double epsilon, psi, gamma;
  double delta, w1, q1, q2;
  double seasonal_term = 0.0;
  const std::vector<double> w2_values, season_change_points;
  std::vector<double> contacts, aging;  // cannot be const as mapped to Tensor
  Rcpp::NumericVector phi_1, phi_2, upsilon_1, upsilon_2;
  Eigen::Tensor<double, 2> contacts_tensor, aging_tensor;
  Eigen::Tensor<double, 1> param_;
  // npi, interv, pop
  explicit norovirus_model(const Rcpp::List &params)
      : rho(params["rho"]),
        b(params["b"]),
        d(params["d"]),
        sigma(params["sigma"]),
        epsilon(params["epsilon"]),
        psi(params["psi"]),
        gamma(params["gamma"]),
        delta(params["D_immun"]),
        w1(params["season_amp"]),
        q1(params["probT_under5"]),
        q2(params["probT_over5"]),
        phi_1(params["phi_1"]),
        phi_2(params["phi_2"]),
        upsilon_1(params["upsilon_1"]),
        upsilon_2(params["upsilon_2"]),
        w2_values(Rcpp::as<std::vector<double>>(params["season_offset"])),
        season_change_points(
            Rcpp::as<std::vector<double>>(params["season_change_points"])),
        contacts(Rcpp::as<std::vector<double>>(params["contacts"])),
        aging(Rcpp::as<std::vector<double>>(params["aging"])) {}

  void init_model() {
    // parameters
    // w2 now processed in operator
    delta = 1.0 / (delta * 365.0);
    w1 = w1 / 100.0;
    q1 = std::exp(q1);
    q2 = std::exp(q2);

    // param Tensor
    param_.setValues({q1, q2, q2, q2});

    // TODO: vaccination and waning

    // Sigma Tensor initialisation, create a diagonal matrix
    sigma_tensor.setZero();
    for (size_t i = 0; i < sigma.size(); i++) {
      sigma_tensor(i, i) = sigma[i];
    }

    // map to tensors; number of age groups is hardcoded to 4 giving 4x4
    // matrices
    contacts_tensor =
        Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
            &contacts[0], 4, 4);
    aging_tensor = Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>>(
        &aging[0], 4, 4);
  }

  // add nolint flags to allow passing by reference
  void operator()(state_type &x,   // NOLINT
                  state_type &dx,  // NOLINT
                  const double t) {
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
        contacts_tensor.contract(groupwise_infected, product_dims) *
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
    Eigen::Tensor<double, 2> recovery_waning = (delta * x_tensor.chip(4, 0));
    dx_tensor.chip(0, 1) = recovery_waning - new_infections;
    // births enter the unvaccinated susceptible, lowest age group only
    dx_tensor(0, 0, 0) = dx_tensor(0, 0, 0) + births;

    // changes in exposed from infections
    dx_tensor.chip(1, 1) = -(epsilon * x_tensor.chip(1, 1)) + new_infections;

    // changes in infectious asymptomatic
    dx_tensor.chip(2, 1) =
        (-psi * x_tensor.chip(2, 1)) +
        (epsilon * x_tensor.chip(2, 1).contract(sigma_tensor, product_dims));

    // changes in infectious symptomatic
    dx_tensor.chip(3, 1) =
        (psi * x_tensor.chip(2, 1)) - (gamma * x_tensor.chip(3, 1)) +
        (epsilon *
         x_tensor.chip(2, 1).contract(1.0 - sigma_tensor, product_dims));
    dx_tensor.chip(3, 1) = dx_tensor.chip(3, 1) + re_infections;

    // changes in recoveries
    dx_tensor.chip(4, 1) = (gamma * x_tensor.chip(3, 1)) - re_infections;

    // calculate background mortality in all epi compartments; uniform mortality
    // rate
    dx_tensor.slice(offsets, extents) = dx_tensor.slice(offsets, extents) -
                                        (x_tensor.slice(offsets, extents) * d);

    // apply aging-related flows to all compartments and all vaccination strata
    // NOTE: dx = x + aging; aging matrix coeffs have appropriate signs
    for (size_t i = 0; i < 3L; i++) {
      dx_tensor.slice(offsets, extents) =
          dx_tensor.slice(offsets, extents) +
          (aging_tensor.contract(x_tensor.slice(offsets, extents),
                                 product_dims));
    }

    // new infections
    dx_tensor.chip(5, 1) = new_infections;

    // reinfections
    dx_tensor.chip(6, 1) = re_infections;
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
