// Copyright 2023 'epidemics' authors. See repository licence in LICENSE.md.

// clang-format off
#include <Rcpp.h>
#include <RcppEigen.h>
#include <epidemics.h>

#include <cmath>
#include <boost/numeric/odeint.hpp>
// clang-format on

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

/// @brief The ODE system as a `struct`, treated as a FunctionObject.
struct norovirus_model {
  // NOTE: discontinuous values as new infections and re-infections are counted
  const std::array<double, 15> epi_indices = {0,  1,  2,  3,  4,  7,  8, 9,
                                              10, 11, 14, 15, 16, 17, 18};
  const double rho, b;
  Eigen::ArrayXd d, phi_1, phi_2;  // phi is age-specific
  Rcpp::NumericVector upsilon, sigma;
  double sigma_v0 = 0.0, sigma_v1 = sigma_v0, sigma_v2 = sigma_v1;
  double upsilon_1 = 0.0, upsilon_2 = 0.0;
  const double epsilon, psi, gamma;
  double delta, w1, w3, q1, q2;
  const std::vector<double> w2_values, season_change_points;
  Eigen::MatrixXd contacts, aging;
  Eigen::ArrayXd param_;

  /// @brief The model FunctionObject constructor.
  /// @param params A list of parameters passed from R. See
  /// [default_parameters()] in R/.
  explicit norovirus_model(const Rcpp::List params)
      : rho(params["rho"]),
        b(params["b"]),
        d(Rcpp::as<Eigen::ArrayXd>(params["d"])),
        phi_1(Rcpp::as<Eigen::ArrayXd>(params["phi_1"])),
        phi_2(Rcpp::as<Eigen::ArrayXd>(params["phi_2"])),
        upsilon(params["upsilon"]),
        sigma(params["sigma"]),
        epsilon(params["epsilon"]),
        psi(params["psi"]),
        gamma(params["gamma"]),
        delta(params["D_immun"]),
        w1(params["season_amp"]),
        w3(params["season_amp_over65"]),
        q1(params["probT_under5"]),
        q2(params["probT_over5"]),
        w2_values(Rcpp::as<std::vector<double> >(params["season_offset"])),
        season_change_points(
            Rcpp::as<std::vector<double> >(params["season_change_points"])),
        contacts(Rcpp::as<Eigen::MatrixXd>(params["contacts"])),
        aging(Rcpp::as<Eigen::MatrixXd>(params["aging"])) {}

  /// @brief Initialise some model parameter values. These operations cannot be
  /// accommodated in the constructor.
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

    // prepare sigma param of prop. infectious reduction
    sigma_v0 = sigma[0];
    sigma_v1 = sigma[1];
    sigma_v2 = sigma[2];

    // prepare upsilon
    upsilon_1 = upsilon[0];
    upsilon_2 = upsilon[1];
  }

  /// @brief Overloaded operator function representing the RHS of an ODE system.
  /// The arguments are in the format suitable for
  /// `boost::odeint::integrate_*()`.
  /// @param x The initial state at any time `t`. The `state_type` is an Eigen
  /// Matrix.
  /// @param dxdt The change between time `t` and `t + increment`. Also an Eigen
  /// Matrix. Cannot be `const` as it needs to be modified.
  /// @param t The current time.
  void operator()(const odetools::state_type &x,
                  odetools::state_type &dxdt,  // NOLINT
                  const double &t) {
    // resize the dxdt vector to the dimensions of x
    dxdt.resize(x.rows(), x.cols());
    dxdt.fill(0.0);

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
    // column indices: unvaccinated: 0:S, 1:E, 2:Is, 3:Ia, 4:R
    // column indices: vax 1 dose: 7:S, 8:E, 9:Is, 10:Ia, 11:R
    // column indices: vax 2 dose: 14:S, 15:E, 16:Is, 17:Ia, 18:R
    Eigen::ArrayXd transmission_rates = param_ * seasonal_term;
    // subset infectious as Is + rho*Ia and calculate infection potential
    // return the age-wise sum using `colwise().sum()`
    Eigen::VectorXd infected =
        (x(Eigen::all, std::array<long, 3>{2, 9, 16}) +
         (x(Eigen::all, std::array<long, 3>{3, 10, 17}) * rho))
            .rowwise()
            .sum();
         
         
    // Adjust the 4th age group (index 3) by multiplying it by w3
   infected[3] *= w3;
         
    // return an array of infection potential
    auto infection_potential =
        transmission_rates * (contacts * infected).array();

    // calculate new infections in all vaccine strata
    auto new_infections = infection_potential * x.col(0).array();
    auto new_infections_v1 = infection_potential * x.col(7).array();
    auto new_infections_v2 = infection_potential * x.col(14).array();

    // calculate re-infections in all vaccine strata
    auto reinfections = infection_potential * x.col(4).array();
    auto reinfections_v1 = infection_potential * x.col(11).array();
    auto reinfections_v2 = infection_potential * x.col(18).array();

    // get current population size as a double and calculate new births
    const double population_size = x(Eigen::all, epi_indices).sum();
    double births_ = (b * population_size);
    Eigen::ArrayXd births(4);  // hardcoded for four age groups
    births << births_, 0.0, 0.0, 0.0;

    // compartmental transitions, refer to R/noromod_r.R
    // NOTE: see branch `array_vax` for an implementation using Eigen::Tensor
    // change in susceptibles
    dxdt.col(0) = births - new_infections + (aging * x.col(0)).array() +
                  (delta * x.col(4).array()) -                // R -> S
                  (phi_1 * x.col(0).array()) +                // S -> SV1
                  (upsilon_1 * x.col(7).array()) +            // SV1 -> S
                  (x.col(11).array() * (delta * upsilon_1));  // RV1 -> S
    // change in exposed
    dxdt.col(1) = new_infections + (aging * x.col(1)).array() -
                  (epsilon * x.col(1).array());  // E -> Is + Ia
    // change in infectious symptomatic
    dxdt.col(2) = (epsilon * sigma_v0 * x.col(1).array()) -  // E -> Is
                  (psi * x.col(2).array()) +                 // Is -> Ia
                  (aging * x.col(2)).array();
    // change in infectious asymptomatic
    dxdt.col(3) = (epsilon * (1.0 - sigma_v0) * x.col(1).array()) +  // E -> Ia
                  (psi * x.col(2).array()) -                         // Is -> Ia
                  (gamma * x.col(3).array()) +                       // Ia -> R
                  reinfections + (aging * x.col(3)).array();
    // change in recovered
    dxdt.col(4) = (gamma * x.col(3).array()) -  // Ia -> R
                  (delta * x.col(4).array()) -  // R -> S
                  reinfections + (aging * x.col(4)).array() +
                  (upsilon_1 * x.col(11).array()) -  // RV1 -> R
                  (phi_1 * x.col(4).array());        // R -> RV1

    // change in vaccinated one dose
    dxdt.col(7) = -new_infections_v1 + (aging * x.col(7)).array() +
                  (phi_1 * x.col(0).array()) +               // S -> SV1
                  (upsilon_2 * x.col(14).array()) +          // SV2 -> SV1
                  (delta * x.col(11).array()) +              // RV1 -> SV1
                  (delta * upsilon_2 * x.col(18).array()) -  // RV2 -> SV1
                  (phi_2 * x.col(7).array()) -               // SV1 -> SV2
                  (upsilon_1 * x.col(7).array());            // SV1 -> S

    // change in exposed
    dxdt.col(8) = new_infections_v1 + (aging * x.col(8)).array() -
                  (epsilon * x.col(8).array());  // EV -> IsV1 + IaV1
    // change in infectious symptomatic
    dxdt.col(9) = (epsilon * sigma_v1 * x.col(8).array()) -  // EV1 -> IsV1
                  (psi * x.col(9).array()) +                 // IsV1 -> IaV1
                  (aging * x.col(9)).array();

    // change in infectious asymptomatic
    dxdt.col(10) =
        (epsilon * (1.0 - sigma_v1) * x.col(8).array()) +  // EV1 -> IaV1
        (psi * x.col(9).array()) -                         // IsV1 -> IaV1
        (gamma * x.col(10).array()) +                      // IaV1 -> RV1
        reinfections_v1 +                                  // RV1 -> IaV1
        (aging * x.col(10)).array();

    // change in recovered
    dxdt.col(11) = (gamma * x.col(10).array()) -  // IaV1 -> RV1
                   (delta * x.col(11).array()) -  // RV1 -> SV1
                   reinfections_v1 + (aging * x.col(11)).array() -
                   (upsilon_1 * x.col(11).array()) -            // RV1 -> R
                   (phi_2 * x.col(11).array()) -                // RV1 -> RV2
                   (x.col(11).array() * (delta * upsilon_1)) +  // RV1 -> S
                   (phi_1 * x.col(4).array()) +                 // R -> RV1
                   (upsilon_2 * x.col(18).array());             // RV2 -> RV1

    // change in vaccinated two doses
    dxdt.col(14) = -new_infections_v2 + (aging * x.col(14)).array() +
                   (phi_2 * x.col(7).array()) -       // SV1 -> SV2
                   (upsilon_2 * x.col(14).array()) +  // SV2 -> SV1
                   (delta * x.col(18).array());       // RV2 -> SV2

    // change in exposed
    dxdt.col(15) = new_infections_v2 + (aging * x.col(15)).array() -
                   (epsilon * x.col(15).array());  // EV2 -> IsV2 + IaV2
    // change in infectious symptomatic
    dxdt.col(16) = (epsilon * sigma_v2 * x.col(15).array()) -  // EV2 -> IsV2
                   (psi * x.col(16).array()) +                 // IsV2 -> IaV2
                   (aging * x.col(16)).array();

    // change in infectious asymptomatic
    dxdt.col(17) =
        (epsilon * (1.0 - sigma_v2) * x.col(15).array()) +  // EV2 -> IaV2
        (psi * x.col(16).array()) -                         // IsV2 -> IaV2
        (gamma * x.col(17).array()) +                       // IaV2 -> RV2
        reinfections_v2 +                                   // RV2 -> IaV2
        (aging * x.col(17)).array();

    // change in recovered
    dxdt.col(18) = (gamma * x.col(17).array()) -  // IaV2 -> RV2
                   (delta * x.col(18).array()) -  // RV2 -> SV2
                   reinfections_v2 + (aging * x.col(18)).array() -
                   (upsilon_2 * x.col(18).array()) -            // RV2 -> Rv1
                   (x.col(18).array() * (delta * upsilon_2)) +  // RV2 -> SV1
                   (phi_2 * x.col(11).array());                 // RV1 -> RV2

    // mortality in all compartments
    dxdt(Eigen::all, epi_indices) =
        dxdt(Eigen::all, epi_indices) * (1.0 - d[0]);
    // change in new infections and re-infections
    dxdt.col(5) = new_infections;
    dxdt.col(6) = reinfections;
    dxdt.col(12) = new_infections_v1;
    dxdt.col(13) = reinfections_v1;
    dxdt.col(19) = new_infections_v2;
    dxdt.col(20) = reinfections_v2;
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
