# noromod 0.1.0

This is an initial minor version of the _noromod_ package, and contains three implementations of a deterministic compartmental model of norovirus transmission with age-stratification.

- `norovirus_model_r()` provides an R-only implementation of the model ODE system which is suitable to be passed to ODE solvers from the _deSolve_ package.

- `norovirus_model_cpp()` provides a version of `norovirus_model_r()` which passes the model parameters, state, and contact matrix to Rcpp for some computational speed gain.

- `noromod_cpp_boost()` provides a performant version of the norovirus model which uses a Boost _odeint_ solver.

All three versions return equivalent results.

The package also provides a set of useful default parameters for the model via the function `default_parameters()`.

_noromod_ relies on the Epiverse-TRACE package _epidemics_ for the state type used in `noromod_cpp_boost()`.
