
<!-- README.md is generated from README.Rmd. Please edit that file -->

# noromod

<!-- badges: start -->

[![R-CMD-check](https://github.com/pratikunterwegs/noromod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pratikunterwegs/noromod/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of noromod is to provide a three implementations of a model of
age-stratified norovirus transmission.

## Installation

You can install the development version of noromod from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
pak::pak("pratikunterwegs/noromod")
```

## Example

This is a basic example which shows you how to use the three
implementations.

Prepare parameters.

``` r
library(noromod)
library(deSolve)
```

``` r
# define parameters
# initial conditions
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
#> Removing participants that have contacts without age information. To change this behaviour, set the 'missing.contact.age' option

# Symmetrical contact matrix
uk_contact_rate_matrix <- as.matrix(UK_structure$matrix)
demography <- UK_structure$demography$population

uk_contact_rate_matrix <- t(t(uk_contact_rate_matrix) / demography)

# add contact matrix to pop
params <- default_parameters()
params[["contacts"]] <- uk_contact_rate_matrix

# time points
times <- seq(0, 11000)
```

### Using the R-only model with deSolve

``` r
data <- as.data.frame(deSolve::lsoda(init, times, norovirus_model_r, params))
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go

plot(rowSums(data[, seq(6, 9)]), type = "l")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Using Boost solvers

``` r
# Using Boost solvers for increased speed
# run model
data <- noromod_cpp_boost(
  initial_conditions = init_mat,
  params = params, time_end = 11000, increment = 1
)

data <- do.call(rbind, lapply(data[["x"]], as.vector))

plot(rowSums(data[, seq(5, 8)]), type = "l")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Speed comparison

``` r
microbenchmark::microbenchmark(
  "noromod_r" = deSolve::lsoda(init, times, norovirus_model_r, params),
  "noromod_cpp" = deSolve::lsoda(init, times, norovirus_model_cpp, params),
  "noromod_cpp_boost" = noromod_cpp_boost(
    initial_conditions = init_mat,
    params = params, time_end = max(times), increment = 1
  )
)
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> DLSODA-  At T (=R1), too much accuracy requested  
#>       for precision of machine..  See TOLSF (=R2) 
#> In above message, R1 = 10771.6, R2 = nan
#> 
#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Excessive
#> precision requested.  scale up `rtol' and `atol' e.g by the factor 10

#> Warning in deSolve::lsoda(init, times, norovirus_model_r, params): Returning
#> early. Results are accurate, as far as they go
#> Unit: milliseconds
#>               expr        min        lq      mean    median        uq       max
#>          noromod_r 1188.64930 1236.1780 1333.8733 1283.4355 1358.7774 2554.5153
#>        noromod_cpp  356.43054  389.9030  435.4173  416.0641  454.5616  814.8951
#>  noromod_cpp_boost   99.30052  102.4702  115.9454  109.7454  117.4761  281.5254
#>  neval
#>    100
#>    100
#>    100
```
