#' Handle data from the Norovirus Boost model
#'
#' @param output An output from the function [noromod_cpp_boost()], which is a
#' list with the elements `x` for the model state, and `time` for the simulation
#' time.
#'
#' @return A data.frame with columns named `<compartment>_<age_group>`, and time
#' for the state at each time point.
#' @export
output_to_df <- function(output) {
  stopifnot(
    "`output` must be a list with the names 'x' and 'time'" =
      (is.list(output) && setequal(c("x", "time"), names(output)))
  )

  state <- lapply(output[["x"]], as.vector)
  state <- do.call(what = "rbind", state)
  state <- as.data.frame(state)

  names_df <- expand.grid(
    age_group = seq(4),
    compartments = c(
      "susceptible", "exposed", "infectious_symp", "infectious_asymp",
      "recovered", "reinfections", "new_infectious"
    )
  )

  colnames(state) <- sprintf(
    "%s_%i",
    names_df$compartments,
    names_df$age_group
  )

  state$time <- output[["time"]]

  # return state
  state
}
