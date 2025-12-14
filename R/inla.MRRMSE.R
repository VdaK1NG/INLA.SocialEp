#' Calculate MRRMSE between fitted values in INLA and simulated values
#'
#' This function estimates Mean Relative Root Mean Prediction Error between fitted values in INLA and simulated values
#'
#' @param fit_values Mean fitted values from INLA models
#' @param sim_values Vector of simulated values
#' @param n.sim Number of simulated datasets
#' @return Value estimated for MRRMSE
#' @export
#'

inla.MRRMSE <- function(fit_values, sim_values, n.sim=1){

  # Make sure fitted values have the same length as simulated values provided
  if(length(fit_values)!=length(sim_values)){stop("Vector of simulated values has different length than the vector of fitted values")}

  # Estimate MRRMSE
  MRRMSE <- median(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, x=adj_risk, y=real_risk, SIMPLIFY=FALSE))/n.sim))

  # Return values
  return(MRRMSE)
}
