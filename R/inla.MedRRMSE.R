#' Calculate MedRRMSE between fitted values in INLA and simulated values
#'
#' This function estimates Median Relative Root Mean Prediction Error between fitted values in INLA and simulated values
#'
#' @param fit_values Mean fitted values from INLA models
#' @param sim_values Vector of simulated values
#' @param n.sim Number of simulated datasets
#' @return Value estimated for MedRRMSE
#' @export
#'

inla.MedRRMSE <- function(fit_values, sim_values, n.sim=1){

  # Make sure fitted values have the same length as simulated values provided
  if(length(fit_values)!=length(sim_values)){stop("Vector of simulated values has different length than the vector of fitted values")}

  # Estimate MedRRMSE
  MedRRMSE <- median(sqrt(Reduce("+",mapply(function(x,y){((x-y)/y)^2}, x=adj_risk, y=real_risk, SIMPLIFY=FALSE))/n.sim))

  # Return values
  return(MedRRMSE)
}
