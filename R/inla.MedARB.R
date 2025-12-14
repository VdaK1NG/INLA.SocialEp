#' Calculate MedARB between fitted values in INLA and simulated values
#'
#' This function estimates Median Absolute Relative Bias between fitted values in INLA and simulated values
#'
#' @param fit_values Mean fitted values from INLA models
#' @param sim_values Vector of simulated values
#' @param n.sim Number of simulated datasets
#' @return Value estimated for MedARB
#' @export
#'

inla.MedARB <- function(fit_values, sim_values, n.sim=1){

  # Make sure fitted values have the same length as simulated values provided
  if(length(fit_values)!=length(sim_values)){stop("Vector of simulated values has different length than the vector of fitted values")}

  # Estimate MedARB
  MedARB <- median(abs(Reduce("+",mapply(function(x,y){(x-y)/y}, x=fit_values, y=sim_values, SIMPLIFY=FALSE)))/n.sim)

  # Return values
  return(MedARB)
}
