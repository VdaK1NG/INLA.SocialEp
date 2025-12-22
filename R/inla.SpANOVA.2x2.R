#' Apply SpANOVA modelization using a wrapper function in INLA
#'
#' This function is a wrapper for multidimensional spatial factor models
#' in INLA, using a sequential shared spatial effects with nested structure
#' as discussed in AÑADIR REFERENCIA. <br/>
#' <br/>
#' <b>DISCLAIMER:</b> Observed and expected values have to be given in an specific order.
#' Consider n the number of areas, the first n values (1:n) of the obs (exp) should be the ones belonging to the
#' FIRST level (the first position of the <b>lev.fac1</b> vector argument) of the
#' FIRST factor (the first position of the <b>fac.names</b> vector argument) and
#' to the FIRST level (the first position of the <b>lev.fac2</b> vector argument)
#' of the SECOND factor (the second position of the <b>fac.names</b> vector argument).
#' The n following values ((n+1):2n) of the obs (exp) should be the ones
#' belonging to the FIRST level (the first position of the <b>lev.fac1</b> vector
#' argument) of the FIRST factor (the first position of the fac.names</b>
#' vector argument) and to the SECOND level (the second position of the
#' <b>lev.fac2</b> vector argument) of the SECOND factor (the second position of
#' the <b>fac.names</b> vector argument).<br/>
#' <br/>
#' The n following values ((2n+1):3n) of the obs (exp) should be the ones
#' belonging to the SECOND level (the second position of the <b>lev.fac1</b> vector
#' argument) of the FIRST factor (the first position of the <b>fac.names</b>
#' vector argument) and to the FIRST level (the first position of the
#' <b>lev.fac2</b> vector argument) of the SECOND factor (the second position of
#' the <b>fac.names</b> vector argument).<br/>
#' <br/>
#' The n following values ((3n+1):4n) of the obs (exp) should be the ones
#' belonging to the SECOND level (the second position of the <b>lev.fac1</b> vector
#' argument) of the FIRST factor (the first position of the <b>fac.names</b>
#' vector argument) and to the SECOND level (the second position of the
#' <b>lev.fac2</b> vector argument) of the SECOND factor (the second position of
#' the <b>fac.names</b> vector argument).<br/>
#' <br/>
#' The first n values are O1, the second n values are O2, the third n values
#' are O3 and the last n values are O4
#'
#' @param obs Vector of observed values
#' @param exp Vector of expected values
#' @param gr Graph for the underlying spatial structure
#' @param fac.names Names of the factors included
#' @param lev.fac1 Levels of the first factor included
#' @param lev.fac2 Levels of the second factor included
#' @param scale.mod Scale copies of random spatial effects or not, default is TRUE
#' @param sp.prior Select prior for the random spatial effect, options are sdunif and pc.prec
#' @param pc.prec.val Define values por the pc prior in case it was chosen. Default values are c(1, 0.01)
#' @param sp.copy.fixed Fix copied values for the random spatial effects, default is TRUE
#' @param save.res Save fitted values or not from the different models, default is TRUE
#' @param save.random Save values adjusted from the random spatial effects or not from the different models, default is TRUE
#' @param save.hyper Save hyperparameter values from each individual model, default is TRUE
#' @param save.mod.data Save modelling data to run the model afterwards
#' @return List with all the models analyzed and a summary table with the most common performance metrics.
#' @export

inla.SpANOVA.2x2 <- function(obs, exp, gr, fac.names = NULL, lev.fac1 = NULL, lev.fac2 = NULL, scale.mod=TRUE, sp.prior="sdunif", pc.prec.val = c(1, 0.01),
                             sp.copy.fixed=TRUE, save.res=TRUE, save.random=TRUE, save.hyper=TRUE, save.mod.data=TRUE) {

  ## Print warnings (warnings are printed as they occur)
  options(warn=1)

  ## Define uniform priors for the standard deviation
  if(sp.prior=="sdunif"){
    prior.prec <- list(prior = "expression: log_dens = 0 - log(2) - theta / 2; return(log_dens);", initial = 0)
  }else if(sp.prior=="pc.prec"){prior.prec <- list(prior = "pc.prec", param = pc.prec.val)}

  ## Define sp model
  sp.model <- "besag"

  ## Define number of levels and number of areas
  n.levels <- c(2, 2)
  n.areas <- gr$n

  ## Get maximum number of groups/diseases
  n.groups <- base::prod(n.levels)

  ## Check if obs and exp are the same
  #if (sum(obs) != as.integer(sum(exp))){
  #  stop("ERROR: Number of total observations and expected values are not the same.")
  #}

  ## Check if data has the proper length
  if (length(obs) != n.areas*n.groups){
    stop("ERROR: The length of the observations does not correspond with the number of groups and areas.")
  }

  if (length(exp) != n.areas*n.groups){
    stop("ERROR: The length of the expected values does not correspond with the number of groups and areas.")
  }

  ## Empty list
  data.INLA <- list(OBS_f1l1_f2l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l1_f2l1 = matrix(NA, nrow = n.areas),
                    OBS_f1l2_f2l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l2_f2l1 = matrix(NA, nrow = n.areas),
                    OBS_f1l1_f2l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l1_f2l2 = matrix(NA, nrow = n.areas),
                    OBS_f1l2_f2l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f1l2_f2l2 = matrix(NA, nrow = n.areas),
                    OBS_f2l1_f1l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l1_f1l1 = matrix(NA, nrow = n.areas),
                    OBS_f2l2_f1l1 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l2_f1l1 = matrix(NA, nrow = n.areas),
                    OBS_f2l1_f1l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l1_f1l2 = matrix(NA, nrow = n.areas),
                    OBS_f2l2_f1l2 = matrix(NA, nrow = n.areas*n.groups, ncol = n.groups), EXP_f2l2_f1l2 = matrix(NA, nrow = n.areas))

  ## Create empty lists where info about the models will be stored
  data.models <- list()

  ## Here we order the Observed and Expected values following the default order: ------------------
  # O1, O2, O3 and O4
  data.INLA$OBS_f1l1_f2l1[1:n.areas, 1] <- obs[1:n.areas] # O1
  data.INLA$OBS_f1l1_f2l1[n.areas + 1:n.areas, 2] <- obs[n.areas + 1:n.areas] # O2
  data.INLA$OBS_f1l1_f2l1[2*n.areas + 1:n.areas, 3] <- obs[2*n.areas + 1:n.areas] # O3
  data.INLA$OBS_f1l1_f2l1[3*n.areas + 1:n.areas, 4] <- obs[3*n.areas + 1:n.areas] # O4
  # E1, E2, E3 and E4
  data.INLA$EXP_f1l1_f2l1[1:n.areas] <- exp[1:n.areas] # E1
  data.INLA$EXP_f1l1_f2l1[n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  data.INLA$EXP_f1l1_f2l1[2*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  data.INLA$EXP_f1l1_f2l1[3*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4

  ## Here we order the Observed and Expected values following this order: ------------------
  # O3, O4, O1 and O2
  data.INLA$OBS_f1l2_f2l1[1:n.areas, 1] <- obs[2*n.areas + 1:n.areas] # O3
  data.INLA$OBS_f1l2_f2l1[n.areas + 1:n.areas, 2] <- obs[3*n.areas + 1:n.areas] # O4
  data.INLA$OBS_f1l2_f2l1[2*n.areas + 1:n.areas, 3] <- obs[1:n.areas] # O1
  data.INLA$OBS_f1l2_f2l1[3*n.areas + 1:n.areas, 4] <- obs[n.areas + 1:n.areas] # O2
  # E3, E4, E1 and E2
  data.INLA$EXP_f1l2_f2l1[1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  data.INLA$EXP_f1l2_f2l1[n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f1l2_f2l1[2*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  data.INLA$EXP_f1l2_f2l1[3*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2

  ## Here we order the Observed and Expected values following this order: ------------------
  # O2, O1, O4 and O3
  data.INLA$OBS_f1l1_f2l2[1:n.areas, 1] <- obs[n.areas + 1:n.areas]  # O2
  data.INLA$OBS_f1l1_f2l2[n.areas + 1:n.areas, 2] <- obs[1:n.areas] # O1
  data.INLA$OBS_f1l1_f2l2[2*n.areas + 1:n.areas, 3] <- obs[3*n.areas + 1:n.areas] # O4
  data.INLA$OBS_f1l1_f2l2[3*n.areas + 1:n.areas, 4] <- obs[2*n.areas + 1:n.areas] # O3
  # E2, E1, E4 and E3
  data.INLA$EXP_f1l1_f2l2[1:n.areas, 1] <- exp[n.areas + 1:n.areas]  # E2
  data.INLA$EXP_f1l1_f2l2[n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  data.INLA$EXP_f1l1_f2l2[2*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f1l1_f2l2[3*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3

  ## Here we order the Observed and Expected values following this order: ------------------
  # O4, O3, O2 and O1
  data.INLA$OBS_f1l2_f2l2[1:n.areas, 1] <- obs[3*n.areas + 1:n.areas] # E4
  data.INLA$OBS_f1l2_f2l2[n.areas + 1:n.areas, 2] <- obs[2*n.areas + 1:n.areas] # E3
  data.INLA$OBS_f1l2_f2l2[2*n.areas + 1:n.areas, 3] <- obs[n.areas + 1:n.areas] # E2
  data.INLA$OBS_f1l2_f2l2[3*n.areas + 1:n.areas, 4] <- obs[1:n.areas] # E1
  # E4, E3, E2 and E1
  data.INLA$EXP_f1l2_f2l2[1:n.areas, 1] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f1l2_f2l2[n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  data.INLA$EXP_f1l2_f2l2[2*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  data.INLA$EXP_f1l2_f2l2[3*n.areas + 1:n.areas] <- exp[1:n.areas] # E1

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

  ## Here the observed are ordered following this order: ------------------
  # O1, O3, O2 and O4
  data.INLA$OBS_f2l1_f1l1[1:n.areas, 1] <- obs[1:n.areas] # O1
  data.INLA$OBS_f2l1_f1l1[n.areas + 1:n.areas, 2] <- obs[2*n.areas + 1:n.areas] # O3
  data.INLA$OBS_f2l1_f1l1[2*n.areas + 1:n.areas, 3] <- obs[n.areas + 1:n.areas] # O2
  data.INLA$OBS_f2l1_f1l1[3*n.areas + 1:n.areas, 4] <- obs[3*n.areas + 1:n.areas]  # O4
  # E1, E3, E2 and E4
  data.INLA$EXP_f2l1_f1l1[1:n.areas] <- exp[1:n.areas] # O1
  data.INLA$EXP_f2l1_f1l1[n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # O3
  data.INLA$EXP_f2l1_f1l1[2*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # O2
  data.INLA$EXP_f2l1_f1l1[3*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas]  # O4

  ## Here the observed are ordered following this order: ------------------
  # O2, O4, O1 and O3
  data.INLA$OBS_f2l2_f1l1[1:n.areas, 1] <- obs[n.areas + 1:n.areas] # O2
  data.INLA$OBS_f2l2_f1l1[n.areas + 1:n.areas, 2] <- obs[3*n.areas + 1:n.areas] # O4
  data.INLA$OBS_f2l2_f1l1[2*n.areas + 1:n.areas, 3] <- obs[1:n.areas] # O1
  data.INLA$OBS_f2l2_f1l1[3*n.areas + 1:n.areas, 4] <- obs[2*n.areas + 1:n.areas] # O3
  # E2, E4, E1 and E3
  data.INLA$EXP_f2l2_f1l1[1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  data.INLA$EXP_f2l2_f1l1[n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f2l2_f1l1[2*n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  data.INLA$EXP_f2l2_f1l1[3*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3

  ## Here the observed are ordered following this order: ------------------
  # O3, O1, O4 and O2
  data.INLA$OBS_f2l1_f1l2[1:n.areas, 1] <- obs[2*n.areas + 1:n.areas] # O3
  data.INLA$OBS_f2l1_f1l2[n.areas + 1:n.areas, 2] <- obs[1:n.areas] # O1
  data.INLA$OBS_f2l1_f1l2[2*n.areas + 1:n.areas, 3] <- obs[3*n.areas + 1:n.areas] # O4
  data.INLA$OBS_f2l1_f1l2[3*n.areas + 1:n.areas, 4] <- obs[n.areas + 1:n.areas] # O2
  # E3, E1, E4 and E2
  data.INLA$EXP_f2l1_f1l2[1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  data.INLA$EXP_f2l1_f1l2[n.areas + 1:n.areas] <- exp[1:n.areas] # E1
  data.INLA$EXP_f2l1_f1l2[2*n.areas + 1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f2l1_f1l2[3*n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2

  ## Here the observed are ordered following this order: ------------------
  # O4, O2, O3 and O1
  data.INLA$OBS_f2l2_f1l2[1:n.areas, 1] <- obs[3*n.areas + 1:n.areas] # O4
  data.INLA$OBS_f2l2_f1l2[n.areas + 1:n.areas, 2] <- obs[n.areas + 1:n.areas] # O2
  data.INLA$OBS_f2l2_f1l2[2*n.areas + 1:n.areas, 3] <- obs[2*n.areas + 1:n.areas] # O3
  data.INLA$OBS_f2l2_f1l2[3*n.areas + 1:n.areas, 4] <- obs[1:n.areas] # O1
  # E4, E2, E3 and E1
  data.INLA$EXP_f2l2_f1l2[1:n.areas] <- exp[3*n.areas + 1:n.areas] # E4
  data.INLA$EXP_f2l2_f1l2[n.areas + 1:n.areas] <- exp[n.areas + 1:n.areas] # E2
  data.INLA$EXP_f2l2_f1l2[2*n.areas + 1:n.areas] <- exp[2*n.areas + 1:n.areas] # E3
  data.INLA$EXP_f2l2_f1l2[3*n.areas + 1:n.areas] <- exp[1:n.areas] # E1

  ## Intercepts for each group ------------------
  # Here the sequence is: "f1l1-f2l1" "f1l1-f2l2" "f1l2-f2l1" "f1l2-f2l2"
  data.INLA$alpha_f1l1_f2l1 <- rep(c(paste0(lev.fac1[1], "-", lev.fac2[1]), paste0(lev.fac1[1], "-", lev.fac2[2]),
                                     paste0(lev.fac1[2], "-", lev.fac2[1]), paste0(lev.fac1[2], "-", lev.fac2[2])),
                                   each = n.areas)
  data.INLA$alpha_f1l1_f2l1 <- as.factor(data.INLA$alpha_f1l1_f2l1)

  # Here the sequence is: "f1l2-f2l1" "f1l2-f2l2" "f1l1-f2l1" "f1l1-f2l2"
  data.INLA$alpha_f1l2_f2l1 <- rep(c(paste0(lev.fac1[2], "-", lev.fac2[1]), paste0(lev.fac1[2], "-", lev.fac2[2]),
                                     paste0(lev.fac1[1], "-", lev.fac2[1]), paste0(lev.fac1[1], "-", lev.fac2[2])),
                                   each = n.areas)
  data.INLA$alpha_f1l2_f2l1 <- as.factor(data.INLA$alpha_f1l2_f2l1)

  # Here the sequence is: "f1l1-f2l2" "f1l1-f2l1" "f1l2-f2l2" "f1l2-f2l1"
  data.INLA$alpha_f1l1_f2l2 <- rep(c(paste0(lev.fac1[1], "-", lev.fac2[2]), paste0(lev.fac1[1], "-", lev.fac2[1]),
                                     paste0(lev.fac1[2], "-", lev.fac2[2]), paste0(lev.fac1[2], "-", lev.fac2[1])),
                                   each = n.areas)
  data.INLA$alpha_f1l1_f2l2 <- as.factor(data.INLA$alpha_f1l1_f2l2)

  # Here the sequence is: "f1l2-f2l2" "f1l2-f2l1" "f1l1-f2l2" "f1l1-f2l1"
  data.INLA$alpha_f1l2_f2l2 <- rep(c(paste0(lev.fac1[2], "-", lev.fac2[2]), paste0(lev.fac1[2], "-", lev.fac2[1]),
                                     paste0(lev.fac1[1], "-", lev.fac2[2]), paste0(lev.fac1[1], "-", lev.fac2[1])),
                                   each = n.areas)
  data.INLA$alpha_f1l2_f2l2 <- as.factor(data.INLA$alpha_f1l2_f2l2)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

  # Here the sequence is: "f2l1-f1l1" "f2l1-f1l2" "f2l2-f1l1" "f2l2-f1l2"
  data.INLA$alpha_f2l1_f1l1 <- rep(c(paste0(lev.fac2[1], "-", lev.fac1[1]), paste0(lev.fac2[1], "-", lev.fac1[2]),
                                     paste0(lev.fac2[2], "-", lev.fac1[1]), paste0(lev.fac2[2], "-", lev.fac1[2])),
                                   each = n.areas)
  data.INLA$alpha_f2l1_f1l1 <- as.factor(data.INLA$alpha_f2l1_f1l1)

  # Here the sequence is: "f2l2-f1l1" "f2l2-f1l2" "f2l1-f1l1" "f2l1-f1l2"
  data.INLA$alpha_f2l2_f1l1 <- rep(c(paste0(lev.fac2[2], "-", lev.fac1[1]), paste0(lev.fac2[2], "-", lev.fac1[2]),
                                     paste0(lev.fac2[1], "-", lev.fac1[1]), paste0(lev.fac2[1], "-", lev.fac1[2])),
                                   each = n.areas)
  data.INLA$alpha_f2l2_f1l1 <- as.factor(data.INLA$alpha_f2l2_f1l1)

  # Here the sequence is: "f2l1-f1l2" "f2l1-f1l1" "f2l2-f1l2" "f2l2-f1l1"
  data.INLA$alpha_f2l1_f1l2 <- rep(c(paste0(lev.fac2[1], "-", lev.fac1[2]), paste0(lev.fac2[1], "-", lev.fac1[1]),
                                     paste0(lev.fac2[2], "-", lev.fac1[2]), paste0(lev.fac2[2], "-", lev.fac1[1])),
                                   each = n.areas)
  data.INLA$alpha_f2l1_f1l2 <- as.factor(data.INLA$alpha_f2l1_f1l2)

  # Here the sequence is: "f2l2-f1l2" "f2l2-f1l1" "f2l1-f1l2" "f2l1-f1l1"
  data.INLA$alpha_f2l2_f1l2 <- rep(c(paste0(lev.fac2[2], "-", lev.fac1[2]), paste0(lev.fac2[2], "-", lev.fac1[1]),
                                     paste0(lev.fac2[1], "-", lev.fac1[2]), paste0(lev.fac2[1], "-", lev.fac1[1])),
                                   each = n.areas)
  data.INLA$alpha_f2l2_f1l2 <- as.factor(data.INLA$alpha_f2l2_f1l2)


  ## Create IDs for each area and group ------------------
  data.INLA$AREAID <- rep(1:n.areas, n.groups)

  ## Create IDs for shared spatial effects - phi ------------------
  data.INLA$ID_g1 <- data.INLA$AREAID
  data.INLA$ID_g1[-c(1:n.areas)] <- NA
  data.INLA$ID_g2 <- data.INLA$AREAID
  data.INLA$ID_g2[-(n.areas + 1:n.areas)] <- NA
  data.INLA$ID_g3 <- data.INLA$AREAID
  data.INLA$ID_g3[-(2 * n.areas + 1:n.areas)] <- NA
  data.INLA$ID_g4 <- data.INLA$AREAID
  data.INLA$ID_g4[-(3 * n.areas + 1:n.areas)] <- NA

  ## Create IDs for IID Effect - omega_j ------------------
  data.INLA$omega_j <- 1:(n.areas*4)

  ####################################################### MODEL 0 #######################################################
  # This scenario considers a different intercept in the linear predictor of each group
  # and one individual IID Effect. (ResMod in the article)

  # Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +

    # Intercept
    alpha_f1l1_f2l1 +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  # Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-8)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$RME <- ResMod$summary.fitted.values$mean
    rm(ResMod)

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M0"
  rm(list_temp)

  cat("-/-/-/-/-/-/-/-/-/- M0 finished -/-/-/-/-/-/-/-/-/- \n")

  ####################################################### MODEL 1 #######################################################
  # This scenario considers only a different intercept and a different spatial
  # effect in the linear predictor of each group. (ResMod in the article)

  # Individual spatial effects - phi_g
  data.INLA$phi_1 <- data.INLA$ID_g1
  data.INLA$phi_2 <- data.INLA$ID_g2
  data.INLA$phi_3 <- data.INLA$ID_g3
  data.INLA$phi_4 <- data.INLA$ID_g4

  # Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +

    # Intercept
    alpha_f1l1_f2l1 +

    # Individual spatial effects - phi_g
    f(phi_1, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_2, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_3, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_4, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  # Run Model
  try(ResMod <- INLA::inla(formula = formula,data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-8)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M1"
  rm(list_temp)

  cat("-/-/-/-/-/-/-/-/-/- M1 finished -/-/-/-/-/-/-/-/-/- \n")

  ####################################################### MODEL 2 #######################################################
  # This scenario considers a different intercept and the same shared spatial
  # effect in the linear predictor of each group. (M2 in the article)
  # Four different combinations have to be fit where each group is considered as
  # a reference group in each of the four models.

  # Define IDs for overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4

  # IID Effect - omega_j
  data.INLA$omega_j <- 1:(n.areas*4)

  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +

    # Intercept
    alpha_f1l1_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M2"
  rm(list_temp)

  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l1 ~ 0 +

    # Intercept
    alpha_f1l2_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l1
    list_temp$obs <- OBS_f1l2_f2l1
    list_temp$exp <- EXP_f1l2_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M2"
  rm(list_temp)

  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l2 ~ 0 +

    # Intercept
    alpha_f1l1_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l2
    list_temp$obs <- OBS_f1l1_f2l2
    list_temp$exp <- EXP_f1l1_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M2"
  rm(list_temp)

  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l2 ~ 0 +

    # Intercept
    alpha_f1l2_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l2
    list_temp$obs <- OBS_f1l2_f2l2
    list_temp$exp <- EXP_f1l2_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    ## Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M2"
  rm(list_temp)

  cat("-/-/-/-/-/-/-/-/-/- M2 finished -/-/-/-/-/-/-/-/-/- \n")

  ####################################################### MODEL 3 #######################################################
  # This scenario considers a different intercept, the same shared spatial
  # effect in the linear predictor of each group and a specific efect for the
  # second category of the first factor (M3 in the article). Two different
  # combinations have to be fit where each of the two levels is considered
  # as a reference category of the first factor.

  # Overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4

  # Define IDs for 2-group shared effect - phi_21
  data.INLA$phi_21  <- data.INLA$ID_g3
  data.INLA$phi_21_g4 <- data.INLA$ID_g4

  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +

    # Intercept
    alpha_f1l1_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M3"
  rm(list_temp)

  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l1 ~ 0 +

    # Intercept
    alpha_f1l2_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l1
    list_temp$obs <- OBS_f1l2_f2l1
    list_temp$exp <- EXP_f1l2_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M3"
  rm(list_temp)

  cat("-/-/-/-/-/-/-/-/-/- M3 finished -/-/-/-/-/-/-/-/-/- \n")

  ####################################################### MODEL 4 #######################################################
  # This scenario considers a different intercept, the same shared spatial
  # effect in the linear predictor of each group and a specific efect for the
  # second category of the second factor (M4 in the article). Two different
  # combinations have to be fit where each of the two levels is considered
  # as a reference category of the second factor.

  # Define IDs for 2-group shared effect - phi_21
  data.INLA$phi_12  <- data.INLA$ID_g2
  data.INLA$phi_12_g4 <- data.INLA$ID_g4

  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +

    # Intercept
    alpha_f1l1_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M4"
  rm(list_temp)

  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l2 ~ 0 +

    # Intercept
    alpha_f1l1_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l2
    list_temp$obs <- OBS_f1l1_f2l2
    list_temp$exp <- EXP_f1l1_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M4"
  rm(list_temp)

  cat("-/-/-/-/-/-/-/-/-/- M4 finished -/-/-/-/-/-/-/-/-/- \n")

  ####################################################### MODEL 5 #######################################################
  # This scenario considers a different intercept, the same shared spatial
  # effect in the linear predictor of each group, a specific spatial effect for
  # the second category of the first factor and a diferent specific spatial
  # effect the second category of the second factor. (M5 in the article)
  # Four different combinations have to be fit.

  pb <- txtProgressBar(min = 0, max = 4, style = 3,  width = 50, char = "=")

  # Define IDs for overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4

  # Define IDs for 2 groups shared effect - phi_12
  data.INLA$phi_21  <- data.INLA$ID_g3
  data.INLA$phi_21_g4 <- data.INLA$ID_g4

  # Define IDs for 2 groups shared effect - phi_21
  data.INLA$phi_12 <- data.INLA$ID_g2
  data.INLA$phi_12_g4 <- data.INLA$ID_g4

  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +
    # Intercept
    alpha_f1l1_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M5"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 1)

  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l1 ~ 0 +
    # Intercept
    alpha_f1l2_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l1
    list_temp$obs <- OBS_f1l2_f2l1
    list_temp$exp <- EXP_f1l2_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M5"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 2)

  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l2 ~ 0 +
    # Intercept
    alpha_f1l1_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l2
    list_temp$obs <- OBS_f1l1_f2l2
    list_temp$exp <- EXP_f1l1_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }


  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M5"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 3)

  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l2 ~ 0 +
    # Intercept
    alpha_f1l2_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_12_g4, copy = "phi_12", fixed = sp.copy.fixed) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l2
    list_temp$obs <- OBS_f1l2_f2l2
    list_temp$exp <- EXP_f1l2_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M5"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 4)

  cat("\n-/-/-/-/-/-/-/-/-/- M5 finished -/-/-/-/-/-/-/-/-/-\n")

  ####################################################### MODEL 6 #######################################################
  # This scenario considers a different intercept, the same shared spatial
  # effect in the linear predictor of each group, a specific spatial effect for
  # the second category of the first factor,  a diferent specific spatial
  # effect the second category of the second factor (just in one group), and
  # another specific spatial effect for the interaction between two factors.
  # (M6 in the article)
  # Eight different combinations have to be fit.

  pb <- txtProgressBar(min = 0, max = 8, style = 3,  width = 50, char = "=")

  # Overall shared spatial effect - phi_11
  data.INLA$phi_11 <- data.INLA$ID_g1
  data.INLA$phi_11_g2 <- data.INLA$ID_g2
  data.INLA$phi_11_g3 <- data.INLA$ID_g3
  data.INLA$phi_11_g4 <- data.INLA$ID_g4

  # # Define IDs for 2 groups shared effect - phi_12
  data.INLA$phi_21  <- data.INLA$ID_g3
  data.INLA$phi_21_g4 <- data.INLA$ID_g4

  # # Define IDs for individual effect 1 - phi_21
  data.INLA$phi_12 <- data.INLA$ID_g2

  # # Define IDs for individual effect 2 - phi_22
  data.INLA$phi_22 <- data.INLA$ID_g4

  ### F1L1-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l1 ~ 0 +
    # Intercept
    alpha_f1l1_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l1
    list_temp$obs <- OBS_f1l1_f2l1
    list_temp$exp <- EXP_f1l1_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 1)

  ### F1L2-F2L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l1 ~ 0 +
    # Intercept
    alpha_f1l2_f2l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l1
    list_temp$obs <- OBS_f1l2_f2l1
    list_temp$exp <- EXP_f1l2_f2l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 2)

  ### F1L1-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l1_f2l2 ~ 0 +
    # Intercept
    alpha_f1l1_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l1_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l1_f2l2
    list_temp$obs <- OBS_f1l1_f2l2
    list_temp$exp <- EXP_f1l1_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 3)

  ### F1L2-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f1l2_f2l2 ~ 0 +
    # Intercept
    alpha_f1l2_f2l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f1l2_f2l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f1l2_f2l2
    list_temp$obs <- OBS_f1l2_f2l2
    list_temp$exp <- EXP_f1l2_f2l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 4)

  ### F2L1-F1L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f2l1_f1l1 ~ 0 +
    # Intercept
    alpha_f2l1_f1l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l1_f1l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
      list_temp$formula <- formula
      list_temp$alpha <- alpha_f2l1_f1l1
      list_temp$obs <- OBS_f2l1_f1l1
      list_temp$exp <- EXP_f2l1_f1l1
    }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 5)

  ### F2L2-F1L1 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f2l2_f1l1 ~ 0 +
    # Intercept
    alpha_f2l2_f1l1 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l2_f1l1,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f2l2_f1l1
    list_temp$obs <- OBS_f2l2_f1l1
    list_temp$exp <- EXP_f2l2_f1l1
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 6)

  ### F2L1-F1L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f2l1_f1l2 ~ 0 +
    # Intercept
    alpha_f2l1_f1l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l1_f1l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f2l1_f1l2
    list_temp$obs <- OBS_f2l1_f1l2
    list_temp$exp <- EXP_f2l1_f1l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 7)

  ### F2L2-F2L2 ---------------------------------------------------------------------------------------------------------------------

  ## Formulas for the model
  formula <- OBS_f2l2_f1l2 ~ 0 +
    # Intercept
    alpha_f2l2_f1l2 +

    # Overall shared spatial effect - phi_11
    f(phi_11, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_11_g2, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g3, copy = "phi_11", fixed = sp.copy.fixed) +
    f(phi_11_g4, copy = "phi_11", fixed = sp.copy.fixed) +

    # 2 groups shared effect - phi_21
    f(phi_21, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +
    f(phi_21_g4, copy = "phi_21", fixed = sp.copy.fixed) +

    # Individual effect 1 - phi_12
    f(phi_12, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # Individual effect 1 - phi_22
    f(phi_22, model = sp.model, graph = gr,
      constr=TRUE, rankdef=1,  # Sum-To-Zero Constrain
      scale.model=scale.mod, # Scale spatial effects
      hyper = list(prec = prior.prec) # Prior Distributions
    ) +

    # IID Effect - omega_j
    f(omega_j, model = "iid")

  ## Run model in INLA
  try(ResMod <- INLA::inla(formula = formula, data = data.INLA, family = rep("poisson", 4), E = data.INLA$EXP_f2l2_f1l2,
                     control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                     control.inla = list(tolerance.step = 1e-5)))

  # Save Data of the Model into list of models
  list_temp <- list()

  if(save.mod.data==TRUE){
    list_temp$formula <- formula
    list_temp$alpha <- alpha_f2l2_f1l2
    list_temp$obs <- OBS_f2l2_f1l2
    list_temp$exp <- EXP_f2l2_f1l2
  }

  if(exists("ResMod")){
    list_temp$DIC <- ResMod$dic
    list_temp$WAIC <- ResMod$waic
    list_temp$CPU <- ResMod$cpu.used

    # LOOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=-1)
    list_temp$loocv <- round(-mean(log(data_temp$cv)), 2)

    # LGOCV
    data_temp <- inla.group.cv(ResMod, num.level.sets=3)
    list_temp$lgocv.m3 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=5)
    list_temp$lgocv.m5 <- round(-mean(log(data_temp$cv)), 2)

    data_temp <- inla.group.cv(ResMod, num.level.sets=10)
    list_temp$lgocv.m10 <- round(-mean(log(data_temp$cv)), 2)

    # Extract Residuals
    if(save.res==TRUE){list_temp$residuals <- ResMod$residuals$deviance.residuals}

    # Extract Random Effects
    if(save.random==TRUE){list_temp$summary.random <- ResMod$summary.random}

    # Extract Hyperparameters
    if(save.hyper==TRUE){list_temp$summary.hyperpar <- ResMod$summary.hyperpar}

    # Extract RME for each group
    list_temp$summary.fitted.values <- ResMod$summary.fitted.values

  } else {
    list_temp$DIC$dic <- NA
    list_temp$WAIC$waic <- NA
    list_temp$CPU <- c(NA, NA, NA, NA)
    list_temp$loocv <- NA
    list_temp$lgocv.m3 <- NA
    list_temp$lgocv.m5 <- NA
    list_temp$lgocv.m10 <- NA
    list_temp$residuals <- NA
    list_temp$RME <- NA
    list_temp$summary.random <- NA
    list_temp$summary.hyperpar <- NA
    list_temp$summary.fitted.values <- NA
  }

  # Clean RAM
  data.models[[length(data.models)+1]] <- list_temp
  names(data.models)[length(data.models)] <- "M6"
  rm(list_temp)

  utils::setTxtProgressBar(pb, 8)

  cat("\n-/-/-/-/-/-/-/-/-/- M6 finished -/-/-/-/-/-/-/-/-/-\n")

  ####################################################### RESULTS #######################################################

  ## Add the combination of the groups in the final order to each model for easier interpretation -------------------------------------
  # ResMod
  data.models[[1]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
  # ResMod
  data.models[[2]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))
  # M2
  data.models[[3]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))

  data.models[[4]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))

  data.models[[5]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))

  data.models[[6]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"),
                               paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))
  # M3
  data.models[[7]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))

  data.models[[8]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))
  # M4
  data.models[[9]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                               paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))

  data.models[[10]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"),
                                paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))
  # M5
  data.models[[11]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                                paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))

  data.models[[12]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"),
                                paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))

  data.models[[13]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"),
                                paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))

  data.models[[14]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"),
                                paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))
  # M6
  data.models[[15]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"),
                                paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"))

  data.models[[16]]$groups <- c(paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"),
                                paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"))

  data.models[[17]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[1], " group"), paste0(lev.fac1[2], "-", lev.fac2[2], " group"),
                                paste0(lev.fac1[1], "-", lev.fac2[1], " group"), paste0(lev.fac1[1], "-", lev.fac2[2], " group"))

  data.models[[18]]$groups <- c(paste0(lev.fac1[2], "-", lev.fac2[2], " group"), paste0(lev.fac1[2], "-", lev.fac2[1], " group"),
                                paste0(lev.fac1[1], "-", lev.fac2[2], " group"), paste0(lev.fac1[1], "-", lev.fac2[1], " group"))

  data.models[[19]]$groups <- c(paste0(lev.fac2[1], "-", lev.fac1[1], " group"), paste0(lev.fac2[1], "-", lev.fac1[2], " group"),
                                paste0(lev.fac2[2], "-", lev.fac1[1], " group"), paste0(lev.fac2[2], "-", lev.fac1[2], " group"))

  data.models[[20]]$groups <- c(paste0(lev.fac2[1], "-", lev.fac1[2], " group"), paste0(lev.fac2[1], "-", lev.fac1[1], " group"),
                                paste0(lev.fac2[2], "-", lev.fac1[2], " group"), paste0(lev.fac2[2], "-", lev.fac1[1], " group"))

  data.models[[21]]$groups <- c(paste0(lev.fac2[2], "-", lev.fac1[1], " group"), paste0(lev.fac2[2], "-", lev.fac1[2], " group"),
                                paste0(lev.fac2[1], "-", lev.fac1[1], " group"), paste0(lev.fac2[1], "-", lev.fac1[2], " group"))

  data.models[[22]]$groups <- c(paste0(lev.fac2[2], "-", lev.fac1[2], " group"), paste0(lev.fac2[2], "-", lev.fac1[1], " group"),
                                paste0(lev.fac2[1], "-", lev.fac1[2], " group"), paste0(lev.fac2[1], "-", lev.fac1[1], " group"))

  ## Add name of the level being adjusted in each shared spatial effect ---------------------------------------------------------------

  # ResMod
  data.models[[2]]$sp_effects <- c(paste0("phi_1(", lev.fac1[1], "-", lev.fac2[1], ")"), paste0("phi_2(", lev.fac1[1], "-", lev.fac2[2], ")"),
                                   paste0("phi_3(", lev.fac1[2], "-", lev.fac2[1], ")"), paste0("phi_4(", lev.fac1[2], "-", lev.fac2[2], ")"))

  # M2
  data.models[[3]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "]"))
  data.models[[4]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "]"))
  data.models[[5]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "]"))
  data.models[[6]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "]"))

  # M3
  data.models[[7]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[2], " effect", ")"))

  data.models[[8]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[1], " effect", ")"))

  # M4
  data.models[[9]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_12(", lev.fac2[2], " effect", ")"))

  data.models[[10]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_12(", lev.fac2[1], " effect", ")"))

  # M5
  data.models[[11]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[2], " effect", ")"))

  data.models[[12]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_21(", lev.fac1[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[1], " effect", ")"))

  data.models[[13]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[2], " effect", ")"))

  data.models[[14]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "])"), paste0("phi_21(", lev.fac1[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[1], " effect", ")"))

  # M6
  data.models[[15]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac1[1], "*", lev.fac2[2], " effect", ")"), paste0("phi_22(", lev.fac1[2], "*", lev.fac2[2], " effect", ")"))

  data.models[[16]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[1], "-", lev.fac2[2], "])"), paste0("phi_21(", lev.fac1[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac1[1], "*", lev.fac2[1], " effect", ")"), paste0("phi_22(", lev.fac1[2], "*", lev.fac2[1], " effect", ")"))

  data.models[[17]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[1], "])"), paste0("phi_21(", lev.fac1[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac1[2], "*", lev.fac2[2], " effect", ")"), paste0("phi_22(", lev.fac1[1], "*", lev.fac2[2], " effect", ")"))

  data.models[[18]]$sp_effects <- c(paste0("phi_11(General[", lev.fac1[2], "-", lev.fac2[2], "])"), paste0("phi_21(", lev.fac1[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac1[2], "*", lev.fac2[1], " effect", ")"), paste0("phi_22(", lev.fac1[1], "*", lev.fac2[1], " effect", ")"))


  data.models[[19]]$sp_effects <- c(paste0("phi_11(General[", lev.fac2[1], "-", lev.fac1[1], "])"), paste0("phi_21(", lev.fac2[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[1], "*", lev.fac1[2], " effect", ")"), paste0("phi_22(", lev.fac2[2], "*", lev.fac1[2], " effect", ")"))

  data.models[[20]]$sp_effects <- c(paste0("phi_11(General[", lev.fac2[1], "-", lev.fac1[2], "])"), paste0("phi_21(", lev.fac2[2], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[1], "*", lev.fac1[1], " effect", ")"), paste0("phi_22(", lev.fac2[2], "*", lev.fac1[1], " effect", ")"))

  data.models[[21]]$sp_effects <- c(paste0("phi_11(General[", lev.fac2[2], "-", lev.fac1[1], "])"), paste0("phi_21(", lev.fac2[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[2], "*", lev.fac1[2], " effect", ")"), paste0("phi_22(", lev.fac2[1], "*", lev.fac1[2], " effect", ")"))

  data.models[[22]]$sp_effects <- c(paste0("phi_11(General[", lev.fac2[2], "-", lev.fac1[2], "])"), paste0("phi_21(", lev.fac2[1], " effect", ")"),
                                    paste0("phi_12(", lev.fac2[2], "*", lev.fac1[1], " effect", ")"), paste0("phi_22(", lev.fac2[1], "*", lev.fac1[1], " effect", ")"))

  ## Add name to each model
  names(data.models)[1] <- "M0"
  names(data.models)[2] <- "M1"

  names(data.models)[3] <- paste0("M2.(", lev.fac1[1],")")
  names(data.models)[4] <- paste0("M2.(", lev.fac1[2],")")
  names(data.models)[5] <- paste0("M2.(", lev.fac2[1],")")
  names(data.models)[6] <- paste0("M2.(", lev.fac2[2],")")

  names(data.models)[7] <- paste0("M3.", fac.names[1], "(", lev.fac1[1],")")
  names(data.models)[8] <- paste0("M3.", fac.names[1], "(", lev.fac1[2],")")

  names(data.models)[9] <- paste0("M4.", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[10] <- paste0("M4.", fac.names[2], "(", lev.fac2[2],")")

  names(data.models)[11] <- paste0("M5.", fac.names[1], "(", lev.fac1[1],")", "+", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[12] <- paste0("M5.", fac.names[1], "(", lev.fac1[2],")", "+", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[13] <- paste0("M5.", fac.names[1], "(", lev.fac1[1],")", "+", fac.names[2], "(", lev.fac2[2],")")
  names(data.models)[14] <- paste0("M5.", fac.names[1], "(", lev.fac1[2],")", "+", fac.names[2], "(", lev.fac2[2],")")

  names(data.models)[15] <- paste0("M6.", fac.names[1], "(", lev.fac1[1],")", "*", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[16] <- paste0("M6.", fac.names[1], "(", lev.fac1[2],")", "*", fac.names[2], "(", lev.fac2[1],")")
  names(data.models)[17] <- paste0("M6.", fac.names[1], "(", lev.fac1[1],")", "*", fac.names[2], "(", lev.fac2[2],")")
  names(data.models)[18] <- paste0("M6.", fac.names[1], "(", lev.fac1[2],")", "*", fac.names[2], "(", lev.fac2[2],")")

  names(data.models)[19] <- paste0("M6.", fac.names[2], "(", lev.fac2[1],")", "*", fac.names[1], "(", lev.fac1[1],")")
  names(data.models)[20] <- paste0("M6.", fac.names[2], "(", lev.fac2[2],")", "*", fac.names[1], "(", lev.fac1[1],")")
  names(data.models)[21] <- paste0("M6.", fac.names[2], "(", lev.fac2[1],")", "*", fac.names[1], "(", lev.fac1[2],")")
  names(data.models)[22] <- paste0("M6.", fac.names[2], "(", lev.fac2[2],")", "*", fac.names[1], "(", lev.fac1[2],")")

  data_summary <- data.frame(
    "NUMBER"=1:22,

    "MODEL"=c("M0", "M1",
              paste0("M2-ind(", lev.fac1[1], "-", lev.fac2[1], ")"), paste0("M2-ind(", lev.fac1[2], "-", lev.fac2[1], ")"),
              paste0("M2-ind(", lev.fac2[1], "-", lev.fac1[2], ")"), paste0("M2-ind(", lev.fac2[2], "-", lev.fac1[2], ")"),

              paste0("M3-", fac.names[1], ".(", lev.fac1[1],")"), paste0("M3-", fac.names[1], ".(", lev.fac1[2],")"),

              paste0("M4-", fac.names[2], ".(", lev.fac2[1],")"), paste0("M4-", fac.names[2], ".(", lev.fac2[2],")"),

              paste0("M5-", fac.names[1], ".(", lev.fac1[1],")", "+", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M5-", fac.names[1], ".(", lev.fac1[2],")", "+", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M5-", fac.names[1], ".(", lev.fac1[1],")", "+", fac.names[2], ".(", lev.fac2[2],")"),
              paste0("M5-", fac.names[1], ".(", lev.fac1[2],")", "+", fac.names[2], ".(", lev.fac2[2],")"),

              paste0("M6-", fac.names[1], ".(", lev.fac1[1],")", "*", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[2],")", "*", fac.names[2], ".(", lev.fac2[1],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[1],")", "*", fac.names[2], ".(", lev.fac2[2],")"),
              paste0("M6-", fac.names[1], ".(", lev.fac1[2],")", "*", fac.names[2], ".(", lev.fac2[2],")"),

              paste0("M6-", fac.names[2], ".(", lev.fac2[1],")", "*", fac.names[1], ".(", lev.fac1[1],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[2],")", "*", fac.names[1], ".(", lev.fac1[1],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[1],")", "*", fac.names[1], ".(", lev.fac1[2],")"),
              paste0("M6-", fac.names[2], ".(", lev.fac2[2],")", "*", fac.names[1], ".(", lev.fac1[2],")")
    ),

    "DIC"=c(data.models[[1]]$DIC$dic, # ResMod
            data.models[[2]]$DIC$dic, # ResMod
            data.models[[3]]$DIC$dic, data.models[[4]]$DIC$dic, data.models[[5]]$DIC$dic, data.models[[6]]$DIC$dic, # M2
            data.models[[7]]$DIC$dic, data.models[[8]]$DIC$dic, # M3
            data.models[[9]]$DIC$dic, data.models[[10]]$DIC$dic, # M4
            data.models[[11]]$DIC$dic, data.models[[12]]$DIC$dic, data.models[[13]]$DIC$dic, data.models[[14]]$DIC$dic, # M5
            data.models[[15]]$DIC$dic, data.models[[16]]$DIC$dic, data.models[[17]]$DIC$dic, data.models[[18]]$DIC$dic, # M6
            data.models[[19]]$DIC$dic, data.models[[20]]$DIC$dic, data.models[[21]]$DIC$dic, data.models[[22]]$DIC$dic # M6
    ),

    "WAIC"=c(data.models[[1]]$WAIC$waic, # ResMod
             data.models[[2]]$WAIC$waic, # ResMod
             data.models[[3]]$WAIC$waic, data.models[[4]]$WAIC$waic, data.models[[5]]$WAIC$waic, data.models[[6]]$WAIC$waic, # M2
             data.models[[7]]$WAIC$waic, data.models[[8]]$WAIC$waic, # M3
             data.models[[9]]$WAIC$waic, data.models[[10]]$WAIC$waic, # M4
             data.models[[11]]$WAIC$waic, data.models[[12]]$WAIC$waic, data.models[[13]]$WAIC$waic, data.models[[14]]$WAIC$waic, # M5
             data.models[[15]]$WAIC$waic, data.models[[16]]$WAIC$waic, data.models[[17]]$WAIC$waic, data.models[[18]]$WAIC$waic, # M6
             data.models[[19]]$WAIC$waic, data.models[[20]]$WAIC$waic, data.models[[21]]$WAIC$waic, data.models[[22]]$WAIC$waic # M6
    ),

    "CPU"=c(data.models[[1]]$CPU[4], # ResMod
            data.models[[2]]$CPU[4], # ResMod
            data.models[[3]]$CPU[4], data.models[[4]]$CPU[4], data.models[[5]]$CPU[4], data.models[[6]]$CPU[4], # M2
            data.models[[7]]$CPU[4], data.models[[8]]$CPU[4], # M3
            data.models[[9]]$CPU[4], data.models[[10]]$CPU[4], # M4
            data.models[[11]]$CPU[4], data.models[[12]]$CPU[4], data.models[[13]]$CPU[4], data.models[[14]]$CPU[4], # M5
            data.models[[15]]$CPU[4], data.models[[16]]$CPU[4], data.models[[17]]$CPU[4], data.models[[18]]$CPU[4], # M6
            data.models[[19]]$CPU[4], data.models[[20]]$CPU[4], data.models[[21]]$CPU[4], data.models[[22]]$CPU[4] # M6
    ),

    "LOOCV"=c(data.models[[1]]$loocv,
              data.models[[2]]$loocv,
              data.models[[3]]$loocv, data.models[[4]]$loocv, data.models[[5]]$loocv, data.models[[6]]$loocv,
              data.models[[7]]$loocv, data.models[[8]]$loocv,
              data.models[[9]]$loocv, data.models[[10]]$loocv,
              data.models[[11]]$loocv, data.models[[12]]$loocv, data.models[[13]]$loocv, data.models[[14]]$loocv,
              data.models[[15]]$loocv, data.models[[16]]$loocv, data.models[[17]]$loocv, data.models[[18]]$loocv,
              data.models[[19]]$loocv, data.models[[20]]$loocv, data.models[[21]]$loocv, data.models[[22]]$loocv
    ),

    "LOGCV-3"=c(data.models[[1]]$lgocv.m3,
                data.models[[2]]$lgocv.m3,
                data.models[[3]]$lgocv.m3, data.models[[4]]$lgocv.m3, data.models[[5]]$lgocv.m3, data.models[[6]]$lgocv.m3,
                data.models[[7]]$lgocv.m3, data.models[[8]]$lgocv.m3,
                data.models[[9]]$lgocv.m3, data.models[[10]]$lgocv.m3,
                data.models[[11]]$lgocv.m3, data.models[[12]]$lgocv.m3, data.models[[13]]$lgocv.m3, data.models[[14]]$lgocv.m3,
                data.models[[15]]$lgocv.m3, data.models[[16]]$lgocv.m3, data.models[[17]]$lgocv.m3, data.models[[18]]$lgocv.m3,
                data.models[[19]]$lgocv.m3, data.models[[20]]$lgocv.m3, data.models[[21]]$lgocv.m3, data.models[[22]]$lgocv.m3
    ),

    "LOGCV-5"=c(data.models[[1]]$lgocv.m5,
                data.models[[2]]$lgocv.m5,
                data.models[[3]]$lgocv.m5, data.models[[4]]$lgocv.m5, data.models[[5]]$lgocv.m5, data.models[[6]]$lgocv.m5,
                data.models[[7]]$lgocv.m5, data.models[[8]]$lgocv.m5,
                data.models[[9]]$lgocv.m5, data.models[[10]]$lgocv.m5,
                data.models[[11]]$lgocv.m5, data.models[[12]]$lgocv.m5, data.models[[13]]$lgocv.m5, data.models[[14]]$lgocv.m5,
                data.models[[15]]$lgocv.m5, data.models[[16]]$lgocv.m5, data.models[[17]]$lgocv.m5, data.models[[18]]$lgocv.m5,
                data.models[[19]]$lgocv.m5, data.models[[20]]$lgocv.m5, data.models[[21]]$lgocv.m5, data.models[[22]]$lgocv.m5
    ),

    "LOGCV-10"=c(data.models[[1]]$lgocv.m10,
                 data.models[[2]]$lgocv.m10,
                 data.models[[3]]$lgocv.m10, data.models[[4]]$lgocv.m10, data.models[[5]]$lgocv.m10, data.models[[6]]$lgocv.m10,
                 data.models[[7]]$lgocv.m10, data.models[[8]]$lgocv.m10,
                 data.models[[9]]$lgocv.m10, data.models[[10]]$lgocv.m10,
                 data.models[[11]]$lgocv.m10, data.models[[12]]$lgocv.m10, data.models[[13]]$lgocv.m10, data.models[[14]]$lgocv.m10,
                 data.models[[15]]$lgocv.m10, data.models[[16]]$lgocv.m10, data.models[[17]]$lgocv.m10, data.models[[18]]$lgocv.m10,
                 data.models[[19]]$lgocv.m10, data.models[[20]]$lgocv.m10, data.models[[21]]$lgocv.m10, data.models[[22]]$lgocv.m10
    )
  )

  data.models[[23]] <- data_summary
  names(data.models)[23] <- "Summary"
  data.models <- INLA.SocialEp::inla.null.sp(data.models)
  print(data_summary)

  # Devolvemos lista con todos los modelos
  return(data.models)
}
