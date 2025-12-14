#' Extract null effects from SR-ANOVA models
#'
#' This function extracts the total number of null spatial effects inside a 
#' SR-ANOVA object, for each one of the different specifications included.
#'
#' @param mod Mod file coming from a SR-ANOVA object
#' @param thres Threshold of sd use for the spatial effect to consider it null
#' @return SR-ANOVA object with the summary table updated to include null spatial effects
#' @export

inla.null.sp <- function(mod, thres = 0.125) {
  # helper: for a given index and a vector of random-field names,
  # return a numeric vector of sds (or NA when not available)
  get_sds <- function(idx, fields) {
    # guard: mod[[idx]] may not exist or may not be a list
    if (length(mod) < idx) return(rep(NA_real_, length(fields)))
    obj <- mod[[idx]]
    if (!is.list(obj)) return(rep(NA_real_, length(fields)))
    sr <- obj$summary.random
    if (is.null(sr) || !is.list(sr)) return(rep(NA_real_, length(fields)))
    sapply(fields, function(f) {
      # check that sr[[f]] exists and has a "mean" column
      if (!is.null(sr[[f]]) && is.data.frame(sr[[f]]) && "mean" %in% names(sr[[f]])) {
        # safe sd (if single value or all NA -> sd gives NA, fine)
        tryCatch(sd(sr[[f]]$mean, na.rm = TRUE), error = function(e) NA_real_)
      } else {
        NA_real_
      }
    })
  }
  
  summary_temp <- mod$Summary
  
  # compute sds for each model block (NA where field not present)
  M_2  <- get_sds(2,  c("phi_1",  "phi_2",  "phi_3",  "phi_4"))
  M_3  <- get_sds(3,  c("phi_11"))
  M_4  <- get_sds(4,  c("phi_11"))
  M_5  <- get_sds(5,  c("phi_11"))
  M_6  <- get_sds(6,  c("phi_11"))
  M_7  <- get_sds(7,  c("phi_11", "phi_12"))
  M_8  <- get_sds(8,  c("phi_11", "phi_12"))
  M_9  <- get_sds(9,  c("phi_11", "phi_21"))
  M_10 <- get_sds(10, c("phi_11", "phi_21"))
  M_11 <- get_sds(11, c("phi_11", "phi_12", "phi_21"))
  M_12 <- get_sds(12, c("phi_11", "phi_12", "phi_21"))
  M_13 <- get_sds(13, c("phi_11", "phi_12", "phi_21"))
  M_14 <- get_sds(14, c("phi_11", "phi_12", "phi_21"))
  M_15 <- get_sds(15, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_16 <- get_sds(16, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_17 <- get_sds(17, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_18 <- get_sds(18, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_19 <- get_sds(19, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_20 <- get_sds(20, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_21 <- get_sds(21, c("phi_11", "phi_12", "phi_21", "phi_22"))
  M_22 <- get_sds(22, c("phi_11", "phi_12", "phi_21", "phi_22"))
  
  # build n_null: first element 0 (as in your original code), then counts
  n_null <- c(
    0,
    sum(M_2  <= thres, na.rm = TRUE),
    sum(M_3  <= thres, na.rm = TRUE),
    sum(M_4  <= thres, na.rm = TRUE),
    sum(M_5  <= thres, na.rm = TRUE),
    sum(M_6  <= thres, na.rm = TRUE),
    sum(M_7  <= thres, na.rm = TRUE),
    sum(M_8  <= thres, na.rm = TRUE),
    sum(M_9  <= thres, na.rm = TRUE),
    sum(M_10 <= thres, na.rm = TRUE),
    sum(M_11 <= thres, na.rm = TRUE),
    sum(M_12 <= thres, na.rm = TRUE),
    sum(M_13 <= thres, na.rm = TRUE),
    sum(M_14 <= thres, na.rm = TRUE),
    sum(M_15 <= thres, na.rm = TRUE),
    sum(M_16 <= thres, na.rm = TRUE),
    sum(M_17 <= thres, na.rm = TRUE),
    sum(M_18 <= thres, na.rm = TRUE),
    sum(M_19 <= thres, na.rm = TRUE),
    sum(M_20 <= thres, na.rm = TRUE),
    sum(M_21 <= thres, na.rm = TRUE),
    sum(M_22 <= thres, na.rm = TRUE)
  )
  
  summary_temp$sp.null <- n_null
  mod$Summary <- summary_temp
  return(mod)
}
