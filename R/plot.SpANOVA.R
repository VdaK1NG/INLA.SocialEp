#' Represent different values from SpANOVA modelling
#'
#' This function represents the desired value from a SpANOVA object or an INLA object after running inla.rerun.SpANOVA
#'
#' @param obj Object to generate the plot from
#' @param obj_type Type of object provided, options are SpANOVA or INLA after running inla.rerun.SpANOVA
#' @param fill_by Values to represent, choosing between Spatial, Heterogeneity and RR
#' @param n_mod Number of the model specification that the user wants to represent
#' @param sp_object Spatial object to be plotted
#' @param breaks Breaks for the color palette in case the user wants to modify the default ones
#' @param fil_scale Vector of colors in case the user wants to modify the default ones
#' @param col_frontiers Colour for the lines that define each polygon, default is black
#' @param scale_name Name for the filling scale of the plot
#' @param sp_null Threshold desired to considered a spatial effect represented null, default is 0.125
#' @param legend.position Position for the legend of the plot
#' @param ncol_fig Number of columns desired for the grid of figures
#' @return Representation of the values desired using ggplot2 as baseline
#' @export

plot.SpANOVA <- function(
    obj,
    obj_type=c("SpANOVA", "INLA"),
    fill_by=c("Spatial", "Heterogeneity", "RR"),
    n_mod,
    sp_object,
    breaks=NA,
    fil_scale=c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B"),
    col_frontiers="black",
    scale_name="Values",
    sp_null=0.125,
    legend.position="right",
    ncol_fig=2
    ){

  # Basic checks
  if(!(fill_by %in% c("Spatial", "Heterogeneity", "RR"))){stop("Please specify a proper fill_by argument.")}
  if(!(obj_type %in% c("SpANOVA", "INLA"))){stop("Please specify a proper object type.")}
  if(length(fil_scale)<2){stop("Please specify at least two colors for the colour scale.")}

  # Prepare breaks
  if(sum(is.na(breaks)) & fill_by %in% c("Spatial", "Heterogeneity")){
    breaks <- c(-5, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 5)
  }else if(sum(is.na(breaks)) & fill_by=="RR"){
    breaks <- c(0, 0.5, 0.9, 1.1, 2, Inf)
  }

  # Prepare colour scale
  CScale <- grDevices::colorRampPalette(fil_scale)

  # Retrieve groups
  groups_names <- obj[[n_mod]]$groups
  sp_names <- obj[[n_mod]]$sp_effects

  # Prepare data to plot
  if(obj_type=="SpANOVA"){
    if(fill_by=="Heterogeneity"){
      sp_object$sp1 <- obj[[n_mod]]$summary.random$omega_j$mean[1:nrow(sp_object)]
      sp_object$sp2 <- obj[[n_mod]]$summary.random$omega_j$mean[(nrow(sp_object)+1):(2*nrow(sp_object))]
      sp_object$sp3 <- obj[[n_mod]]$summary.random$omega_j$mean[((2*nrow(sp_object)+1)):(3*nrow(sp_object))]
      sp_object$sp4 <- obj[[n_mod]]$summary.random$omega_j$mean[(3*nrow(sp_object)+1):(4*nrow(sp_object))]
    }else if(fill_by=="RR"){
      sp_object$sp1 <- obj[[n_mod]]$summary.fitted.values$mean[1:nrow(sp_object)]
      sp_object$sp2 <- obj[[n_mod]]$summary.fitted.values$mean[(nrow(sp_object)+1):(2*nrow(sp_object))]
      sp_object$sp3 <- obj[[n_mod]]$summary.fitted.values$mean[((2*nrow(sp_object)+1)):(3*nrow(sp_object))]
      sp_object$sp4 <- obj[[n_mod]]$summary.fitted.values$mean[(3*nrow(sp_object)+1):(4*nrow(sp_object))]
    }else if(fill_by=="Spatial"){
      if(n_mod == 1){
        sp_object$sp1 <- NA
        sp_object$sp2 <- NA
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod == 2){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_1$mean
        sp_object$sp2 <- obj[[n_mod]]$summary.random$phi_2$mean
        sp_object$sp3 <- obj[[n_mod]]$summary.random$phi_3$mean
        sp_object$sp4 <- obj[[n_mod]]$summary.random$phi_4$mean
      }else if(n_mod %in% c(3, 4, 5, 6)){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_11$mean
        sp_object$sp2 <- NA
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(7, 8)){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_11$mean
        sp_object$sp2 <- obj[[n_mod]]$summary.random$phi_21$mean
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(9, 10)){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_11$mean
        sp_object$sp2 <- obj[[n_mod]]$summary.random$phi_12$mean
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(11, 12, 13, 14)){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_11$mean
        sp_object$sp2 <- obj[[n_mod]]$summary.random$phi_12$mean
        sp_object$sp3 <- obj[[n_mod]]$summary.random$phi_21$mean
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(15, 16, 17, 18, 19, 20, 21, 22)){
        sp_object$sp1 <- obj[[n_mod]]$summary.random$phi_11$mean
        sp_object$sp2 <- obj[[n_mod]]$summary.random$phi_12$mean
        sp_object$sp3 <- obj[[n_mod]]$summary.random$phi_21$mean
        sp_object$sp4 <- obj[[n_mod]]$summary.random$phi_22$mean
      }
    }
  }else if(obj_type=="INLA"){
    if(fill_by=="Heterogeneity"){
      sp_object$sp1 <- obj$summary.random$omega_j$mean[1:nrow(sp_object)]
      sp_object$sp2 <- obj$summary.random$omega_j$mean[(nrow(sp_object)+1):(2*nrow(sp_object))]
      sp_object$sp3 <- obj$summary.random$omega_j$mean[((2*nrow(sp_object)+1)):(3*nrow(sp_object))]
      sp_object$sp4 <- obj$summary.random$omega_j$mean[(3*nrow(sp_object)+1):(4*nrow(sp_object))]
    }else if(fill_by=="RR"){
      sp_object$sp1 <- obj$summary.fitted.values$mean[1:nrow(sp_object)]
      sp_object$sp2 <- obj$summary.fitted.values$mean[(nrow(sp_object)+1):(2*nrow(sp_object))]
      sp_object$sp3 <- obj$summary.fitted.values$mean[((2*nrow(sp_object)+1)):(3*nrow(sp_object))]
      sp_object$sp4 <- obj$summary.fitted.values$mean[(3*nrow(sp_object)+1):(4*nrow(sp_object))]
    }else if(fill_by=="Spatial"){
      if(n_mod == 1){
        sp_object$sp1 <- NA
        sp_object$sp2 <- NA
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod == 2){
        sp_object$sp1 <- obj$summary.random$phi_1$mean
        sp_object$sp2 <- obj$summary.random$phi_2$mean
        sp_object$sp3 <- obj$summary.random$phi_3$mean
        sp_object$sp4 <- obj$summary.random$phi_4$mean
      }else if(n_mod %in% c(3, 4, 5, 6)){
        sp_object$sp1 <- obj$summary.random$phi_11$mean
        sp_object$sp2 <- NA
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(7, 8)){
        sp_object$sp1 <- obj$summary.random$phi_11$mean
        sp_object$sp2 <- obj$summary.random$phi_21$mean
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(9, 10)){
        sp_object$sp1 <- obj$summary.random$phi_11$mean
        sp_object$sp2 <- obj$summary.random$phi_12$mean
        sp_object$sp3 <- NA
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(11, 12, 13, 14)){
        sp_object$sp1 <- obj$summary.random$phi_11$mean
        sp_object$sp2 <- obj$summary.random$phi_12$mean
        sp_object$sp3 <- obj$summary.random$phi_21$mean
        sp_object$sp4 <- NA
      }else if(n_mod %in% c(15, 16, 17, 18, 19, 20, 21, 22)){
        sp_object$sp1 <- obj$summary.random$phi_11$mean
        sp_object$sp2 <- obj$summary.random$phi_12$mean
        sp_object$sp3 <- obj$summary.random$phi_21$mean
        sp_object$sp4 <- obj$summary.random$phi_22$mean
      }
    }
  }

  # Estimate standard deviation for each effect
  sd_sp1 <- round(sd(pull(sp_object, sp1)), 2)
  sd_sp2 <- round(sd(pull(sp_object, sp2)), 2)
  sd_sp3 <- round(sd(pull(sp_object, sp3)), 2)
  sd_sp4 <- round(sd(pull(sp_object, sp4)), 2)

  # Create titles
  if(fill_by=="Spatial"){
    if(n_mod == 1){
      title1 <- "No Effect"
      title2 <- "No Effect"
      title3 <- "No Effect"
      title4 <- "No Effect"
    }else if(n_mod == 2){
      title1 <- bquote("Spatial Effect | "~phi[1])
      title2 <- bquote("Spatial Effect | "~phi[2])
      title3 <- bquote("Spatial Effect | "~phi[3])
      title4 <- bquote("Spatial Effect | "~phi[4])
    }else if(n_mod %in% c(3, 4, 5, 6)){
      title1 <- bquote("Spatial Effect | "~phi[11])
      title2 <- "No Effect"
      title3 <- "No Effect"
      title4 <- "No Effect"
    }else if(n_mod %in% c(7, 8)){
      title1 <- bquote("Spatial Effect | "~phi[11])
      title2 <- bquote("Spatial Effect | "~phi[12])
      title3 <- "No Effect"
      title4 <- "No Effect"
    }else if(n_mod %in% c(9, 10)){
      title1 <- bquote("Spatial Effect | "~phi[11])
      title2 <- bquote("Spatial Effect | "~phi[21])
      title3 <- "No Effect"
      title4 <- "No Effect"
    }else if(n_mod %in% c(11, 12, 13, 14)){
      title1 <- bquote("Spatial Effect | "~phi[11])
      title2 <- bquote("Spatial Effect | "~phi[12])
      title3 <- bquote("Spatial Effect | "~phi[21])
      title4 <- "No Effect"
    }else if(n_mod %in% c(15, 16, 17, 18, 19, 20, 21, 22)){
      title1 <- bquote("Spatial Effect | "~phi[11])
      title2 <- bquote("Spatial Effect | "~phi[12])
      title3 <- bquote("Spatial Effect | "~phi[21])
      title4 <- bquote("Spatial Effect | "~phi[22])
    }
  }else if(fill_by=="RR"){
    title1 <- paste0("Relative Risk Adjusted | G1")
    title2 <- paste0("Relative Risk Adjusted | G2")
    title3 <- paste0("Relative Risk Adjusted | G3")
    title4 <- paste0("Relative Risk Adjusted | G4")
  }else if(fill_by=="Heterogeneity"){
    title1 <- bquote("Heterogeneity Effect | "~omega[1])
    title2 <- bquote("Heterogeneity Effect | "~omega[2])
    title3 <- bquote("Heterogeneity Effect | "~omega[3])
    title4 <- bquote("Heterogeneity Effect | "~omega[4])
  }

  # Prepare Figures
  if(fill_by %in% c("Spatial", "Heterogeneity")){
    # Create fig1
    if(sum(is.na(sp_object$sp1))==0){
      fig_fig1 <- convert_col(
        data = sp_object$sp1,
        breaks = breaks,
        pal_fun = CScale,
        right = TRUE,
        include.lowest = FALSE,
        na.col = NA
      )

      fig_values1 <- fig_fig1$colors[which(fig_fig1$tags %in% cut(pull(sp_object, sp1),  breaks=breaks))]

      fig1 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=cut(sp1,  breaks=breaks)), colour=col_frontiers) +
        geom_sf(aes(), fill=fig_fig1$fill_by, colour=col_frontiers) +
        ggtitle(title1, subtitle = ifelse(sd_sp1 < sp_null, paste0("Standard Deviation = ", sd_sp1, "*"), paste0("Standard Deviation = ", sd_sp1))) +
        scale_fill_manual(values = fig_values1, name=scale_name) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= ifelse(sd_sp1 < sp_null, "red", "grey13"), hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }else{
      fig1 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=sp1), colour=col_frontiers) +
        scale_fill_manual(name=scale_name) +
        ggtitle(title1, subtitle = paste0("No Effect")) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= "red", hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }


    # Create fig2
    if(sum(is.na(sp_object$sp2))==0){
      fig_fig2 <- convert_col(
        data = sp_object$sp2,
        breaks = breaks,
        pal_fun = CScale,
        right = TRUE,
        include.lowest = FALSE,
        na.col = NA
      )

      fig_values2 <- fig_fig2$colors[which(fig_fig2$tags %in% cut(pull(sp_object, sp2),  breaks=breaks))]

      fig2 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=cut(sp2,  breaks=breaks)), colour=col_frontiers) +
        geom_sf(aes(), fill=fig_fig2$fill_by, colour=col_frontiers) +
        ggtitle(title2, subtitle = ifelse(sd_sp2 < sp_null, paste0("Standard Deviation = ", sd_sp2, "*"), paste0("Standard Deviation = ", sd_sp2))) +
        scale_fill_manual(values = fig_values2, name=scale_name) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= ifelse(sd_sp2 < sp_null, "red", "grey13"), hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }else{
      fig2 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=sp2), colour=col_frontiers) +
        scale_fill_manual(name=scale_name) +
        ggtitle(title2, subtitle = paste0("No Effect")) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= "red", hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }

    # Create fig3
    if(sum(is.na(sp_object$sp3))==0){
      fig_fig3 <- convert_col(
        data = sp_object$sp3,
        breaks = breaks,
        pal_fun = CScale,
        right = TRUE,
        include.lowest = FALSE,
        na.col = NA
      )

      fig_values3 <- fig_fig3$colors[which(fig_fig3$tags %in% cut(pull(sp_object, sp3),  breaks=breaks))]

      fig3 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=cut(sp3,  breaks=breaks)), colour=col_frontiers) +
        geom_sf(aes(), fill=fig_fig3$fill_by, colour=col_frontiers) +
        ggtitle(title3, subtitle = ifelse(sd_sp3 < sp_null, paste0("Standard Deviation = ", sd_sp3, "*"), paste0("Standard Deviation = ", sd_sp3))) +
        scale_fill_manual(values = fig_values3, name=scale_name) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= ifelse(sd_sp3 < sp_null, "red", "grey13"), hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }else{
      fig3 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=sp3), colour=col_frontiers) +
        scale_fill_manual(name=scale_name) +
        ggtitle(title3, subtitle = paste0("No Effect")) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= "red", hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }

    # Create fig4
    if(sum(is.na(sp_object$sp4))==0){
      fig_fig4 <- convert_col(
        data = sp_object$sp4,
        breaks = breaks,
        pal_fun = CScale,
        right = TRUE,
        include.lowest = FALSE,
        na.col = NA
      )

      fig_values4 <- fig_fig4$colors[which(fig_fig4$tags %in% cut(pull(sp_object, sp4),  breaks=breaks))]

      fig4 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=cut(sp4,  breaks=breaks)), colour=col_frontiers) +
        geom_sf(aes(), fill=fig_fig4$fill_by, colour=col_frontiers) +
        ggtitle(title4, subtitle = ifelse(sd_sp4 < sp_null, paste0("Standard Deviation = ", sd_sp4, "*"), paste0("Standard Deviation = ", sd_sp4))) +
        scale_fill_manual(values = fig_values4, name=scale_name) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= ifelse(sd_sp4 < sp_null, "red", "grey13"), hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }else{
      fig4 <- ggplot2::ggplot(data=sp_object) +
        geom_sf(aes(fill=sp4), colour=col_frontiers) +
        scale_fill_manual(name=scale_name) +
        ggtitle(title4, subtitle = paste0("No Effect")) +
        theme(
          plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
          plot.subtitle = element_text(size=11, face= "bold", colour= "red", hjust = 0),
          axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.75, "cm"),
          strip.text.x = element_text(size = 8, face = "bold.italic"),
          legend.position = legend.position
        )
    }
  }else if(fill_by == "RR"){
    # Create fig1
    fig_fig1 <- convert_col(
      data = sp_object$sp1,
      breaks = breaks,
      pal_fun = CScale,
      right = TRUE,
      include.lowest = FALSE,
      na.col = NA
    )

    fig_values1 <- fig_fig1$colors[which(fig_fig1$tags %in% cut(pull(sp_object, sp1),  breaks=breaks))]

    fig1 <- ggplot2::ggplot(data=sp_object) +
      geom_sf(aes(fill=cut(sp1,  breaks=breaks)), colour=col_frontiers) +
      geom_sf(aes(), fill=fig_fig1$fill_by, colour=col_frontiers) +
      ggtitle(title1, subtitle = groups_names[1]) +
      scale_fill_manual(values = fig_values1, name=scale_name) +
      theme(
        plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
        plot.subtitle = element_text(size=11, face= "bold", colour= "grey13", hjust = 0),
        axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.75, "cm"),
        strip.text.x = element_text(size = 8, face = "bold.italic"),
        legend.position = legend.position
      )

    # Create fig2
    fig_fig2 <- convert_col(
      data = sp_object$sp2,
      breaks = breaks,
      pal_fun = CScale,
      right = TRUE,
      include.lowest = FALSE,
      na.col = NA
    )

    fig_values2 <- fig_fig2$colors[which(fig_fig2$tags %in% cut(pull(sp_object, sp2),  breaks=breaks))]

    fig2 <- ggplot2::ggplot(data=sp_object) +
      geom_sf(aes(fill=cut(sp2,  breaks=breaks)), colour=col_frontiers) +
      geom_sf(aes(), fill=fig_fig2$fill_by, colour=col_frontiers) +
      ggtitle(title2, subtitle = groups_names[2]) +
      scale_fill_manual(values = fig_values2, name=scale_name) +
      theme(
        plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
        plot.subtitle = element_text(size=11, face= "bold", colour= "grey13", hjust = 0),
        axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.75, "cm"),
        strip.text.x = element_text(size = 8, face = "bold.italic"),
        legend.position = legend.position
      )

    # Create fig3
    fig_fig3 <- convert_col(
      data = sp_object$sp3,
      breaks = breaks,
      pal_fun = CScale,
      right = TRUE,
      include.lowest = FALSE,
      na.col = NA
    )

    fig_values3 <- fig_fig3$colors[which(fig_fig3$tags %in% cut(pull(sp_object, sp3),  breaks=breaks))]

    fig3 <- ggplot2::ggplot(data=sp_object) +
      geom_sf(aes(fill=cut(sp3,  breaks=breaks)), colour=col_frontiers) +
      geom_sf(aes(), fill=fig_fig3$fill_by, colour=col_frontiers) +
      ggtitle(title3, subtitle = groups_names[3]) +
      scale_fill_manual(values = fig_values3, name=scale_name) +
      theme(
        plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
        plot.subtitle = element_text(size=11, face= "bold", colour= "grey13", hjust = 0),
        axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.75, "cm"),
        strip.text.x = element_text(size = 8, face = "bold.italic"),
        legend.position = legend.position
      )

    # Create fig4
    fig_fig4 <- convert_col(
      data = sp_object$sp4,
      breaks = breaks,
      pal_fun = CScale,
      right = TRUE,
      include.lowest = FALSE,
      na.col = NA
    )

    fig_values4 <- fig_fig4$colors[which(fig_fig4$tags %in% cut(pull(sp_object, sp4),  breaks=breaks))]

    fig4 <- ggplot2::ggplot(data=sp_object) +
      geom_sf(aes(fill=cut(sp4,  breaks=breaks)), colour=col_frontiers) +
      geom_sf(aes(), fill=fig_fig4$fill_by, colour=col_frontiers) +
      ggtitle(title4, subtitle = groups_names[4]) +
      scale_fill_manual(values = fig_values4, name=scale_name) +
      theme(
        plot.title = element_text(size=14, face= "bold", colour= "black", hjust = 0),
        plot.subtitle = element_text(size=11, face= "bold", colour= "grey13", hjust = 0),
        axis.text.x=element_blank(),  axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.75, "cm"),
        strip.text.x = element_text(size = 8, face = "bold.italic"),
        legend.position = legend.position
      )
  }

  # Create final figure
  gridExtra::grid.arrange(fig1, fig2, fig3, fig4, ncol=ncol_fig)
}
