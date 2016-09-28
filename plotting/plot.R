library(ggplot2)
library(R.matlab)
library(gridExtra)
library(reshape2)
library(scales)

#setwd("~/Dropbox/seunglab/plotting_retina_figures/")
setwd("~/../Dropbox/seunglab/plotting_retina_figures/")


#=====================================================================
#Nice top-level variables to have
all_gc_types <- c(
  "1ni","1no","1ws","1wt",
  "25","27","28","2aw","2an","2i","2o",
  "37c","37d","37r","37v","3i","3o",
  "4i","4on","4ow",
  "51","5si","5so","5ti","5to",
  "63","6sn","6sw","6t",
  "72","73","7ir","7id","7iv","7o",
  "81i","81o","82n","82wi","82wo","85","8n","8w",
  "91","915","9n","9w")

alpha_type_list = c("1wt","4ow","8w","6sw")
mini_alpha_candidate_and_compared_list <- c("4ow","4on","6sw","6sn","4i")

#Strat profile binning param
default_binwidth = 3

#=====================================================================
source("data_IO.R")
source("id_utils.R")
source("data_manipulation.R")
source("aesthetics.R")
source("significance_tests.R")
#=====================================================================


plot_soma_size_hist <- function(type_l){
  
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  vols <- soma_sizes_for_ids( ids )
  
  d <- data.frame("vol"=vols, "cell_id"=ids, "type"=types)
  
  d$alphatype <- 'Other'
  #d$alphatype[d$type == "1wt"] = "1wt,4ow,8w"
  #d$alphatype[d$type == "4ow"] = "1wt,4ow,8w"
  #d$alphatype[d$type == "8w"] = "1wt,4ow,8w"
  d$alphatype[d$type == "1wt"] = "Classic Alpha"
  d$alphatype[d$type == "4ow"] = "Classic Alpha"
  d$alphatype[d$type == "8w"] = "Classic Alpha"
  d$alphatype <- as.factor(d$alphatype)
  
  qplot( d$vol, geom="histogram", fill=d$alphatype) +
    labs(x = expression(paste( "Soma Volume (\u03BC","m"^"3",")")),
         y = "# Cells") +
    scale_fill_manual(breaks="Classic Alpha",
                      values=c("red", color_palette_iwh_somahist[6])) +
    scale_x_continuous(limits=c(0,3000),expand=c(0,0), breaks=seq(0,2500,by=500)) +
    scale_y_continuous(limits=c(0,90),expand=c(0,0)) +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA),
          legend.position = c(1,1),
          legend.justification=c(1,1),
          legend.title = element_blank(),
          text = element_text(size=30, family="Arial"),
          axis.text = element_text(size=30)
    )
}

plot_arbor_vectors <- function( type_l, cell_info, id_mapping, asym, asym_ids, mult_factor ){
  
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  soma_centroids <- soma_proj_for_ids( ids, cell_info, id_mapping, warped=T  )
  asym_vectors <- arbor_vectors_for_ids( ids, asym, asym_ids) 
  
  asym_vectors <- asym_vectors * mult_factor
  
  d <- data.frame("cell_id"=ids, "cell_type"=types,
                  "y"=soma_centroids[,1], "z"=soma_centroids[,2], 
                  "yend"=soma_centroids[,1] + asym_vectors[,1],
                  "zend"=soma_centroids[,2] + asym_vectors[,2])
  
  ymax = 5376
  zmax = 3456 * 23 / 16.5
  
  d$cell_id <- as.factor(d$cell_id)
  
  #flipping y the old-fashioned way
  d$z    <- (zmax) - d$z
  d$zend <- (zmax) - d$zend
  
  #d = d[d$cell_id == 17028,]
  
  ggplot(d, aes(x=y, y=z)) + 
    theme_void() +
    coord_fixed() +
    scale_x_continuous(limits=c(0,ymax),expand=c(0,0)) +
    scale_y_continuous(limits=c(0,zmax),expand=c(0,0)) +
    theme(legend.text=element_text(size=20)) +
    geom_vline(xintercept=ymax, size=3) + 
    geom_vline(xintercept=0, size=3) +
    geom_hline(yintercept=zmax, size=3) + 
    geom_hline(yintercept=0, size=3) +
    guides(colour=F) +
    #geom_segment( aes(x=y, y=z, xend=yend, yend=zend, colour=cell_id),
    #              size=1, arrow=arrow()) +
    scale_colour_manual(values=color_palette_iwh2) +
    geom_segment( aes(x=y, y=z, xend=yend, yend=zend, colour="gray50"),
                  size=6, arrow=arrow()) +
    geom_point(size=15, colour="black") +
    geom_point(size=13, mapping=aes(colour=cell_id)) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA))
}


plot_ca_tuning_curve <- function( ids, tuning, ca_id_map ){
  
  ca_ids <- fetch_ca_ids( ids, ca_id_map )
  cat(ca_ids)
  id_curves <- extract_tuning( tuning, ca_ids )
  cat(id_curves)
  
  d <- data.frame()
  if( length(ids) > 1 ){
    for(i in 1:length(ids)){
      cell_d <- data.frame( "cell_id"=ids[i], "tuning"=id_curves[,i],
                            "theta"=stim_order )
      
      d <- rbind(d, cell_d)
    }  
  } else {
    d <- data.frame( "cell_id"=ids[1], "tuning"=id_curves,
                     "theta"=stim_order )
  }
  
  responses360 <- d[d$theta == 0,]
  responses360$theta <- 2*pi
  d <- rbind(d, responses360)
  
  d$cell_id <- as.factor(d$cell_id)
  angle_order <- order(d$theta)
  d <- d[angle_order,]
  
  #transforming angles to standard coords
  d$theta <- (d$theta + (pi/2)) %% (2*pi)
  d$x <- cos(d$theta) * d$tuning
  d$y <- sin(d$theta) * d$tuning
  
  theta <- seq(from=0,by=.01,to=2*pi)
  r = max(d$tuning)
  circle <- data.frame(x=r*cos(theta), y=r*sin(theta))
  
  ggplot(d, aes(x=x,y=y)) +
    #frame
    geom_segment(x=r,y=0,xend=-r,yend=0, colour="gray60") +
    geom_segment(x=0,y=r,xend=0,yend=-r, colour="gray60") +
    geom_path(data=circle, size=3) +
    annotate("text",label="0",x=r+1,y=0, size=5) +
    annotate("text",label="90",x=0,y=r+1, size=5) +
    annotate("text",label="180",x=-r-1,y=0, size=5) +
    annotate("text",label="270",x=0,y=-r-1, size=5) +
    #data
    geom_path(colour=color_palette_iwh[1], size=3) + 
    geom_point(color=color_palette_iwh[1], size=8) +
    #theme
    theme_void() 
}


plot_ca_trials <- function( ids, roi_sums_all, ca_id_map, num_trials_to_plot ){
  
  ca_ids <- fetch_ca_ids(ids, ca_id_map )
  id_traces <- extract_traces( roi_sums_all, ca_ids )
  
  trial_vector <- assign_trials()
  time_points <- assign_timepoints()
  
  d <- data.frame()
  
  if( length(ids) > 1 ){
    for(i in 1:length(ids)){
      cell_d <- data.frame( "cell_id"=ids[i], "Ca_signal"=id_traces[,i], 
                            "trial"=trial_vector, "time"=time_points,
                            "stimulus" = stimuli)
    
      cell_d <- remove_extra_trace_frames(cell_d)
      d <- rbind(d, cell_d)  
    }
  } else {
    d <- data.frame( "cell_id"=ids[1], "Ca_signal"=id_traces,
                     "trial"=trial_vector, "time"=time_points)
  }
  
  d$cell_id <- as.factor(d$cell_id)
  
  d <- d[d$trial <= (num_trials_to_plot + 1) & d$trial > 1,]
  d$trial   <- as.factor(d$trial)
  
  cat(max(d$Ca_signal))
  cat("\n")
  cat(min(d$Ca_signal))
  
  plot <- ggplot(d, aes(x=time,y=Ca_signal,colour=cell_id))
  for( i in 1:length(stim_times) ){
    plot <- plot + geom_vline(xintercept=stim_times[i], colour="gray40", linetype="dashed", alpha=0.8)
  }
  
  plot +
    geom_line(size=1) +
    facet_wrap( ~ trial, nrow=num_trials_to_plot) +
    guides(colour=F) +
    labs(x="Time (s)",y="delF/F") +
    scale_x_continuous(expand=c(0,0)) +
    scale_colour_manual(values= color_palette_iwh) +
    theme(panel.background = element_rect(fill = "transparent",colour = NA),
          text = element_text(size=30),
          panel.border = element_rect(fill=NA),
          text = element_text(size=30),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
}

plot_ca_trials_by_stim <- function( ids, roi_sums_all, ca_id_map, num_trials_to_plot ){
  
  ca_ids <- fetch_ca_ids(ids, ca_id_map )
  id_traces <- extract_traces( roi_sums_all, ca_ids )
  
  trial_vector <- assign_trials()
  time_points <- assign_timepoints()
  stimuli <- assign_stimuli()
  
  d <- data.frame()
  
  if( length(ids) > 1 ){
    for(i in 1:length(ids)){
      cell_d <- data.frame( "cell_id"=ids[i], "Ca_signal"=id_traces[,i], 
                            "trial"=trial_vector, "time"=time_points,
                            "stimulus" = stimuli)
      
      cell_d <- remove_extra_trace_frames(cell_d)
      d <- rbind(d, cell_d)  
    }
  } else {
    d <- data.frame( "cell_id"=ids[1], "Ca_signal"=id_traces,
                     "trial"=trial_vector, "time"=time_points,
                     "stimulus" = stimuli)
  }
  
  d$cell_id <- as.factor(d$cell_id)
  
  d <- d[d$trial <= (num_trials_to_plot + 1) & d$trial > 1,]
  d$trial   <- as.factor(d$trial)
  
  plots <- list()
  for(i in 1:8){
    
    d_stim <- d[d$stimulus == i,]
    
    stim_on_time_index = (i-1)*2+1
    stim_off_time_index = (i-1)*2+2
    
    new_plot <- ggplot(d_stim, aes(x=time, y=Ca_signal, colour=cell_id)) + 
      geom_vline(xintercept=stim_times[stim_on_time_index],  colour="gray40", linetype="dashed", alpha=0.8) + 
      geom_vline(xintercept=stim_times[stim_off_time_index], colour="gray40", linetype="dashed", alpha=0.8) +
      geom_line(size=1) +
      facet_wrap( ~ stimulus, scales="free", nrow=1) +
      guides(colour=F) +
      labs(x="Time (s)",y="Ca") +
      scale_x_continuous(expand=c(0,0)) +
      scale_colour_manual(values= color_palette_iwh) +
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid = element_blank(),
            text = element_text(size=10),
            panel.border = element_rect(fill=NA),
            text = element_text(size=30),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      #TEMPORARY for fig1f
      scale_y_continuous( limits=c(-14.7, 25))
    cat("here")
    plots[[length(plots) + 1]] <- new_plot
  }
  plots
}


plot_id_strats <- function( ids, skel_strat, binwidth=default_binwidth ){
  
  fake_types <- rep('0',length(ids))
  cat(ids)
  
  strats <- strats_for_ids( ids, fake_types, skel_strat, binwidth )
  
  ggplot(strats, aes(x=Depth,y=Stratification,colour=cell_id)) +
    geom_vline(xintercept=0.45, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_vline(xintercept=0.28, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_vline(xintercept=0.62, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_path(size=2, alpha=0.7) +
    scale_x_continuous(breaks=seq(0,1,by=0.2), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,4.75), breaks=c(1,4), expand=c(0,0)) +
    scale_colour_manual(values= color_palette_iwh, guide=guide_legend(ncol=2)) +
    labs(colour="Cell ID", x="IPL Depth", y="Skeleton Density") +
    #theme_classic() +
    guides(colour=F) +
    theme(panel.border = element_rect(fill=NA),
          text = element_text(size=30),
          axis.text = element_text(size=30),
          legend.position = c(0,1),
          legend.justification=c(0,1),
          #legend.text = element_text(size=32),
          legend.title = element_blank())
}


plot_type_strats <- function( type_l, skel_strat, binwidth=default_binwidth ){
  type_l <- sort(type_l)
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  plot_id_strats( ids, skel_strat, binwidth )
}


plot_type_strats_from_nothing <- function( type_l ){
  skel_strat <- readMat("skel_strat.mat")
  plot_type_strats( type_l, skel_strat)
}


plot_avg_type_strats <- function( type_l, skel_strat, binwidth=default_binwidth ){
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  strats <- strats_for_ids( ids, types, skel_strat, binwidth )
  
  ggplot(strats, aes(x=Depth,y=Stratification, colour=cell_type)) +
    stat_summary(fun.y="mean", geom="line", size=2, alpha=0.8) +
    geom_vline(xintercept=0.45, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_vline(xintercept=0.28, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_vline(xintercept=0.62, size=2, colour="gray50", linetype="dashed", alpha=0.6) +
    scale_x_continuous(breaks=seq(0,1,by=0.2), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,4.75), expand=c(0,0)) +
    scale_colour_manual(values= color_palette_iwh) +#, guide=guide_legend(ncol=3)) +
    theme_classic() +
    labs(colour="Cell Type", x="IPL Depth", y="Skeleton Density") +
    theme(text = element_text(size=30, family="sans"),
          axis.text = element_text(size=30),
           panel.border = element_rect(fill=NA),
           legend.position = c(0,1),
           legend.justification=c(0,1),
           legend.text = element_text(size=30),
           legend.title = element_blank())
}


plot_id_traces <- function( ids, ca_mapping_array, ca_avgs ){
  
  types <- rep('0',length(ids))
  
  ca_ids <- fetch_ca_ids( ids, ca_mapping_array )
  
  traces <- overall_avg_ca_traces_for_ids( ca_ids, types, ca_avgs, ids )
  traces$cell_id <- as.factor(traces$cell_id)
  
  means <- tapply( traces$Trace, list(traces$Time, traces$cell_id), mean, na.rm=TRUE)
  shiftmeans <- t(means) - apply(means,2,min) #shift each min to 0
  normmeans <- melt(t( shiftmeans / apply(t(shiftmeans),2,max)))
  colnames(normmeans) <- c("Time","cell_id","Trace")
  normmeans$cell_id <- as.factor(normmeans$cell_id)
  
  ggplot(normmeans, aes(x=Time,y=Trace,colour=cell_id)) +
    geom_line(size=3) +
    scale_colour_manual(values=color_palette_iwh) +
    theme(text = element_text(size=35))+#, face='bold'))
    labs(colour="Cell Type", x="Time (s)", y="Avg Ca2+ Trace") +
    theme(legend.position = c(0,1),
          legend.justification=c(0,1),
          legend.text = element_text(size=35),
          legend.title = element_text(size=35, face="bold"))
}


plot_avg_type_traces <- function( type_l, ca_mapping_array, ca_avgs ){
  ids <- fetch_ids( type_l )
  ids <- ids[ids$ids == 20233,] #temporary for fig1f
  types <- ids$ret_types
  ids <- ids$ids
  
  ca_ids <- fetch_ca_ids( ids, ca_mapping_array )
  
  cat(ca_ids)
  cat("\n")
  cat(types)
  
  traces <- overall_avg_ca_traces_for_ids( ca_ids, types, ca_avgs, ids )
  
  #return(traces)
  
  means <- data.frame(tapply( traces$Trace, list(traces$Time, traces$cell_type), mean, na.rm=TRUE))
  shiftmeans <- t(means) - apply(means,2,min) #shift each min to 0
  normmeans <- melt(t( shiftmeans / apply(t(shiftmeans),2,max)))
  colnames(normmeans) <- c("Time","cell_type","Trace")
  normmeans$cell_type <- as.factor(normmeans$cell_type)
  
  
  
  ggplot(normmeans, aes(x=Time,y=Trace, colour=cell_type)) +
    geom_vline(xintercept=1, size=1, colour="gray60", linetype="dashed", alpha=0.6) +
    geom_vline(xintercept=2, size=1, colour="gray60", linetype="dashed", alpha=0.6) +
    geom_line(size=1) +
    #geom_abline(slope=0) +
    guides(colour=F) +
    scale_colour_manual(values= color_palette_iwh) +
    theme(text = element_text(size=35))+#, face='bold'))
    labs(colour="Cell Type", x="Time (s)", y="Normalized Avg Ca2+ Trace") +
    theme(legend.position = c(0,1),
          legend.justification=c(0,1),
          legend.text = element_text(size=55),
          legend.title = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size=30),
          panel.border = element_rect(fill=NA),
          text = element_text(size=30),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    #TEMPORARY for fig1f
    scale_y_continuous( limits=c(-0.352, 1.733))
}


plot_type_diams <- function( type_l, cell_info, id_lookup ){
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  diams <- diameters_for_ids( ids, cell_info, id_lookup )
  
  d <- data.frame(ids, types, diams)
  colnames(d) <- c("cell_id","cell_type","max_diam")
  
  ggplot(d, aes(x=cell_type, y=max_diam, colour=cell_type)) +
    geom_boxplot() +
    geom_jitter(alpha=0.5, size=3, width=0.2) +
    labs(x="Cell Type",y="Maximum Diameter (nm?)", colour="Cell Type") +
    guides(colour=FALSE)
  #qplot( d$alpha, d$radius, geom='boxplot', colour=d$type) + geom_point(alpha=0.7) +
  #  labs(title="Type vs. Computed Radius",x="Alpha Type", y="Computed Radius",legend="Cell Type")
}


plot_type_boxplot <- function( type_l , for_ids_fn, metric_name, 
                               mult_factor=NULL, savename=NULL,
                               ... ){
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  metric <- for_ids_fn( ids, ... )
  
  d <- data.frame(ids, types, metric)
  colnames(d) <- c("cell_id","cell_type","metric")
  
  if(!is.null(mult_factor)){
    d$metric <- d$metric * mult_factor
  }
  
  if(!is.null(savename)){
    write.csv(d, savename, row.names=F)
  }
  
  plot_themed_boxplot( d, metric_name)
}

plot_type_file_boxplot <- function( filename, metric_name, mult_factor=NULL ){
  
  d <- read.csv(filename)
  
  if(!is.null(mult_factor)){
    d$metric <- d$metric * mult_factor
  }
  
  plot_themed_boxplot( d, metric_name )
}


plot_themed_boxplot <- function( d, ylabel ){
  
  ggplot(d, aes(x=cell_type, y=metric, colour=cell_type)) +
    geom_boxplot() +
    geom_point( alpha = 0.5, size=3 ) +
    labs(x="Cell Type",
         y= ylabel,
         colour="Cell Type") +
    guides(colour=FALSE) +
    theme_classic() +
    #scale_y_continuous(limits=c(-1,1), expand=c(0.01,0)) + #for asym index only
    scale_y_continuous(expand=c(0.01,0)) + #for SAC contact plots
    coord_flip() +
    #geom_hline(yintercept=0, size=1, colour="gray60", linetype="dashed", alpha=0.6) +
    #scale_y_continuous(limits=c(0,500)) +
    theme(text=element_text(size=30),#, family="Arial"),
          axis.text = element_text(size=16),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))
}

plot_asym_angles <- function( type_l, gcsac, num_ids, denom_ids, uniform=T) {#cell_info, cell_info_ids, uniform=T ){
  
  ids <- fetch_ids( type_l )
  #types <- ids$ret_types
  ids <- ids$ids
  
  #angles <- dss_angles_for_ids( ids, gcsac, num_ids, denom_ids )
  #rs <- dssi_for_ids(ids, gcsac, num_ids, denom_ids )
  angles <- arbor_angles_for_ids( ids, cell_info, cell_info_ids )
  rs <- arbor_magnitudes_for_ids( ids, cell_info, cell_info_ids )
  
  if(uniform){
    rs <- rep(1,length(rs))
  }
  
  d <- data.frame("cell_id"=ids, "r"=rs, "theta"=angles)
  intercepts <- data.frame("r"=0,"theta"=d$theta,"cell_id"=ids)
  d <- rbind(d, intercepts)
  d$cell_id <- as.factor(d$cell_id)
  
  #conversion to degrees [0,360)
  d$theta <- (d$theta * (180/pi)) %% 360
  
  
  ggplot(d, aes(x=theta,y=r,colour=cell_id)) +
    geom_line( size=1, alpha=0.5, lineend="round") +
    scale_x_continuous(lim=c(0,360), breaks=round(seq(0,360,90))) +
    coord_polar(theta="x", start=-(pi/2), direction=-1) +
    #guides(colour=F) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
    theme(title=element_text(size=8))
}

plot_soma_fn_scatter <- function( for_ids_fn, y_name, mult_factor=NULL, ... ){
  
  ids <- fetch_ids( all_gc_types )
  types <- ids$ret_types
  ids <- ids$ids
  
  xs <- soma_sizes_for_ids( ids )
  ys <- for_ids_fn( ids, ... )
  
  if(!is.null(mult_factor)){
    ys <- ys * mult_factor
  }
  
  d <- data.frame( "cell_id"=ids, "cell_type"=types, "x"=xs, "y"=ys )
  
  plot_themed_scatter(d, y_name)
}

plot_soma_file_scatter <- function( y_filename, y_name, mult_factor=NULL ){
  
  ids <- fetch_ids( all_gc_types )
  types <- ids$ret_types
  ids <- ids$ids
  
  xs <- soma_sizes_for_ids( ids )
  ys <- read.csv(y_filename)$metric
  
  if(!is.null(mult_factor)){
    ys <- ys * mult_factor
  }
  
  d <- data.frame( "cell_id"=ids, "cell_type"=types, "x"=xs, "y"=ys )
  
  plot_themed_scatter(d, y_name)
}

plot_pl_corr <- function(){
  
  initial_plot <- plot_soma_file_scatter( "pl.csv", "Path Length (mm)", 1/1000 )
  
  pl_data <- read.csv("pl.csv")
  pl_data$metric <- pl_data$metric / 1000 #scaling to mm
  
  pl_examples <- pl_data[pl_data$cell_id %in% c(20118,26062),]
  pl_examples$soma <- soma_sizes_for_ids( pl_examples$cell_id )
  
  soma_examples <- pl_data[pl_data$cell_id %in% c(26079,20156),]
  soma_examples$soma <- soma_sizes_for_ids( soma_examples$cell_id )
  
  initial_plot +
    geom_text(data=pl_examples, mapping=aes(x=soma, y=metric-0.3, label=cell_id), size=5) +
    geom_text(data=soma_examples, mapping=aes(x=soma-5, y=metric-0.3, label=cell_id), size=5)
}

plot_themed_scatter <- function( d, y_name ){
  
  corr_test <- cor.test(d$x,d$y)
  cat(corr_test$estimate)
  cat("\n")
  cat(corr_test$p.value)
  cat("\n")
  
  ggplot( d, aes(x=x,y=y)) +
    geom_point(size=2, alpha=0.4) +
    geom_smooth(size=2, alpha=0.6,
                method="lm", se=F, colour=color_palette_iwh_somahist[5]) +
    geom_smooth(size=2, alpha=0.6,
                span=2, se=F, colour=color_palette_iwh_somahist[4]) +
    labs(x="Soma Size (um^3)", y=y_name) +
    scale_x_continuous(expand=c(0.01,0)) +
    scale_y_continuous(expand=c(0.01,0)) +
    theme_classic() +
    theme(text=element_text(size=30),#, family="Arial"),
          axis.text = element_text(size=16),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))
}

plot_type_median_scatter <- function( type_l, x_for_ids, y_for_ids, x_name, y_name, ... ){
  
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  #not putting ... here
  # I'm cheating bc I know that x is going to be soma size for now
  xs <- x_for_ids( ids )
  ys <- y_for_ids( ids, ... )
  
  d <- data.frame( "cell_id"=ids, "cell_type"=types, "x"=xs, "y"=ys )
  
  #x_medians <- melt(tapply( d$x, d$cell_type, median, na.rm=T ))
  #y_medians <- melt(tapply( d$y, d$cell_type, median, na.rm=T ))
  
  #colnames(x_medians) <- c("cell_type","x")
  #colnames(y_medians) <- c("cell_type","y")
  
  #d <- cbind(x_medians, y_medians$y)
  #colnames(d) <- c("cell_type","x","y")
  
  ggplot( d, aes(x=x,y=y)) +
    geom_point(size=3, alpha=0.2) +
    labs(x=x_name, y=y_name)
}


#Useful for polar selectivity plots
multiplot <- function(plots, file, cols=1, layout=NULL) {
  library(grid)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col)) +
        theme(panel.grid = element_blank())
    }
  }
}


plot_polar_selectivity <- function( type_l, dsos, dsos_ids, ds_not_os, on_not_off, normalize=F ){
  
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  if(ds_not_os){
    r <- dsos$ds.r
    theta <- dsos$ds.theta
  } else {
    r <- dsos$os.r
    theta <- dsos$os.theta
  }
  
  if(on_not_off){
    r <- r[,1]
    theta <- theta[,1]
    r_norm <- dsos$r.mean[,1]
  } else {
    r <- r[,2]
    theta <- theta[,2]
    r_norm <- dsos$r.mean[,2]
  }
  
  if(normalize){
    r <- r / r_norm
  }
  
  indices <- dsos_ids[as.character(ids)]
  d <- data.frame("r"=r, "theta"=theta, "cell_id"=dsos$omni.id )
  d$theta <- (d$theta - min(d$theta)) * 360 / (2*pi) #conversion to degrees
  d$theta <- (180 - d$theta) %% 360 #flipping to dvrc system (kinda for now)
  d <- d[indices,]
  d$types <- types
  intercepts <- data.frame("r"=0,"theta"=d$theta,"cell_id"=ids)
  intercepts$types <- types
  
  d <- rbind(d, intercepts)
  d$cell_id <- as.factor(d$cell_id)
  
  plots <- list()
  for( i in 1:length(type_l)){
    
    type_d <- d[d$types == type_l[i],]
    
    new_plot <- ggplot(type_d, aes(x=theta,y=r,colour=cell_id)) +
      geom_line( size=1, alpha=0.5, lineend="round") +
      scale_x_continuous(lim=c(0,360), breaks=round(seq(0,330,180))) +
      coord_polar(theta="x", start=0) +
      #guides(colour=F) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
      labs(title=type_l[i]) +
      theme(title=element_text(size=8))
    
    if(normalize){
      ymax <- round(max( c(type_d$r, 0.49), na.rm=T ),2) + 0.01
      new_plot <- new_plot + scale_y_continuous(limits=c(0,ymax),breaks=c(0,ymax))
    }
    
    plots[[length(plots) + 1]] <- new_plot
      
  }
  #ggplot(d, aes(x=theta,y=r,colour=cell_id)) +
    #geom_line( size=2, alpha=0.5, lineend="round") +
    #scale_x_continuous(lim=c(0,360), breaks=round(seq(0,330,30))) +
    #scale_y_continuous(lim=c(0,max_y),breaks=c()) +
    #facet_wrap( ~ types, nrow=3) +
    #coord_polar(theta="x", start=0) +
    #guides(colour=F) +
    #labs(title=title) +
    #theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  plots
}

plot_polar_arbor_vectors <- function( type_l, asym, asym_ids, normalize=F ){
  
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  r <- arbor_magnitudes_for_ids( ids, asym, asym_ids )
  theta <- arbor_angles_for_ids( ids, asym, asym_ids )
  
  d <- data.frame("r"=r, "theta"=theta, "cell_id"=ids, "cell_type"=types)
  d$theta <- d$theta * (180/pi) #conversion to degrees
  d$theta <- (-d$theta) %% 360 #flipping to dvrc system
  intercepts <- data.frame("r"=0,"theta"=d$theta,
                           "cell_id"=ids, "cell_type"=types)
  
  d <- rbind(d, intercepts)
  d$cell_id <- as.factor(d$cell_id)
  
  plots <- list()
  for( i in 1:length(type_l)){
    
    type_d <- d[d$cell_type == type_l[i],]
    
    new_plot <- ggplot(type_d, aes(x=theta,y=r,colour=cell_id)) +
      geom_line( size=1, alpha=0.5, lineend="round") +
      scale_x_continuous(lim=c(0,360), breaks=round(seq(0,330,180))) +
      coord_polar(theta="x", start=(pi/2), direction=1) +
      guides(colour=F) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
      labs(title=type_l[i]) +
      theme(title=element_text(size=8))
    
    if(normalize){
      ymax <- round(max( c(type_d$r, 0.49), na.rm=T ),2) + 0.01
      new_plot <- new_plot + scale_y_continuous(limits=c(0,ymax),breaks=c(0,ymax))
    }
    
    plots[[length(plots) + 1]] <- new_plot
    
  }
  
  plots
}


plot_type_hist_overlays <- function( for_ids_fn, type_l=all_gc_types ){
  1 #stub
}


plot_1d_id_hist <- function(metric, ids, metric_name){
  
  d <- data.frame("metric"=metric, "cell_id"=ids)
  d$cell_id <- as.factor(d$cell_id)
  
  ggplot(d, aes(x=metric, fill=cell_id)) +
    geom_histogram(bins=8) +
    geom_point(y=10, size=5, mapping=aes(colour=cell_id)) +
    labs(x=metric_name, title=sprintf("8w by %s", metric_name)) +
    ylim(c(0,10.5))  
}

#useful fn for pval plot (from StackOverflow 
#http://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#plots results of arbor orientation significance test
plot_pvals <- function( pvals, sci=T ){
  
  types <- names(pvals)
  
  d <- data.frame("types"=types, "pval"=pvals)
  d$types <- as.factor(d$types)
  
  d2 <- data.frame(x=1, y=types, label=types)
  if(sci){
    d2 <- rbind(d2, data.frame(x=2, y=types, label=format(pvals, scientific=T)))
  } else {
    d2 <- rbind(d2, data.frame(x=2, y=types, label=format(round(pvals, 3))))
  }
  
  xstart <- 0.625
  xend <- 2.625
  ystart <- 0.5
  yend <- (1+length(pvals))-0.5
  
  horizontals <- data.frame( "xstart"=xstart, "xend" = xend, "ys" = ystart:yend)
  verticals <- data.frame( "xs"=c(xstart,xend,1.375), "ystart"=ystart, "yend"=yend)
  
  pval_table <- ggplot(d2, aes(x=x, y=y)) + 
    geom_text(aes(label=label), size=5.5) +
    scale_x_continuous(limits=c(0,4)) +
    theme_void() +
    theme(panel.background = element_rect(fill=NA, colour="black", size=0.5),
          text=element_text(size=30),
          axis.text= element_text(size=16))
  
  pval_table <- pval_table + 
    geom_segment(aes(x=xstart, xend=xend, y=ys, yend=ys), data=horizontals, size=0.5) +
    geom_segment(aes(x=xs, xend=xs, y=ystart, yend=yend), data=verticals, size=0.5)
  
  #first hist
  pval_by_type <- ggplot(d, aes(types)) + 
    geom_hline(yintercept=(0.01/length(pvals)), size=1, colour="gray50", linetype="dashed", alpha=0.6) +
    geom_bar(stat="identity",aes(y=pvals, fill=types)) +
    guides(fill=F) +
    labs(x="Cell Type",y="Rayleigh Test P-Value") +
    scale_y_continuous(expand=c(0,.01), trans=reverselog_trans(10), breaks=c(1e-2,1e-4,1e-6,1e-8,1e-10)) +
    #scale_x_reverse()+
    coord_flip() +
    theme_classic() +
    theme(panel.background = element_rect(fill=NA, colour="black", size=0.5),
          text=element_text(size=30),
          axis.text = element_text(size=16) #unsure about this change
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank()
          )
  
  #pval_dist <- ggplot(d, aes(x=pvals)) +
  #  geom_histogram() +
  #  scale_x_log10() +
  #  scale_y_discrete(expand=c(0,0), breaks=NULL) +
  #  geom_blank(data=data.frame('x'=1,'y'=52), mapping=aes(x=x,y=y)) +
  #  labs(x=NULL, y=NULL, title=NULL) +
  #  #layer(geom="text",mapping=aes(x=pvals, y=types, label=types), stat="identity", position=position_nudge(y=9)) +
  #  theme_classic() +
  #  theme(panel.background=element_rect(fill=NA)) +#,
  #        #axis.text.y = element_blank(),
  #        #axis.ticks.y = element_blank()) +
  #  coord_flip()
  
  # grid.arrange(grobs=list(pval_by_type, pval_dist), nrow=1, ncol=2, width=c(10,1), heights=c(1))
  pval_by_type
  #pval_table
}

plot_complexity_comparison <- function(ids, fills){
  
  d <- read.csv('bpppl.csv')
  
  d$cell_id <- as.factor(d$cell_id)
  d <- d[d$cell_id %in% ids,]
  
  ggplot(d, aes(x=cell_id, y=metric*1e6)) + 
    geom_bar(stat='identity', colour="black", fill=fills) +
    labs(x="",y="Arbor Complexity (mm-2)") +
    theme_classic() +
    scale_y_continuous(expand=c(0,0),limits=c(0,max(d$metric)*1e6+5)) +
    theme(text=element_text(size=30),#, family="Arial"),
          axis.text = element_text(size=16),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))
        
}

plot_diam_dist_distribution <- function(id, cell_info, id_mapping, percentile=0.1 ){
  
  mat <- load_sample_skeleton_mat(id)
  
  diams <- mat$rad
  
  dists <- read_node_distances(id, cell_info, id_mapping, mat)
  num_points <- length(dists)
  
  value_ranking <- order(order(dists))
  diams <- diams[ value_ranking < percentile * length(diams) ]
  dists <- dists[ value_ranking < percentile * length(dists) ]
  #diams <- diams[ dists < 500 ]
  #dists <- dists[ dists < 500 ]
  
  #cat(length(diams) / num_points)
  
  d <- data.frame(Distance=dists, Diameter=diams) 
  
  ggplot(d, aes(x=Distance, y=Diameter)) +
    geom_point(alpha=0.05) +
    geom_smooth() +
    labs(title=id)
}
