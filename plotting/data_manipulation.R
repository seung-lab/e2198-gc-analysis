library(ggplot2)
library(R.matlab)
library(reshape2)

#=====================================================================

#Ca imaging variables
timepoints_per_trial_per_direction = 31
num_trials = 5
num_directions = 8
short_direction = 6 #direction with shorter number of frames
num_repeated_frames = 3
stim_order <- c(180,0,225,45,270,90,315,135) * (2*pi)/360
stim_times <- c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30)


#=====================================================================
source("data_IO.R")
source("id_utils.R")

#Takes a for_ids fn (see below) and runs it over the ids within
# a type list
metric_for_type <- function( type_l, for_ids_fn, ... ){
  ids <- fetch_ids( type_l )
  types <- ids$ret_types
  ids <- ids$ids
  
  metric <- for_ids_fn( ids, ... )
  
  list("metric"=metric, "ids"=ids, "types"=types)
}


#rolled my own binning fn, bc this seems missing from ggplot
# also normalizes the output to a density height
# assumes that steps between xs are evenly spaced
bin_avg <- function( xs, ys, bin_width ){
    
  new_xs <- c()
  new_ys <- c()
  b <- 1
  e <- b + bin_width - 1
  num_values <- length(xs)
  
  if( num_values == 0 ){
    return(list("xs"=c(),"ys"=c()))
  }
  
  while( e < num_values ){
    new_xs <- c(new_xs, mean(xs[b:e]))
    new_ys <- c(new_ys, mean(ys[b:e]))
    
    b = b + bin_width
    e = e + bin_width
  }
  
  #last window
  new_xs <- c(new_xs, mean(xs[b:num_values]))
  new_ys <- c(new_ys, mean(ys[b:num_values]))
  
  if( length(new_xs) > 1){
    spacing <- abs(new_xs[2] - new_xs[1])
    new_ys <- new_ys / spacing
  }
  res <- list("xs"=new_xs, "ys"=new_ys)
  res
}

#=====================================================================
#For-ids fns (IMPORTANT)
# exchanges a list of cell_ids + supplementary data structures
# for the desired metric values
#
# If reading is quick (.csv), it's done every time the fn's called, but otherwise
# requires loading the supplementary data beforehand (.mat's)

#Results in cubic microns
soma_sizes_for_ids <- function( ids ){
  
  soma_d <- read.csv('completed_051016_measurements_t3.csv',sep=';',header=F)
  colnames(soma_d) <- c("ID","type","radius","coords","vol","type2","type3")
  #doing this in a rush...going to be ugly
  
  volumes <- NULL
  for( i in 1:length(ids) ){
    new_vol <- soma_d$vol[soma_d$ID == ids[i]]
    volumes <- c(volumes, new_vol)
  }
  
  volumes * (16.5*16.5*23)*(10*10*10)/(1e9) #num_voxels -> spatial volume (um^3)
}

#fetches the stratification profiles from a list of lists of lists of arrays
# (which is what happens when you read in this skel_strat.mat file)
# results are binned according to the binwidth
#Supplementary data file: skel_strat.mat
strats_for_ids <- function( ids, types, skel_strat, binwidth=default_binwidth ){
  
  strats <- data.frame()
  
  for( i in 1:length(ids) ){
    
    cell_id <- ids[i]
    cell_type <- types[i]
    #no clue why this many list indices are necessary
    cell_strat <- skel_strat[[1]][[cell_id]][[1]] 
    
    binned <- bin_avg(cell_strat[,1],cell_strat[,2],binwidth)
    # cat(sum(binned$ys * (binned$xs[1] - binned$xs[2]) * binwidth))
    
    strats <- rbind(strats, data.frame("xs"=binned$xs, "ys"=binned$ys, cell_id, cell_type))
  }
  
  colnames(strats) <- c("Depth","Stratification","cell_id", "cell_type")
  
  strats$Depth <- as.numeric(strats$Depth)
  
  #truncating IPL depth from [0,100], scaling this
  # range to [0,1], and scaling the strat heights
  # to match the change in scale for a density interpretation
  strats$Depth <- strats$Depth / 100 
  strats$Stratification <- strats$Stratification * 100 
  strats <- strats[ strats$Depth > -0.1 & strats$Depth < 1.1,]
  
  strats$Stratification <- as.numeric(strats$Stratification)
  strats$cell_id <- as.factor(strats$cell_id)
  
  strats
}

#Avg Ca trace over directions, trials (roi_sums_all)
overall_avg_ca_traces_for_ids <- function( ids, types, avgs, cell_ids ){
  
  traces <- data.frame()
  
  time <- 1:(dim(avgs)[1]/8) * 0.128
  time <- rep(time, 8)
  for( i in 1:length(ids) ){
    new_trace <- avgs[,ids[i]]
    new_type <- types[i]
    
    new_rows <- data.frame( time, new_trace, ids[i], types[i], cell_ids[i] )
    traces <- rbind( traces, new_rows )
  }
  colnames(traces) <- c("Time","Trace","Ca_id","cell_type","cell_id")
  
  traces
}

#---------------------------------------------------------------------
#Ca tuning curve fns

extract_tuning <- function( tuning, ca_ids ){
  tuning[,,ca_ids]
}

extract_traces <- function( roi_sums, ca_ids ){
  roi_sums[,ca_ids]
}

remove_extra_trace_frames <- function( trace_dframe ){
  
  #within a trial, where are the repeated frames?
  repeat_offset <- timepoints_per_trial_per_direction * 6
  
  frames_to_remove <- NULL
  for( i in 1:num_trials ){
    trial_offset <- (i-1) * timepoints_per_trial_per_direction * num_directions
    full_offset <- trial_offset + repeat_offset + 1:3
    
    frames_to_remove <- c(frames_to_remove, full_offset)
  }
  
  if( length(dim(trace_dframe)) == 1 ){
    trace_dframe[-frames_to_remove]
  } else {
    trace_dframe[-frames_to_remove,]
  }
}

assign_trials <- function(){
  res <- NULL
  timepoints_per_trial <- timepoints_per_trial_per_direction * num_directions
  for( i in 1:num_trials ){
    res <- c(res, rep(i,timepoints_per_trial))
  }
  res
}

assign_timepoints <- function(){
  timepoints_per_trial <- timepoints_per_trial_per_direction * num_directions
  rep(0.128 * 1:timepoints_per_trial, num_trials)
}

assign_stimuli <- function(){
  rep(rep(rep(c(1,2,3,4,5,6,7,8),each=timepoints_per_trial_per_direction), 8),5)
}

#---------------------------------------------------------------------
#Cell_info based metrics 

#Generalized fn to fetch a metric from the cell_info data
# USE WITH CARE - sapply can be finicky
fetch_ci_metric_for_ids <- function( ids, cell_info, id_mapping, metric_index ){
  
  column_indices <- id_mapping[as.character(ids)]
  sapply(cell_info[metric_index, column_indices, ], function(x){x})
}

#Max Diameter over warped 2D convex hull (in nm)
diameters_for_ids <- function( ids, cell_info, id_mapping ){
  fetch_ci_metric_for_ids( ids, cell_info, id_mapping, 14 ) * 66 #scaling to nm
}

#Area over warped 2D convex hull (in nm)
hull_areas_for_ids <- function( ids, cell_info, id_mapping ){
  fetch_ci_metric_for_ids( ids, cell_info, id_mapping, 15 ) * 66*66 #scaling to nm
}

#Original Volume-Based Centroid Asymmetry (in log scale)
centroid_asym_for_ids <- function( ids, cell_info, id_mapping ){
  log(fetch_ci_metric_for_ids( ids, cell_info, id_mapping, 17 )) #log scaling
}

#Original Volume-Based Dorsoventral (2an direction) Asymmetry (in log scale)
dv_asym_for_ids <- function( ids, cell_info, id_mapping ){
  log(fetch_ci_metric_for_ids( ids, cell_info, id_mapping, 18 ))
}

#Vector from soma center -> 2D warped hull centroid
# NOTE: this fn has to be a bit different bc we manipulate a 2d cell_info
# metric a bit
asym_vectors_for_ids <- function( ids, cell_info, id_mapping ){
  
  column_indices <- id_mapping[as.character(ids)]
  
  vectors <- NULL
  for( i in 1:length(ids) ){
    #using asymm_index
    new_vector <- cell_info[11,column_indices[i],][[1]]
    #originally scaled by radius of convex hull
    new_vector <- new_vector * sqrt( cell_info[,column_indices[i],]$area.hull / pi)[[1]]
    #scaling z
    new_vector <- new_vector * c(1, 25/16)
    
    vectors <- rbind(vectors, new_vector)
  }
  
  vectors
}


#Angle of the vector from soma center -> 2D warped hull centroid
asymmetry_angle_for_ids <- function( ids, cell_info, id_mapping ){
  column_indices <- id_mapping[as.character(ids)]
  
  #metric_index = 11
  vector_list <- cell_info[11, column_indices, ]
  
  angles <- NULL
  mean_vec <- c(0,0)
  for( i in 1:length(vector_list)){
    vec <- vector_list[[i]]
    mean_vec <- mean_vec + vec
    
    angles <- c(angles, atan2(vec[2],vec[1]))
  }
  
  #Can be nice to print
  cat("mean angle\n")
  cat((atan2(mean_vec[2],mean_vec[1])) * (180/pi) %% 360)
  
  angles
}

#Magnitude of the vector from soma center -> 2D warped hull centroid
# OUT OF DATE - doesn't take hull radius scaling into acct atm
asymmetry_magnitudes_for_ids <- function( ids, cell_info, id_mapping ){
  column_indices <- id_mapping[as.character(ids)]
  
  #metric_index = still 11
  vector_list <- cell_info[11,column_indices, ]
  
  magnitudes <- NULL
  for( i in 1:length(vector_list) ){
    vec <- vector_list[[i]]
    
    magnitudes <- c(magnitudes, norm(vec,"2"))
  }
  
  magnitudes
}

#Warped soma coordinate in YZ
soma_proj_for_ids <- function( ids, cell_info, id_mapping, warped=F ){
  
  column_indices <- id_mapping[as.character(ids)]
  
  soma_coords <- NULL
  for( i in 1:length(ids) ){
    
    if(warped){
      #using soma_coords_warped_mip2_zscaled
      new_centroid <- cell_info[13,column_indices[i],][[1]]
      
    } else {
      #using soma_coord
      new_centroid <- cell_info[9,column_indices[i],][[1]]
      
    }
    
    soma_coords <- rbind(soma_coords, new_centroid)
  }
  
  if(warped){
    #NOTE warped coordinates are in YZX
    soma_coords[,1:2]
  } else {
    soma_coords[,2:3]
  }
}


#soma w/o projection to yz plane
soma_for_ids <- function( ids, cell_info, id_mapping, warped=F ){
  
  column_indices <- id_mapping[as.character(ids)]
  
  soma_coords <- NULL
  for( i in 1:length(ids) ){
    
    if(warped){
      #using soma_coords_warped_mip2_zscaled
      new_centroid <- cell_info[13,column_indices[i],][[1]]
      
    } else {
      #using soma_coord
      new_centroid <- cell_info[9,column_indices[i],][[1]]
      
    }
    
    soma_coords <- rbind(soma_coords, new_centroid)
  }
  
  if(warped){
    #NOTE warped coordinates are in YZX
    soma_coords[,c(3,1,2)]
  } else {
    soma_coords[,c(1,2,3)]
  }
}


#---------------------------------------------------------------------
#Skeleton file based metrics
# these are mostly self-explanatory, so documentation is minimal

# = num_edges
read_skel_path_length <- function( id, mat=NULL ){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
  }
  
  sources <- mat$n[ mat$e[,1], ]
  dests   <- mat$n[ mat$e[,2], ]
  
  coord_dists <- dests - sources
  
  phys_dists <- sweep(coord_dists, MARGIN=2, 4*c(23,16.5,16.5), '*')
  edge_lengths <- sqrt(rowSums(phys_dists^2))
  
  sum(edge_lengths)
}

read_skel_num_nodes <- function( id, mat=NULL ){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
  }
  
  result <- dim(mat$n)[1]
  
  if( is.null(result) ){
    result <- NA
  }
  
  result
}

read_skel_branch_nodes <- function( id, mat=NULL ){
  
   if( is.null(mat) ){
     mat <- load_skeleton_mat(id)
   }
  
  result <- dim(mat$bn)[1]
  
  if( is.null(result) ){
    result <- NA
  }
  
  result
}

read_median_branch_width <- function( id, mat=NULL ){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
  }
  
  result <- median(mat$rad, na.rm=TRUE)
  
  if( is.null(result) ){
    result <- NA
  }
  
  result
}

read_mean_branch_width <- function( id, mat=NULL ){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
  }
  
  result <- mean(mat$rad, na.rm=TRUE)
  
  if( is.null(result) ){
    result <- NA
  }
  
  result
}


read_max_branch_width <- function( id, mat=NULL ){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
  }
  
  result <- max(mat$rad, na.rm=TRUE)
  
  if( is.null(result) ){
    result <- NA
  }
  
  result
}

max_branch_width_for_ids <- function( ids ){
  mbws <- NULL
  for( i in 1:length(ids) ){
    writeLines(sprintf("%d",ids[i]))
    mbws <- c(mbws, read_max_branch_width(ids[i]))
  }
  mbws
}

median_branch_width_for_ids <- function( ids ){
  mbws <- NULL
  for( i in 1:length(ids) ){
    writeLines(sprintf("%d",ids[i]))
    mbws <- c(mbws, read_median_branch_width(ids[i]))
  }
  mbws
}

mean_branch_width_for_ids <- function( ids ){
  mean_bws <- NULL
  for( i in 1:length(ids) ){
    writeLines(as.character(ids[i]))
    mean_bws <- c(mean_bws, read_mean_branch_width(ids[i]))
  }
  mean_bws
}

path_lengths_for_ids <- function( ids ){
  
  path_lengths <- c()
  for( i in 1:length(ids) ){
    writeLines(sprintf("Cell %d",ids[i]))
    path_lengths <- c(path_lengths, read_skel_path_length(ids[i]))
  }
  
  path_lengths
}

bp_per_path_length_for_ids <- function( ids ){
  
  quotients <- NULL
  for( i in 1:length(ids) ){
    writeLines(sprintf("%d",ids[i]))
    
    mat <- load_skeleton_mat(ids[i])
    bp <- read_skel_branch_nodes( ids[i], mat )
    pl <- read_skel_path_length( ids[i], mat )
    
    quotients <- c(quotients, bp/pl)
  }
  
  quotients
}

bp_for_ids <- function( ids ){
  
  bp_counts <- c()
  for( i in 1:length(ids) ){
    writeLines(sprintf("%d",ids[i]))
    bp_counts <- c(bp_counts, read_skel_branch_nodes( ids[i] ) )
  }
  
  bp_counts
}

bp_per_hull_area_for_ids <- function( ids, cell_info, id_mapping ){
  
  bp_counts <- bp_for_ids(ids)
  quotients <- bp_counts / hull_areas_for_ids( ids, cell_info, id_mapping )
  
  quotients
}

num_nodes_for_ids <- function( ids ){
  node_counts <- NULL
  for( i in 1:length(ids) ){
    writeLines(sprintf("%d",ids[i]))
    node_counts <- c( node_counts, read_skel_num_nodes(ids[i]))
  }
  node_counts
}

#hybrid skeleton ci metrics
nodes_per_area_for_ids <- function( ids, cell_info, id_mapping ){
  
  node_counts <- num_nodes_for_ids( ids )
  quotients <- node_counts / hull_areas_for_ids( ids, cell_info, id_mapping )
  
  quotients
}

read_node_distances <- function( id, cell_info, id_mapping, mat=NULL){
  
  if( is.null(mat) ){
    mat <- load_skeleton_mat(id)
    #mat <- load_sample_skeleton_mat(id)
  }
  
  #node coords in ZYX
  node_coords <- mat$n
  node_coords <- node_coords[,c(3,2,1)] #flipping back to XYZ
  node_coords[,3] <- node_coords[,3] * 23 / 16.5 #scaling to match soma
  
  #soma's returned as XYZ
  soma_coord <- soma_for_ids( id, cell_info, id_mapping, T )
  
  #euc distance
  sqrt(colSums(apply(node_coords,1,'-',soma_coord)^2))
}

primary_dendrite_width_for_ids <- function(ids, cell_info, id_mapping, percentile=0.1){
  
  widths <- NULL
  for( id in ids ){
    
    cat(id)
    cat("\n")
    
    mat <- load_skeleton_mat(id)
    
    diams <- mat$rad
    
    dists <- read_node_distances(id, cell_info, id_mapping, mat)
    num_points <- length(dists)
    
    value_ranking <- order(order(dists))
    diams <- diams[ value_ranking < percentile * length(diams) ]
    
    widths <- c(widths, mean(diams))
  }
  
  widths
}

#---------------------------------------------------------------------
#gcsac fns

circshift <- function(x,k){
  k <- k %% length(x)
  c( tail(x,k), head(x,-k) )
}

bin_percentages <- function( ps ){
  binwidth <- 45
  
  #shifting to match axes?
  ps <- circshift(ps, 23)
  
  b <- 1
  e <- b + binwidth - 1
  len <- length(ps)
  new_ps <- NULL
  while( e < len ){
    new_ps <- c( new_ps, sum(ps[b:e], na.rm=T) )
    b <- b + binwidth
    e <- e + binwidth
  }
  
  #last_window
  new_ps <- c(new_ps, sum(ps[b:len], na.rm=T))
  
  new_ps
}

form_vectors <- function( radii, angles ){
  
  vectors <- NULL
  for( i in 1:length(radii) ){
    new_vec <- c( cos(angles[i])*radii[i], sin(angles[i])*radii[i] )
    vectors <- rbind( vectors, new_vec)
  }
  
  vectors
}

fetch_dss_vector <- function( id, gcsac, num_ids, denom_ids, mode="both" ){
  
  #mode
  #"both" = on+off
  #"on"   = on, 
  #"of"   = off 
  #(too lazy to do this better)
  
  num_index <- num_ids[as.character(id)]
  denom_index <- denom_ids[as.character(id)]
  
  numerator <- gcsac$gc.num.vals[[num_index]][[1]]
  denominator <- gcsac$gc.denom.vals[[denom_index]][[1]]
  
  numerator <- list(
    "ON"=numerator[,1],
    "OFF"=numerator[,2]
  )
  denominator <- list(
    "ON"=denominator[,1],
    "OFF"=denominator[,2]
  )
  
  num_binned <- lapply( numerator, bin_percentages )
  denom_binned <- lapply( denominator, bin_percentages )
  
  #for now, summing across ON and OFF
  if( mode == "both" ){
    num_binned <- num_binned$ON + num_binned$OFF
    denom_binned <- denom_binned$ON + denom_binned$OFF
  } else if( mode == "on" ){
    num_binned <- num_binned$ON
    denom_binned <- denom_binned$ON
  } else { #assumed off
    num_binned <- num_binned$OFF
    denom_binned <- denom_binned$OFF
  }
  
  percentages <- num_binned / denom_binned
  percentages[is.nan(percentages)] <- 0
  
  #angles <- seq(0,2*pi,by=2*pi/length(percentages[[1]]))
  angles <- seq(0,2*pi,by=2*pi/length(percentages))
  
  #vectors <- lapply(binned, form_vectors, angles)
  vectors <- form_vectors(percentages, angles)
  
  #lapply(vectors, colSums) / lapply(binned, sum)
  colSums(vectors) / sum(percentages)
}

dss_vectors_for_ids <- function(ids, gcsac, num_ids, denom_ids, mode="both"){
  
  vectors <- NULL
  for( id in ids ){
    vectors <- rbind( vectors, fetch_dss_vector(id, gcsac, num_ids, denom_ids, mode))
  }
  
  vectors
}

dssi_for_ids <- function( ids, gcsac, num_ids, denom_ids, mode ){
  apply( dss_vectors_for_ids(ids, gcsac, num_ids, denom_ids, mode), 1, norm, "2" )
}

dss_angles_for_ids <- function( ids, gcsac, num_ids, denom_ids, mode ){
  
  vecs <- dss_vectors_for_ids(ids, gcsac, num_ids, denom_ids, mode)
  
  angles <- NULL
  for( i in 1:dim(vecs)[1] ){
    new_vec <- vecs[i,]
    
    angles <- c(angles, atan2(new_vec[2],new_vec[1]))
  }
  
  angles
}

#Normalized SAC Contact (contact voxels over nearby area)
# NOTE this currently uses a ratio of averages approach
# NOTE#2 this should ONLY be used for overall ratios, not directionality
fetch_normalized_sac_contact <- function( id, gcsac, num_ids, denom_ids, combine=F ){
  
  num_index <- num_ids[as.character(id)]
  denom_index <- denom_ids[as.character(id)]
  
  numerator <- colSums(gcsac$gc.num.vals[[num_index]][[1]])
  denominator <- colSums(gcsac$gc.denom.vals[[denom_index]][[1]])
  if(combine){
    numerator <- sum(numerator)
    denominator <- sum(denominator)
  }
  
  numerator / denominator
}

sac_contact_for_ids <- function( ids, gcsac, num_ids, denom_ids, onoff_index ){
  
  contacts <- NULL
  for( i in 1:length(ids) ){
    sac_contact <- fetch_normalized_sac_contact( ids[i], gcsac, num_ids, denom_ids )
    contacts <- c(contacts, sac_contact[onoff_index])
  }
  
  contacts
}

on_contact_for_ids <- function( ids, gcsac, num_ids, denom_ids ){
  sac_contact_for_ids(ids, gcsac, num_ids, denom_ids, 1)
}

off_contact_for_ids <- function( ids, gcsac, num_ids, denom_ids ){
  sac_contact_for_ids(ids, gcsac, num_ids, denom_ids, 2)
}

both_contact_for_ids <- function( ids, gcsac, num_ids, denom_ids ){
  contacts <- NULL
  for( i in 1:length(ids) ){
    sac_contact <- fetch_normalized_sac_contact( ids[i], gcsac, num_ids, denom_ids, combine=T )
    contacts <- c(contacts, sac_contact)
  }
  contacts
}

#---------------------------------------------------------------------
#DSI OSI fns

fetch_normalized_selectivity_index <- function( id, dsos, id_lookup, ds_not_os, combine=F ){
  
  if( ds_not_os ){
    radii <- dsos$ds.r
    theta <- dsos$ds.theta
  } else {
    radii <- dsos$os.r
    theta <- dsos$os.theta
  }
  
  index <- id_lookup[as.character(id)]
  radius <- radii[index,]
  angle  <- theta[index,]
  denominator <- dsos$r.mean[index,]
  
  if(combine && !is.na(radius[1])){ #sum over ON and OFF values
    on_num_v  <- c( radius[1] * cos(angle[1]), radius[1] * sin(angle[1]) )
    off_num_v <- c( radius[2] * cos(angle[2]), radius[2] * sin(angle[2]) )
    
    numerator <- norm( on_num_v + off_num_v, "2" )
    denominator <- sum(denominator)
  }
  
  if(!is.na(radius[1])){
    numerator / denominator
  } else {
    NA
  }
}

selectivity_index_for_ids <- function( ids, dsos, id_lookup, ds_not_os, onoff_index ){
  
  indices <- NULL
  for(i in 1:length(ids)){
    
    on_and_off <- fetch_normalized_selectivity_index( ids[i], dsos, id_lookup, ds_not_os )
    indices <- c(indices, on_and_off[onoff_index])
  }
  indices
}

overall_selectivity_index_for_ids <- function( ids, dsos, id_lookup, ds_not_os ){
  indices <- NULL
  for(i in 1:length(ids)){
    indices <- c(indices,
                 fetch_normalized_selectivity_index( ids[i], dsos, id_lookup, ds_not_os, T))
  }
  indices
}

on_dsi_for_ids <- function( ids, dsos, id_lookup ){
  selectivity_index_for_ids( ids, dsos, id_lookup, T, 1)
}

off_dsi_for_ids <- function( ids, dsos, id_lookup ){
  selectivity_index_for_ids( ids, dsos, id_lookup, T, 2)
}

overall_dsi_for_ids <- function( ids, dsos, id_lookup ){
  overall_selectivity_index_for_ids( ids, dsos, id_lookup, T)
}

on_osi_for_ids <- function( ids, dsos, id_lookup ){
  selectivity_index_for_ids( ids, dsos, id_lookup, F, 1 )
}

off_osi_for_ids <- function( ids, dsos, id_lookup ){
  selectivity_index_for_ids( ids, dsos, id_lookup, F, 2)
}

overall_osi_for_ids <- function( ids, dsos, id_lookup ){
  overall_selectivity_index_for_ids( ids, dsos, id_lookup, F)
}

dsi_angles_for_ids <- function( ids, dsos, id_lookup, ds_not_os=T ){
  
  if( ds_not_os ){
    radii <- dsos$ds.r
    theta <- dsos$ds.theta
  } else {
    radii <- dsos$os.r
    theta <- dsos$os.theta
  }
  
  index <- id_lookup[as.character(ids)]
  radius <- radii[index,]
  angle  <- theta[index,]
  
  if(length(index) > 1){
    on_v  <- cbind( radius[,1] * cos(angle[,1]), radius[,1] * sin(angle[,1]) )
    off_v <- cbind( radius[,2] * cos(angle[,2]), radius[,2] * sin(angle[,2]) )
  } else {
    on_v  <- cbind( radius[1] * cos(angle[1]), radius[1] * sin(angle[1]) )
    off_v <- cbind( radius[2] * cos(angle[2]), radius[2] * sin(angle[2]) )
  }
    
  overall_v <- on_v + off_v
  
  -atan2(overall_v[,2],overall_v[,1])
}

#---------------------------------------------------------------------
#Skeleton-based asymmetry

#Fetches signed/unsigned versions of asymmetry index
# unsigned = (pos_projections) / (#skeleton_nodes)
# signed   = (pos_projections-neg_projections) / (#skeleton_nodes)
asym_index_for_ids <- function(ids, asym, asym_ids, dv = T, signed=T){
  
  indices <- asym_ids[as.character(ids)]
  
  if(dv){
    selected_index <- asym$asym_2an_proj
  } else {
    selected_index <- asym$asym_type
  }
  
  if(signed){
    2* selected_index[indices] - 1
  } else {
    selected_index[indices]
  }
}

#Fetches the arbor vector in the mip2 warped (and zscaled) space
arbor_vectors_for_ids <- function(ids, asym, asym_ids){
  
  indices <- asym_ids[as.character(ids)]
  
  vectors <- NULL
  for(i in 1:length(ids)){
    new_vector_str <- as.character(asym$arbor_vector[indices[i]])
    vectors <- rbind( vectors, map_string_to_vector(new_vector_str) )
  }
  
  vectors
}

arbor_angles_for_ids <- function(ids, asym, asym_ids){
  
  arbor_vectors <- arbor_vectors_for_ids(ids, asym, asym_ids)
  angles <- NULL
  #mean_vector <- c(0,0)
  
  for( i in 1:length(ids) ){
    arbor_vector <- arbor_vectors[i,]
    #mean_vector <- mean_vector + c(arbor_vector[1], arbor_vector[2])
    # -1 flips to DVRC orientation
    angles <- c(angles,-atan2(arbor_vector[2], arbor_vector[1])) 
  }
  
  
  #Can be nice to print
  #cat("mean angle\n")
  #cat((-atan2(mean_vector[2],mean_vector [1])) * (180/pi) %% 360)
  
  angles
}

arbor_magnitudes_for_ids <- function(ids, asym, asym_ids){
  
  arbor_vectors <- arbor_vectors_for_ids(ids, asym, asym_ids)
  magnitudes <- NULL
  
  for( i in 1:length(ids) ){
    arbor_vector <- arbor_vectors[i,]
    magnitudes <- c(magnitudes, norm(arbor_vector,"2"))
  }
  
  magnitudes
}

#Maps a string in the format "[a,b]" to c(a,b) as a numeric
map_string_to_vector <- function( vector_str ){
  
  values <- substr( vector_str, 2,nchar(vector_str)-1 )
  values <- strsplit( values, ",")[[1]]
  
  as.numeric(values)
}
