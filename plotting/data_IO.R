library(ggplot2)
library(R.matlab)
library(reshape2)

#Location of completed skeleton mat files
skel_mat_fmt = "~/seungmount/research/Alex/e2198/Cells_done/skel_%d.mat"
sample_skel_mat_fmt = "~/Dropbox/seunglab/plotting_retina_figures/sample_cells/skel_%d.mat"

#=====================================================================
#Data IO

#Creates a lookup array for cell_ids to make my life easier
make_id_lookup <- function( id_list ){
  
  id_list <- as.vector(id_list)
  id_lookup <- 1:length(id_list)
  names(id_lookup) <- id_list
  
  id_lookup
}

#Contains skeleton-based asymmetry scores
# and arbor vectors
load_asym <- function(){
  asym <- read.csv("asymmetry.csv",header=F,sep=";")
  #only for testing the effect of cutoffs on asym measures
  #asym <- read.csv("asymmetry_no_numerous_cutoffs.csv",header=F,sep=";")
  asym <- asym[-ncol(asym)]
  
  colnames(asym) <- c("cell_id",
                      "arbor_vector","skel_centroid","soma_coord",
                      "asym_cent", "asym_2an_proj", "asym_type")
  
  ids <- make_id_lookup(asym$cell_id)
  
  list("asym"=asym, "ids"=ids)
}

#Contains GC-SAC contact numerators and denominators
# binned to each degree for each cell (indexed by id_lookups)
load_gcsac <- function() {
  gcsac <- readMat("gcsac_20160615.mat")
  
  num_id_lookup <- 1:length( gcsac$gc.num.keys )
  names(num_id_lookup) <- gcsac$gc.num.keys
  
  denom_id_lookup <- 1:length( gcsac$gc.denom.keys )
  names(denom_id_lookup) <- gcsac$gc.denom.keys
  
  list("gcsac"=gcsac, "num_id"=num_id_lookup, "denom_id"=denom_id_lookup)
}

#Loads the most recent cell_info mat file containing
# lots of stuff
load_cell_data <- function() {
  cell_info <- readMat("cell_info_clustering.20160623warpedSomaCorrection.mat")
  
  cell_info <- cell_info[[1]]
  
  cell_ids <- as.vector(cell_info[1,,])
  id_lookup <- 1:length(cell_ids)
  names(id_lookup) <- cell_ids
  
  list("cell_info"=cell_info, "id_lookup"=id_lookup)
}

#Contains the mean response vector (magnitude and direction) 
# averaged over either direction or orientation
load_dsos <- function(){
  dsos <- readMat("ca_dsos_forR3.mat")
  
  id_lookup <- make_id_lookup( dsos$omni.id )
  
  list("dsos"=dsos, "id_lookup"=id_lookup)
}

#Contains skeleton-based stratification profiles
# for each cell
load_strat <- function(){
  readMat("skel_strat.mat")
}

load_avg_ca <- function(){
  
  mat1 <- readMat("avg_ca.mat")
  detrended_mat <- readMat("roi_sums_means_flatten_detrended.mat")
  mat1$roi.sums.means.flatten <- detrended_mat$roi.sums.means.flatten
  
  mat1
}

load_roi_sums <- function(){
  #mat <- readMat("roi_data.mat")
  mat <- readMat("coeffs16.20160822.mat")
  mat$roi.sums.all
}

load_tuning <- function(){
  readMat("tuning.mat")[[1]]
}

load_skeleton_mat <- function( id ){
  fname <- sprintf(skel_mat_fmt, id)
  readMat(fname)
}

load_sample_skeleton_mat <- function( id ){
  fname <- sprintf(sample_skel_mat_fmt, id)
  readMat(fname)
}
