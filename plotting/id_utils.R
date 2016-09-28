library(ggplot2)
library(R.matlab)
library(reshape2)

#Function to convert the gc variable (imported into R) of Alex's 
# final_classification
make_id_csv <- function(gc, output_name){
  
  num_types <- dim(gc)[2]
  
  type_df <- data.frame()
  for( i in 1:num_types ){
    type_info <- gc[,i,]
    
    new_data <- data.frame("type"=type_info$name[[1]], 
                           "cell_id"=as.vector(type_info$cells))
    
    type_df <- rbind(type_df, new_data)
  }
  cat("here")

  write.csv(type_df, output_name, row.names=F) 
}

#Function to fetch ids for a given iterable of types
# utilizing the final_classification.csv
fetch_ids <- function(types){
  
  classification <- read.csv("final_final_classification.csv")
  
  ids <- c()
  ret_types <- c()
  for(i in 1:length(types)){
    new_ids = classification$cell_id[classification$type == types[i]]
    
    ret_types <- c(ret_types, rep(types[i], length(new_ids)))
    ids = c(ids, new_ids)
  }
  
  rm(classification)
  
  data.frame(ids, ret_types)
}

fetch_types <- function(ids){
  
  classification <- read.csv("final_classification.csv")
  
  ret_types <- c()
  for( i in 1:length(ids)){
    new_type <- classification$type[classification$cell_id == as.character(ids[i])]
    ret_types <- c(ret_types, new_type)
  }
  
  ret_types <- levels(classification$type)[ret_types]
  rm(classification)
  data.frame(ids, ret_types)
}

#translates from types to ids without extra fluff
type_ids <- function(type_l){
  fetch_ids(type_l)$ids
}

id_types <- function(id_l){
  fetch_types(id_l)$ret_types
}

soma_id_to_omni_id <- function(ids){
  
  id_map <- read.csv("omni_id_to_soma_id.csv", sep=";")
  
  ret_ids <- c()
  for(i in 1:length(ids)){
    ret_ids <- c(ret_ids, id_map$omni_id[id_map$soma_id == ids[i]])
  }
  
  ret_ids
}

omni_id_to_soma_id <- function(ids){
  
  id_map <- read.csv("omni_id_to_soma_id.csv", sep=";")
  
  ret_ids <- c()
  for(id in ids){
    ret_ids <- c(ret_ids, id_map$soma_id[id_map$omni_id == id])
  }
  ret_ids
}

soma_id_types <- function(ids){
  fetch_types(soma_id_to_omni_id(ids))$ret_types
}
#Translates the omni ids into ca recording ids using a mapping array
# (this can be loaded using avg_ca.mat / load_avg_ca() )
fetch_ca_ids <- function( ids, ca_mapping_array ){
  
  ca_indices <- c()
  for(i in 1:length(ids)){
    ca_indices <- c( ca_indices, which( ca_mapping_array[,2] == ids[i] )[1] )
  }
  
  ca_mapping_array[ca_indices,1]
}
