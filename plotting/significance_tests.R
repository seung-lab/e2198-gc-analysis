
library(circular)

source("data_manipulation.R")

#		log10.p.value <- min((-z)/log(10) + log10(temp), 0)

type_asymmetry_significance_test <- function( type_l, asym, asym_ids ){
  
  pvals <- NULL
  for( i in 1:length(type_l) ){
    ids <- type_ids(type_l[i])
    angles <- arbor_angles_for_ids( ids, asym, asym_ids )
    
    res <- rayleigh.test(circular(angles))
    pvals <- c(pvals, res$p.value)
  }
  
  names(pvals) <- type_l
  pvals
}

type_dsi_significance_test <- function( type_l, dsos, dsos_ids ){
  
  pvals <- NULL
  for( type in type_l ){
    
    ids <- type_ids(type)
    angles <- dsi_angles_for_ids( ids, dsos, dsos_ids )
    
    res <- rayleigh.test(circular(angles))
    pvals <- c(pvals, res$p.value)
  }
  
  names(pvals) <- type_l
  pvals
}

#felt pretty aggressive with these
cutoffs <- list("2an"=c(26129, 26190, 26041, 26082, 26147, 17130, 17177, 17062, 20024, 26172),
                "2aw"=c(26038, 26189, 17075, 26163, 26150, 26110, 26131, 26095, 17200, 17061, 26055, 17060, 26018),
                "63" =c(26148, 26057, 26068, 26089, 26191, 26141, 26027, 26028, 26125),
                "51" =c(26177, 26154, 26122, 26136, 26039, 17098, 26085, 26113, 17035),
                "25" =c(26175, 26145, 17176, 26117, 26066, 26040, 26134, 26099, 26031),
                "6sw"=c(17083, 26020),
                "82n"=c(26080,26072,26076),
                "6sn"=c(26035,17082, 26171, 20198),
                "3o" =c(17076, 26155))


type_asymmetry_significance_test_wo_cutoffs <- function( type_l, asym, asym_ids ){
  
  pvals <- NULL
  for( i in 1:length(type_l) ){
    ids <- type_ids(type_l[i])
    
    #Defining a filtering fn
    not_in_cutoffs_list <- function(x) {
      if( x %in% cutoffs[type_l[i]][[1]] ) return(FALSE) else return(TRUE)
    }
    filter <- Vectorize(not_in_cutoffs_list)
    
    ids <- ids[filter(ids)]
    angles <- arbor_angles_for_ids( ids, asym, asym_ids )
    
    res <- rayleigh.test(circular(angles))
    pvals <- c(pvals, res$p.value)
  }
  
  names(pvals) <- type_l
  pvals
}

type_dssi_significance_test <- function( type_l, gcsac, num_ids, denom_ids, mode="both" ){
  
  pvals <- NULL
  for( type in type_l ){
    
    ids <- type_ids(type)
    angles <- dss_angles_for_ids( ids, gcsac, num_ids, denom_ids, mode )
    
    res <- rayleigh.test(circular(angles))
    pvals <- c(pvals, res$statistic)
  }
  
  names(pvals) <- type_l
  pvals
}
