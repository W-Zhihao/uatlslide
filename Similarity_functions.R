# Case-based reasoning method

# This method is based on Wang et al.,(2022) in Geoscientific Model Development (https://doi.org/10.5194/gmd-15-8765-2022)

# It includes Similarity Functions and Landslide Susceptibility Modeling.

# Author: Zhihao Wang - zhihao.wang@uni-jena.de

#----------------------Similarity Functions---------------------------------------------------------------------

# Similarity functions for Geological characteristics, Data characteristics and Characteristics of study area

# @Description of parameters
# # source_area ---- a `list`, values of attributes for study area
# 
# # target_area ---- a `list`, values of attributes for study area


### Spatial Resolution

RESOLUTION_SIMILARITY_FUN <- function(source_area, target_area){

  if (source_area['spatial resolution'] >= target_area['spatial resolution']) {

    sim_resolution <- 2 ^ -(2 * abs(log10(target_area['spatial resolution'])-log10(source_area['spatial resolution']))) ^ 0.5

  }else{

    sim_resolution <- 1

  }

  return(sim_resolution)

}

### Total relief

TOTALRELIEF_SIMILARITY_FUN <- function(source_area, target_area){

  sim_totalrelief <- 1-(abs(target_area['total relief'] - source_area['total relief']) 
                        / max(8848-target_area['total relief'], target_area['total relief']))

  return(sim_totalrelief)

}

### Mean Slope

MEANSLOPE_SIMILARITY_FUN <- function(source_area, target_area, limit_slope){

  # Args:
  #   limit_slope: should cover mean slope of each study area, in our paper is 40 degrees.

  sim_meanslope <- 1 - (abs(target_area['mean slope'] - source_area['mean slope']) 
                        / max(limit_slope-target_area['mean slope'], target_area['mean slope']))

  return(sim_meanslope)

}

### STANDARD DEVIATION OF SLOPE

STDSLOPE_SIMILARITY_FUN <- function(source_area, target_area){

  sim_stdslope <-  2 ^ -(2*abs(log10(target_area['std slope'])-log10(source_area['std slope']))) ^ 0.5

  return(sim_stdslope)

}

### GEOLOGICAL UNITS

GEOLOGICALUNITS_SIMILARITY_FUN <- function(source_area, target_area){
  
  

  Igneous_sim = identical(source_area['igneous'], target_area['igneous'])

  Sedimentary_sim = identical(source_area['sedimentary'], target_area['sedimentary'])

  Metamorphic_sim = identical(source_area['metamorphic'], target_area['metamorphic'])

  sim_geounits = mean(c(Igneous_sim, Sedimentary_sim, Metamorphic_sim))

  return(sim_geounits)

}











