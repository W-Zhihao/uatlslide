########################################################
### Weighting function
########################################################
### This function is to simulate the relation between 
### the contribution (i.e., weight) of the model trained
### by the data from the target area and the number of
### data from the target area.
### 
### Input:
### - c_weight (numeric): 
###   control the speed of weight value increasing for
###   the model trained by the data from the target area
### - batch_size (numeric): 
###   the number of selected training set
### - epoch (numeric):
###   the decision boundary
### Output:
### - weight_al (numeric):
###   weight for the model trained by the data 
###   from the target area
########################################################
### (c) 2023 Zhihao Wang
########################################################


weight_AL_function <- function(c_weight, batch_size, epoch){
  
  stopifnot(is.numeric(c_weight))
  stopifnot(is.numeric(epoch))
  stopifnot(is.numeric(batch_size))
  
  weight_al <- 1 - 1 / (1 + c_weight * batch_size * epoch)
  
  return(weight_al)
  
}

