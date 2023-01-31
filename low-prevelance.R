########################################################
### Low-prevalence Margin Sampling (low-prevalence MS)
########################################################
### This active learning strategy mainly focus on 
### selecting the points, which spread around 
### the decision boundary bases on the
### posterior probabilities obtained by the active 
### learning learner.
### Input:
### - proba (data frame): 
###  posterior probabilities of unlabeled data set
### - batch_size (numeric): 
###   the number of selected training set
### - p (numeric, from 0 to 1):
###   the decision boundary
### Output:
### - select_loc (vector):
###   locations of "informative" points from unlabeled data set
########################################################
### (c) 2023 Zhihao Wang
########################################################


low_pre_MS <- function(proba, batch_size, p){
  
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(batch_size))
  
  select_loc <- order(abs(proba[,'TRUE'] - p))[1:batch_size]
  
  return(select_loc)
  
}