########################################################
### Adaptive unsupervised active-transfer learning 
### (adaptive UAT)
########################################################
### This function is to re-estimate the c_weight in  
### weight_function.R in each epoch based on the data
### acquired from the target area during the previous 
### epoch.
### 
### Input:
### - c_weight_range (numeric): 
###   the range of c_weight values
### - f_target (formula): 
###   the formula of the target area
### - train_set (data frame):
###   the data from the target area
### - clf.case (list): 
###   the model trained by the "suitable" source area
### - batch_size (numeric): 
###   the number of selected training set
### - epoch (numeric):
###   the decision boundary
### Output:
### - optimal.cwe (numeric):
###   the optimal c_weight for the model trained by the 
###   data from the target area in the current epoch
########################################################
### (c) 2023 Zhihao Wang
########################################################
library(mgcv)
library(sperrorest)
library(ROCR)
library(dplyr)

adaptive_UAT <- function(c_weight_range, f_target, train_set, clf.case, batch_size, epoch)
{
  
# Cite: Machine Learning for Geospatial Modelling,
# GEO 408B FSU Jena,Alexander Brenning and Jason Goetz
## A function for plotting ROC curves and
## calculating the area under the ROC curve (AUROC)
## using the ROCR package:
  
  auroc <- function(pred, obs, plot=FALSE) {
    stopifnot(is.logical(obs) | is.factor(obs))
    stopifnot(is.numeric(pred))
    stopifnot(length(pred) == length(obs))
    if (is.factor(obs)) stopifnot(nlevels(obs)==2)
    require(ROCR)
    predobj <- prediction( pred, obs )
    if (plot) plot(performance(predobj, "tpr", "fpr"))
    auroc <- performance( predobj, measure="auc", fpr.stop=0.1)@y.values[[1]]
    return(auroc)
  }
  
## change logit to probability
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }

## calculate the predictive accuracy for all c_weight values
  tunning <- function(c_weight, clf.target=NULL, clf.case=NULL, test=NULL, batch_size=NULL,
                     epoch=NULL){
    
    w.target <- 1 - 1 / (1 + c_weight * batch_size* epoch)
    
    if(w.target == 1){
      
      pred_o <- as.data.frame(predict(clf.target, newdata = test, type = "response"))
      
    }else{
      
      # standardize for combing the model trained by the "suitable" source area and
      # the model trained by the selected data from the target area
      
      pred_c_initial <- as.data.frame(predict(clf.case, newdata = test, type = "link"))
      pred_o_initial <- as.data.frame(predict(clf.target, newdata = test, type = "link"))
      pred_standarize <- pred_c_initial %>% mutate_all(~(scale(.) %>% as.vector))
      pred_c <- logit2prob(pred_standarize)
      pred_o_standarize <- pred_o_initial %>% mutate_all(~(scale(.) %>% as.vector))
      pred_o <- logit2prob(pred_o_standarize)
      
    }
    
    pred_f <- (1-w.target) * pred_c + w.target * pred_o
    pred <- data.frame(pred_f, 1-pred_f)
    names(pred) = c("TRUE","FALSE")
    
    auc = auroc(pred[ , "TRUE"], test$slides == "TRUE", plot = F)
    
    meta <- list(auc=auc, cw_target = c_weight)
    
    return(meta)
  }
  
  a <- c()
  
  # stratified sampling
  
  for (rep in c(1:5)) { # in our paper is 50 times
    
      resamp_F <- partition_kmeans(train_set[which(train_set$slides==F),],
                                   coords=c("x","y"), nfold = 10, seed1 = 123+rep)
      resamp_T <- partition_kmeans(train_set[which(train_set$slides==T),],
                                   coords=c("x","y"), nfold = 10, seed1 = 123+rep)
      
      ## construct the train set (70%) and the test set (30%)
      ss <- sample(1:10,3)
      sel.test_F <- c(resamp_F[[1]][[ss[1]]]$test,resamp_F[[1]][[ss[2]]]$test,
                      resamp_F[[1]][[ss[3]]]$test)
      sel.test_T <- c(resamp_T[[1]][[ss[1]]]$test,resamp_T[[1]][[ss[2]]]$test,
                      resamp_T[[1]][[ss[3]]]$test)
      
      train <- rbind(train_set[which(train_set$slides==F),][-sel.test_F,],
                     train_set[which(train_set$slides==T),][-sel.test_T,])
      
      test <- rbind(train_set[which(train_set$slides==F),][sel.test_F,],
                    train_set[which(train_set$slides==T),][sel.test_T,])
      
      ## construct the model trained by train
      clf.target <- mgcv::gam(f_target, data = train, method = "REML",family = "binomial")
                
      ## get the predictive accuracy of each c_weight
      a[[rep]] <- lapply(c_weight_range, tunning, clf.target = clf.target,
                       clf.case = clf.case, test = test,
                       batch_size=batch_size, epoch=epoch)
    
  }
  
  ## find the optimal c_weight
  cc <- a %>% bind_rows()
  mean.we <- sapply(c_weight_range, function(i) mean(cc$auc[which(cc$cw_target == i)]))
  optimal.cwe <- c_weight_range[which.max(mean.we)]
  
  return(optimal.cwe)
}
