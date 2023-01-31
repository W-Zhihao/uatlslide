
# Unsupervised active-transfer learning for automated landslide mapping

Zhihao Wang

## Preface

This vignette walks you through most of the analyses performed for the
paper that introduces strategies used in paper. Please refer to that
paper for conceptual and formal details, and cite it when using or
referring to the methods and results presented herein.

> Wang, Z. H., Brenning, A. (submitted) Unsupervised active-transfer
> learning for automated landslide mapping. Submitted to Computers and
> Geosciences, 2023.

Here is the example to elaborate how to select "informative" data from
the target area and how we combine the model trained by "suitable"
existing source area with the model trained by selected data from the
target area.

## Getting started

### Work environment

Make sure that all required packages and their dependencies are
installed and up-to-date. In addition, you will need the packages the
packages `mgcv`, `ROCR`,`dplyr`, `sperrorest` for running this new
framework based on r.

### Case study and data preparation

We use Palu as the target area and Reuleut as the source area.

data_1.rds are from the study area named Palu.
data_2.rds are from the study area named Reuleut.

Note: these data are subset of whole data sets for quick going through
methods used. Please find all data sets in Data availability in paper.

Let's get started by preparing the data set:

``` r

data1 <- readRDS("~/data_1.rds") # Palu
data2 <- readRDS("~/data_2.rds") # Reuleut


summary(data1)
summary(data2)

```

### Case-based reasoning

Case-based reasoning method:

    Similarity_functions.R ------- calculating the similarity between
    the source area and the target area

``` r
# calculating overall similarities

source("~/Similarity_functions.R") # similarity functions of all attributes

attributes <- readRDS("~/simfactors.rds") # load attributes of two study areas

# > attributes$reuleut
# mean slope      std slope           igneous        sedimentary        metamorphic        
# 15.1              9.5                1.0                1.0                1.0    
# total relief    spatial resolution 
# 1849.0            8.0 

# > attributes$palu
# mean slope      std slope            igneous       sedimentary       metamorphic 
# 21.8              12.4                1.0               1.0                1.0 
# total relief    spatial resolution 
# 2100.0            8.0 

## the target area: palu

### the source area: reuleut

sim_meanslope <- MEANSLOPE_SIMILARITY_FUN(attributes$reuleut, attributes$palu, 60) 
# result: 0.825
sim_resolution <- RESOLUTION_SIMILARITY_FUN(attributes$reuleut, attributes$palu)
# result: 1
sim_totalrelief <- TOTALRELIEF_SIMILARITY_FUN(attributes$reuleut, attributes$palu)
# result: 0.963
sim_stdslope <- STDSLOPE_SIMILARITY_FUN(attributes$reuleut, attributes$palu)
# result: 0.716
sim_geounits <- GEOLOGICALUNITS_SIMILARITY_FUN(attributes$reuleut, attributes$palu)
# result: 1

overall_sim <- min(sim_meanslope,sim_resolution,sim_totalrelief,sim_stdslope,sim_geounits)
# result: 0.72

```

The CBR model:

 The model trained by data2


``` r

# Formula of source area
f_source <- fo <- as.formula(slides ~ s(slope, k = 6)+s(plancurv, k =6) +
                  s(profcurv, k = 6) + s(log.carea, k = 6) + s(cslope, k = 6) + 
                   s(TWI, k = 6)+ s(dem, k = 6))
cbr_model <- mgcv::gam(f_source, data = data2, method = "REML",family = "binomial")

```

### Initial training set selected

``` r

source('~/low-prevelance.R')

# calculate the posterior probabilities of the target area

pred_c_initial <- as.data.frame(predict(cbr_model, newdata = data1, type = "response"))
pred_c <- data.frame(pred_c_initial, 1-pred_c_initial)
names(pred_c) = c("TRUE","FALSE")

select_loc <- low_pre_MS(pred_c, 100, 0.9)
initial_ts <- data1[select_loc,]

# > summary(initial_ts)
#        x                  y                dem             carea               cslope           slope     
#  Min.   :13349114   Min.   :-168744   Min.   : 508.4   Min.   :    71.76   Min.   : 7.461   Min.   : 2.844
#  1st Qu.:13357592   1st Qu.:-154974   1st Qu.: 865.7   1st Qu.:   426.83   1st Qu.:21.841   1st Qu.:15.508
#  Median :13359049   Median :-147317   Median : 918.4   Median :   890.90   Median :32.298   Median :21.340
#  Mean   :13362412   Mean   :-149021   Mean   : 909.0   Mean   :  6323.30   Mean   :29.467   Mean   :23.890
#  3rd Qu.:13362266   3rd Qu.:-144856   3rd Qu.:1005.7   3rd Qu.:  3512.43   3rd Qu.:38.527   3rd Qu.:35.021
#  Max.   :13388203   Max.   :-127543   Max.   :1247.8   Max.   :105664.23   Max.   :39.385   Max.   :47.521  
#     plancurv            profcurv              TWI           slides     log.carea    
#  Min.   :-0.134790   Min.   :-0.019090   Min.   : 5.708   FALSE:41   Min.   :1.856  
#  1st Qu.:-0.022725   1st Qu.:-0.002750   1st Qu.: 6.662   TRUE :59   1st Qu.:2.630  
#  Median :-0.001705   Median : 0.002690   Median : 7.427              Median :2.949  
#  Mean   :-0.007396   Mean   : 0.002161   Mean   : 7.840              Mean   :3.098  
#  3rd Qu.: 0.015995   3rd Qu.: 0.007627   3rd Qu.: 8.759              3rd Qu.:3.546  
#  Max.   : 0.127290   Max.   : 0.017290   Max.   :12.374              Max.   :5.024  

# remove these data from the data1

unlabel_set <- data1[-select_loc,]
    
```

### Adaptive unsupervised active-transfer learning

Here, we show the example for "informative" data selection from the target area.
Repeating the following steps can get more "informative" data from the target area.


Step one: getting the optimal c_weight value for AL model in epoch 1.
See Adaptive unsupervised active-transfer learning.R for more details

``` r
source("~/Adaptive unsupervised active-transfer learning.R")

# inputs 
c_weight_range <- c(0.005,0.007)
# Formula for target area:
f_target <- as.formula(slides ~ s(slope, k = 6)+s(plancurv, k =6) +
                  s(profcurv, k = 6) + s(log.carea, k = 6) + s(cslope, k = 6) + 
                   s(TWI, k = 6)+ s(dem, k = 6))
batch_size <- 25
epoch <- 1

# output the optimal c_weight
opti.cwe <- adaptive_UAT(c_weight_range, f_target, initial_ts, cbr_model, batch_size, epoch)
# 0.007
    
```

Step two: getting the optimal weight value for AL and CBR models in epoch 1.
See weight_function.R for more details

``` r
source("~/weight_function.R")

weight_al <- weight_AL_function(opti.cwe, batch_size, epoch)
weight_cbr <- 1 - weight_al
    
```
Step three: selecting the "informative" data set from the unlabeled set.

``` r

# calculate the posterior probabilities of the target area based on AL and CBR models

## change logit to probability

  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }

clf.target <- mgcv::gam(f_target, data = initial_ts, method = "REML",family = "binomial") 

## standardize and combine

pred_c_initial <- as.data.frame(predict(cbr_model, newdata = unlabel_set, type = "link"))
pred_o_initial <- as.data.frame(predict(clf.target, newdata = unlabel_set, type = "link"))
pred_standarize <- pred_c_initial %>% mutate_all(~(scale(.) %>% as.vector))
pred_c <- logit2prob(pred_standarize) # CBR model
pred_o_standarize <- pred_o_initial %>% mutate_all(~(scale(.) %>% as.vector))
pred_o <- logit2prob(pred_o_standarize) # AL model

pred_f <- weight_cbr*pred_c + weight_al*pred_o
pred <- data.frame(pred_f, 1-pred_f)
names(pred) = c("TRUE","FALSE")

## new selection

select_loc <- low_pre_MS(pred, 25, 0.9)
initial_ts <- rbind(initial_ts, unlabel_set[select_loc,])
unlabel_set <- unlabel_set[-select_loc,]
    
```

### Regular unsupervised active-transfer learning

This method is to use the c_weight value obtained in the first epoch (e.g., 0.007)
for the rest of epochs.

The process of "informative" data selection is similar with adaptive 
unsupervised active-transfer learning, except the optimal c_weight selection.
