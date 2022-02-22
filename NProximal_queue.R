library(caret)
library(corpcor)
library(rdetools)
library(expm)
library(MASS)
source("DGP2.R")
source("NEval.R")

dim_x <- 5
dim_z <- 2
dim_w <- 2
type <- "parametric"
kernel_sigma_list <- c(5, 15, 35, 60, 100, 150)  # default 15
CV_K_fold = 5
CF_K_fold = 5



N <- 400





## candidates of hyperparameters
lm_Hh1_list <- c(0.0001, 0.001, 0.01, 0.1)
lm_Qh1_list <- c(0.1, 1, 10, 100)
lm_Hh0_list <- c(0.0001, 0.001, 0.01, 0.1)
lm_Qh0_list <- c(0.1, 1, 10, 100)
lm_Hq1_list <- c(0.1, 1, 10, 100)
lm_Qq1_list <- c(0.0001, 0.001, 0.01, 0.1)
lm_Hq0_list <- c(0.1, 1, 10, 100)
lm_Qq0_list <- c(0.0001, 0.001, 0.01, 0.1)



# old simulation: lm_Qq <- 0.0025/(N^.8), lm_Hq <- 25/(N^.8), no n^2 in alpha
data <- DGP(N, dim_x, dim_z, dim_w, type)

## use cross validation to determine the hyperparamter
result <- CV_eval(data, kernel_sigma_list, theta = 0, lm_Hh1_list, lm_Qh1_list, lm_Hh0_list, lm_Qh0_list,
                           lm_Hq1_list, lm_Qq1_list, lm_Hq0_list, lm_Qq0_list, CV_K_fold, Nystroem = TRUE, appro_rate = 0.05)

kernel_sigma <- result$kernel_sigma
hyperparameters <- result$hyperparameters



lm_Hh1 <- hyperparameters[1]
lm_Qh1 <- hyperparameters[2]
lm_Hh0 <- hyperparameters[3]
lm_Qh0 <- hyperparameters[4]
lm_Hq1 <- hyperparameters[5]
lm_Qq1 <- hyperparameters[6]
lm_Hq0 <- hyperparameters[7]
lm_Qq0 <- hyperparameters[8]
## use cross fitting to get the empirical bias
est <- CF_eval(data, kernel_sigma, lm_Hh1, lm_Qh1, lm_Hh0, lm_Qh0, lm_Hq1, lm_Qq1, lm_Hq0, lm_Qq0, CF_K_fold, Nystroem = TRUE, appro_rate = 0.05)
est <- Re(est)

write.csv(est, file= paste0("Res/", type, "_", N, "_", njob,".csv"), row.names = FALSE)






