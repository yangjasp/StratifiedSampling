########
#### Stratified vs. Poisson sampling simulatioms
########

#### Author: Jasper Yang

####
## Load packages
####
library(dplyr)
library(tidyr)
library(optimall)
library(MASS)
library(survey)
library(parallel)
library(rprojroot)
library(optparse)

#####
## Source functions
#####
source("SimulationFunctions.R")
source("ossat_functions.R")

######
### Preliminaries--------------
######

# And little utils here
expit<-function(eta) exp(eta)/(1+exp(eta))
logit<-function(p) log(p/(1-p))
l2<-function(v) sqrt(sum(v*v))
tr<-function(M) sum(diag(M))

inf_fun_logit <- function(fit){
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) /nrow(dm)
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
}

######
### Scenario 1: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 200,
### n = 800
######
set.seed(1)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 800, nreps = 1000, r1 = 200,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1200n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1200n800, 
        "results/resultsThreeCovsPilotLowErrorR1200n800.rds")

######
### Scenario 2: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 200,
### n = 1200
######
set.seed(2)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1200, nreps = 1000, r1 = 200,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1200n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1200n1200, 
        "results/resultsThreeCovsPilotLowErrorR1200n1200.rds")

######
### Scenario 3: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 200,
### n = 1600
######
set.seed(3)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1600, nreps = 1000, r1 = 200,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1200n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1200n1600, 
        "results/resultsThreeCovsPilotLowErrorR1200n1600.rds")

######
### Scenario 4: Logistic regression with 3 covariates with pilot samples - High error, r1 = 200,
### n = 800
######
set.seed(4)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 800, nreps = 1000, r1 = 200,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1200n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1200n800, 
        "results/resultsThreeCovsPilotHighErrorR1200n800.rds")

######
### Scenario 5: Logistic regression with 3 covariates with pilot samples - High error, r1 = 200,
### n = 1200
######
set.seed(5)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1200, nreps = 1000, r1 = 200,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1200n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1200n1200, 
        "results/resultsThreeCovsPilotHighErrorR1200n1200.rds")

######
### Scenario 6: Logistic regression with 3 covariates with pilot samples - High error, r1 = 200,
### n = 1600
######
set.seed(6)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1600, nreps = 1000, r1 = 200,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1200n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1200n1600, 
        "results/resultsThreeCovsPilotHighErrorR1200n1600.rds")

######
### Scenarios 7-9: Data generation outside of loop w/ catX - low error, r1 = 200,
### n = 800, 1200, 1600
######
set.seed(7)

my_ns <- c(800,1200,1600)

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsDataOutside", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4", "my_ns"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, my_ns, function(my_n) {
  data.frame(
    n = my_n,
    ThreeCovsDataOutside(N = N, n = my_n, nreps = 1000, r1 = 200,
                         error = "low", scenario = "catX")
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsCatXPilotLowErrorR1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsCatXPilotLowErrorR1200, 
        "results/resultsThreeCovsCatXPilotLowErrorR1200.rds")

######
### Scenarios 10-12: Data generation outside of loop w/ catX - high error, r1 = 200,
### n = 800, 1200, 1600
######
set.seed(10)

my_ns <- c(800,1200,1600)

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsDataOutside", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4", "my_ns"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, my_ns, function(my_n) {
  data.frame(
    n = my_n,
    ThreeCovsDataOutside(N = N, n = my_n, nreps = 1000, r1 = 200,
                         error = "high", scenario = "catX")
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsCatXPilotHighErrorR1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsCatXPilotHighErrorR1200, 
        "results/resultsThreeCovsCatXPilotHighErrorR1200.rds")

######
### Scenario 13: Data generation outside of loop w/ catX - low error, r1 = 600,
### n = 800, 1200, 1600
######
set.seed(13)

my_ns <- c(800,1200,1600)

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsDataOutside", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4", "my_ns"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, my_ns, function(my_n) {
  data.frame(
    n = my_n,
    ThreeCovsDataOutside(N = N, n = my_n, nreps = 1000, r1 = 600,
                         error = "low", scenario = "catX")
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsCatXPilotLowErrorR1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsCatXPilotLowErrorR1600, 
        "results/resultsThreeCovsCatXPilotLowErrorR1600.rds")

######
### Scenario 14: Data generation outside of loop w/ catX - high error, r1 = 600,
### n = 800, 1200, 1600
######
set.seed(14)

my_ns <- c(800,1200,1600)

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsDataOutside", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4", "my_ns"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, my_ns, function(my_n) {
  data.frame(
    n = my_n,
    ThreeCovsDataOutside(N = N, n = my_n, nreps = 1000, r1 = 600,
                         error = "high", scenario = "catX")
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsCatXPilotHighErrorR1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsCatXPilotHighErrorR1600, 
        "results/resultsThreeCovsCatXPilotHighErrorR1600.rds")


#####
#####
##### Start supplement

######
### Scenario 15: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 600,
### n = 800
######
set.seed(15)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 800, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1600n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1600n800, 
        "results/resultsThreeCovsPilotLowErrorR1600n800.rds")

######
### Scenario 16: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 600,
### n = 1200
######
set.seed(16)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1200, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1600n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1600n1200, 
        "results/resultsThreeCovsPilotLowErrorR1600n1200.rds")

######
### Scenario 17: Logistic regression with 3 covariates with pilot samples - Low error, r1 = 600,
### n = 1600
######
set.seed(17)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1600, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotLowErrorR1600n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotLowErrorR1600n1600, 
        "results/resultsThreeCovsPilotLowErrorR1600n1600.rds")

######
### Scenario 18: Logistic regression with 3 covariates with pilot samples - High error, r1 = 600,
### n = 800
######
set.seed(18)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 800, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1600n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1600n800, 
        "results/resultsThreeCovsPilotHighErrorR1600n800.rds")

######
### Scenario 19: Logistic regression with 3 covariates with pilot samples - High error, r1 = 600,
### n = 1200
######
set.seed(19)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1200, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1600n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1600n1200, 
        "results/resultsThreeCovsPilotHighErrorR1600n1200.rds")

######
### Scenario 20: Logistic regression with 3 covariates with pilot samples - High error, r1 = 600,
### n = 1600
######
set.seed(20)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovs(N = N, n = 1600, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsPilotHighErrorR1600n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsPilotHighErrorR1600n1600, 
        "results/resultsThreeCovsPilotHighErrorR1600n1600.rds")

######
### Scenario 21: Logistic regression with 7 covariates with pilot samples - Low error, r1 = 600,
### n = 800
######
set.seed(21)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 800, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotLowErrorR1600n800 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotLowErrorR1600n800, 
        "results/resultsSevenCovsPilotLowErrorR1600n800.rds")

######
### Scenario 22: Logistic regression with 3 covariates with pilot samples - High error, r1 = 600,
### n = 800
######
set.seed(22)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 800, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotHighErrorR1600n800 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotHighErrorR1600n800, 
        "results/resultsSevenCovsPilotHighErrorR1600n800.rds")

######
### Scenario 23: Logistic regression with 7 covariates with pilot samples - Low error, r1 = 600,
### n = 1200
######
set.seed(23)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 1200, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotLowErrorR1600n1200 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotLowErrorR1600n1200, 
        "results/resultsSevenCovsPilotLowErrorR1600n1200.rds")

######
### Scenario 24: Logistic regression with 3 covariates with pilot samples - High error, r1 = 600,
### n = 1200
######
set.seed(24)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 1200, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotHighErrorR1600n1200 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotHighErrorR1600n1200, 
        "results/resultsSevenCovsPilotHighErrorR1600n1200.rds")

######
### Scenario 25: Logistic regression with 7 covariates with pilot samples - Low error, r1 = 600,
### n = 1600
######
set.seed(25)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 1600, nreps = 1000, r1 = 600,
              error = "low", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotLowErrorR1600n1600 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotLowErrorR1600n1600, 
        "results/resultsSevenCovsPilotLowErrorR1600n1600.rds")

######
### Scenario 26: Logistic regression with 7 covariates with pilot samples - High error, r1 = 600,
### n = 1600
######
set.seed(26)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("SevenCovs", "N", "expit", "l2", 
                    "inf_fun_logit", "logit", "tr",
                    "getMLE.firth", "weighted.model.set1" ,
                    "weighted.model.seq4"
))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
  library(logistf)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    SevenCovs(N = N, n = 1600, nreps = 1000, r1 = 600,
              error = "high", scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsSevenCovsPilotHighErrorR1600n1600 <- do.call(rbind, results_list)
saveRDS(resultsSevenCovsPilotHighErrorR1600n1600, 
        "results/resultsSevenCovsPilotHighErrorR1600n1600.rds")
