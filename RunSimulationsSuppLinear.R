########
#### Stratified vs. Poisson sampling simulations - supplemental linear model
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
source("SimulationFunctionsSuppLinear.R")

######
### Preliminaries--------------
######

l2<-function(v) sqrt(sum(v*v))
tr<-function(M) sum(diag(M))

N <- 10000

######
### Scenario 1: Linear regression with 3 covariates with pilot samples, r1 = 200,
### n = 800
######
set.seed(1)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 800, nreps = 1000, r1 = 200,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1200n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1200n800, 
        "results/resultsThreeCovsLinearPilotR1200n800.rds")

######
### Scenario 2: Linear regression with 3 covariates with pilot samples, r1 = 200,
### n = 1200
######
set.seed(2)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 1200, nreps = 1000, r1 = 200,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1200n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1200n1200, 
        "results/resultsThreeCovsLinearPilotR1200n1200.rds")

######
### Scenario 3: Linear regression with 3 covariates with pilot samples, r1 = 200,
### n = 1600
######
set.seed(3)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 1600, nreps = 1000, r1 = 200,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1200n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1200n1600, 
        "results/resultsThreeCovsLinearPilotR1200n1600.rds")

######
### Scenario 4: Linear regression with 3 covariates with pilot samples, r1 = 600,
### n = 800
######
set.seed(4)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 800, nreps = 1000, r1 = 600,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1600n800 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1600n800, 
        "results/resultsThreeCovsLinearPilotR1600n800.rds")

######
### Scenario 5: Linear regression with 3 covariates with pilot samples, r1 = 600,
### n = 1200
######
set.seed(5)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 1200, nreps = 1000, r1 = 600,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1600n1200 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1600n1200, 
        "results/resultsThreeCovsLinearPilotR1600n1200.rds")

######
### Scenario 6: Linear regression with 3 covariates with pilot samples, r1 = 600,
### n = 1600
######
set.seed(6)

scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")

#cluster
ncores <- detectCores() - 1
cl <- makeCluster(ncores)
clusterExport(cl, c("ThreeCovsLinear", "N", "l2", "tr"))
parallel::clusterEvalQ(cl, {
  library(MASS)
  library(survey)
  library(dplyr)
  library(optimall)
})

# run sims
results_list <- parLapply(cl, scenarios, function(scen) {
  data.frame(
    Scenario = scen,
    ThreeCovsLinear(N = N, n = 1600, nreps = 1000, r1 = 600,
                    scenario = scen)
  )
})
stopCluster(cl)

# combine into one data.frame
resultsThreeCovsLinearPilotR1600n1600 <- do.call(rbind, results_list)
saveRDS(resultsThreeCovsLinearPilotR1600n1600, 
        "results/resultsThreeCovsLinearPilotR1600n1600.rds")
