######
######
###### Helper functions for simulations

## This version should be same as final, but has some extra strategies

######
###### Full data available

ThreeCovsLinear <- function(N, n, nreps = 1000, scenario, r1,
                      split_at = c(0.2, 0.8),
                      methods = c("SRS", "Case-control", "Stratified",
                                  "Optimal Poisson prob", "Case-control Surrogate",
                                  "Stratified with Pilot",
                                  "Error-prone Optimal Individual prob")){
  # Setup
  expit<-function(eta) exp(eta)/(1+exp(eta))
  logit<-function(p) log(p/(1-p))
  l2<-function(v) sqrt(sum(v*v))
  inf_fun_lm <- function(fit){
    dm <- model.matrix(fit)
    Ihat <- crossprod(dm) / nrow(dm)
    infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  }
  
  srs_B0 <- numeric(0)
  srs_B1 <- numeric(0)
  srs_B2 <- numeric(0)
  srs_B3 <- numeric(0)
  cc_B0 <- numeric(0)
  cc_B1 <- numeric(0)
  cc_B2 <- numeric(0)
  cc_B3 <- numeric(0)
  strat_B0 <- numeric(0)
  strat_B1 <- numeric(0)
  strat_B2 <- numeric(0)
  strat_B3 <- numeric(0)
  optp_B0 <- numeric(0)
  optp_B1 <- numeric(0)
  optp_B2 <- numeric(0)
  optp_B3 <- numeric(0)
  ccs_B0 <- numeric(0)
  ccs_B1 <- numeric(0)
  ccs_B2 <- numeric(0)
  ccs_B3 <- numeric(0)
  strats_B0 <- numeric(0)
  strats_B1 <- numeric(0)
  strats_B2 <- numeric(0)
  strats_B3 <- numeric(0)
  opt1_B0 <- numeric(0)
  opt1_B1 <- numeric(0)
  opt1_B2 <- numeric(0)
  opt1_B3 <- numeric(0)
  
  # Setting
  B0 <- ifelse(scenario == "Exp", -0.5, 0.5)
  B1 <- 0.5
  B2 <- 0.5
  B3 <- 0.5
  
  # Loop over nreps (include data generation in loop)
  for(i in 1:nreps){
    my_i <<- i
    set.seed(i)
    ###
    ### Generate data
    if (!(scenario %in% c("zeroMean", "unequalVar", "rareEvent",
                          "mixNormal", "T3", "Exp", "catX"))){
      stop("invalid scenario")
    } else if(scenario == "zeroMean"){
      pop <- data.frame(id = 1:N, mvrnorm(n = N, 
                                          mu = c(0, 0, 0),
                                          Sigma = matrix(c(1, 0.5, 0.5,
                                                           0.5, 1, 0.5,
                                                           0.5, 0.5, 1), nrow = 3)))
    } else if(scenario == "unequalVar"){
      pop <- data.frame(id = 1:N, mvrnorm(n = N, 
                                          mu = c(0, 0, 0),
                                          Sigma = matrix(c(1, 0.5/2, 0.5/3,
                                                           0.5/2, 1/4, 0.5/6,
                                                           0.5/3, 0.5/6, 1/9), nrow = 3)))
    } else if(scenario == "rareEvent"){
      pop <- data.frame(id = 1:N, mvrnorm(n = N, 
                                          mu = c(-1.6, -1.6, -1.6),
                                          Sigma = matrix(c(1, 0.5, 0.5,
                                                           0.5, 1, 0.5,
                                                           0.5, 0.5, 1), nrow = 3)))
    } else if(scenario == "mixNormal"){
      pop <- 0.5*data.frame(id = 1:N, mvrnorm(n = N, 
                                              mu = c(1, 1, 1),
                                              Sigma = matrix(c(1, 0.5, 0.5,
                                                               0.5, 1, 0.5,
                                                               0.5, 0.5, 1), nrow = 3))) + 
        0.5*data.frame(id = 1:N, mvrnorm(n = N, 
                                         mu = c(-1, -1, -1),
                                         Sigma = matrix(c(1, 0.5, 0.5,
                                                          0.5, 1, 0.5,
                                                          0.5, 0.5, 1), nrow = 3)))
    } else if(scenario == "T3"){
      pop <- data.frame(mvtnorm::rmvt(n = N, df = 3,
                                      sigma = matrix(c(1, 0.5, 0.5,
                                                       0.5, 1, 0.5,
                                                       0.5, 0.5, 1), nrow = 3)
      ))/10
      pop <- data.frame(id = 1:N, pop)
    } else if(scenario == "Exp"){
      pop <- data.frame(id = 1:N, X1 = rexp(N, rate = 2),
                        X2 = rexp(N, rate = 2),
                        X3 = rexp(N, rate = 2))
      B0 <- -0.5
    } else if(scenario == "catX"){
      pop <- data.frame(id = 1:N, mvrnorm(n = N, 
                                          mu = c(0, 0, 0),
                                          Sigma = matrix(c(1, 0.2, 0.5,
                                                           0.2, 1, 0.7,
                                                           0.5, 0.7, 1), nrow = 3)))
      pop$X1 <- ifelse(pop$X1 < quantile(pop$X1, 0.33), 1, 0)
      pop$X2 <- ifelse(pop$X2 < quantile(pop$X2, 0.2), 1, 0)
      pop$X3 <- ifelse(pop$X3 > quantile(pop$X3, 0.75), 1, 0)
    }
    
    pop$mu <- B0 + B1 * pop$X1 + B2 * pop$X2 + B3 * pop$X3
    pop$y <- pop$mu + rnorm(N)
    popmodel <- lm(y ~ X1 + X2 + X3, data=pop)
    
    ####
    #### Add measurement error
    pop$error_sd <- ifelse(pop$y < quantile(pop$y, 0.3), 1.5, 1)
    pop$y_obs <- pop$y + rnorm(N, sd = pop$error_sd)
    
    pop.raw <- pop
    
    # Compute error-prone pop model
    popmodel_obs <- lm(y_obs ~ X1 + X2 + X3, data=pop)
    
    ###
    ### Strategy 1 : SRS
    if ("SRS" %in% methods) {
    srs <- sample(pop$id,size = n, replace = FALSE)
    srs_data <- pop[pop$id %in% srs,]
    srsmodel <- lm(y ~ X1 + X2 + X3, data=srs_data)
    srs_B0[i] <- coef(srsmodel)[1]
    srs_B1[i] <- coef(srsmodel)[2]
    srs_B2[i] <- coef(srsmodel)[3]
    srs_B3[i] <- coef(srsmodel)[4]
    }
    
    ####
    #### Strategy 2: Case-control sample
    if ("Case-control" %in% methods) {
    
    ## A: Compute strata
    pop$strat <- ifelse(pop$y <= median(pop$y), 0, 1)
    
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(n/2, n/2), 
                         npop = c(table(pop$strat)))
    design$probs <- design$stratum_size/design$npop
  
    sampled_data <- optimall::sample_strata(pop, strata = "strat",
                                            id = "id", design_data = design,
                                            design_strata = "strat",
                                            probs = "probs",
                                            n_allocated = "stratum_size")
    cc_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strat, probs = ~sampling_prob, 
                           data = cc_data)
    ccmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    cc_B0[i] <- coef(ccmodel)[1]
    cc_B1[i] <- coef(ccmodel)[2]
    cc_B2[i] <- coef(ccmodel)[3]
    cc_B3[i] <- coef(ccmodel)[4]
    }
      
    ####
    #### Strategy 3: Stratified sample (A-optimal)
    if ("Stratified" %in% methods) {
    
    pop <- pop.raw
    
    ## A: Compute influence functions
    inflB0 <- inf_fun_lm(popmodel)[,1]
    pop$inflB0 <- inflB0
    inflB1 <- inf_fun_lm(popmodel)[,2]
    pop$inflB1 <- inflB1
    inflB2 <- inf_fun_lm(popmodel)[,3]
    pop$inflB2 <- inflB2
    inflB3 <- inf_fun_lm(popmodel)[,4]
    pop$inflB3 <- inflB3
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- ifelse(pop$y <= median(pop$y), 0, 1)
    if(scenario != "catX"){
      pop <- optimall::split_strata(pop,  strata = "strata",
                                    split_var = "inflB1",
                                    type = "global quantile",
                                    split_at = split_at) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB2",
                               type = "global quantile",
                               split_at = split_at) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB3",
                               type = "global quantile",
                               split_at = split_at) %>%
        dplyr::rename(strata = new_strata)
    } else {
      # for catX, only split on inflB2 and inflB3 since inflB1 is discrete
      pop$strata <- interaction(pop$y, pop$X1, pop$X2, pop$X3)
    }
    
    levels(pop$strata)<-1:length(table(pop$strata))
    
    ## Merge strata if any are size 1
    stratum_sizes <- table(pop$strata)
    small_strata <- names(stratum_sizes[stratum_sizes <= 1])
    
    if (length(small_strata) > 0) {
      # Get sizes of remaining strata (excluding size-1)
      larger_strata <- stratum_sizes[!(names(stratum_sizes) %in% small_strata)]
      # Identify the next smallest non-tiny stratum to merge into
      merge_target <- names(larger_strata)[which.min(larger_strata)]
      
      # Perform the merge: reassign tiny strata labels
      pop <- pop %>%
        dplyr::mutate(strata = ifelse(as.character(strata) %in% small_strata, 
                                      merge_target, as.character(strata))) %>%
        dplyr::mutate(strata = droplevels(factor(strata)))  # drop unused levels
    }
    nstrata <- length(unique(pop$strata))
    design <- optimum_allocation(pop, strata = "strata", 
                                 y = c("inflB0", "inflB1", "inflB2", "inflB3"),
                                 weights = c(0.25, 0.25, 0.25, 0.25), nsample = n - 2*nstrata, 
                                 method = "Neyman", ndigits = 10)
    design$stratum_size <- pmin(design$npop, design$stratum_size + 2)
    design$probs <- design$stratum_size/design$npop
    
    sampled_data <- optimall::sample_strata(pop, strata = "strata",
                                            id = "id", design_data = design,
                                            probs = "probs",
                                            n_allocated = "stratum_size")
    strat_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strata, probs = ~sampling_prob, 
                           data = strat_data)
    stratmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    strat_B0[i] <- coef(stratmodel)[1]
    strat_B1[i] <- coef(stratmodel)[2]
    strat_B2[i] <- coef(stratmodel)[3]
    strat_B3[i] <- coef(stratmodel)[4]
    }
    
    ####
    ### Strategy 4: Optimal Poisson sampling (uses full data MLE)
    if ("Optimal Poisson prob" %in% methods) {
    
    pop <- pop.raw
    
    ## A: Compute preliminaries
    Mx.inv<- vcov(popmodel)*N
    X <- model.matrix(popmodel)
    Mx_inv_x <- X %*% t(Mx.inv)
    pop$Mx_inv_x_norm <- sqrt(rowSums(Mx_inv_x^2))
    # or equivalently
    # Ihat <- (t(dm) %*% (dm * popmodel$fitted.values * (1 - popmodel$fitted.values))) /nrow(dm)
    # Mx.in <- solve(Ihat)
    pop$y.hat <- fitted(popmodel)
    #pop$optp_raw <- with(pop, abs(y-y.hat)*l2(t(Mx.inv%*%rbind(1,pop$X1, pop$X2, pop$X3))))
    pop$optp_raw <- with(pop, abs(y - y.hat)*Mx_inv_x_norm)
    pop$optp <- n*pop$optp_raw / sum(pop$optp_raw)
    pop$optp <- pmin(pop$optp, 1) # cap at 1
    pop$wt<- sum(pop$optp)/n/pop$optp
    
    # in_sample<-sample(1:nrow(pop),1000, prob=pop$optp)
    # pop$optp <- pop$optp / sum(pop$optp)*n # already did this above
    
    # Poisson: This is actually slightly better than with replacement
    in_sample <- rbinom(N, 1, pop$optp)
    the_sample<-pop[in_sample == 1,]
    
    # with replacement
    # in_sample <- sample(seq(1,N,1), n, replace=TRUE, prob=pop$optp)
    # the_sample <- pop[in_sample,]
    
    svydesign <- svydesign(id = ~id, probs = ~optp, 
                           data = the_sample)
    optmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    optp_B0[i] <- coef(optmodel)[1]
    optp_B1[i] <- coef(optmodel)[2]
    optp_B2[i] <- coef(optmodel)[3]
    optp_B3[i] <- coef(optmodel)[4]
    }
     
    ####
    #### Strategy 5: Case-control sampling using surrogate 
    if ("Case-control Surrogate" %in% methods) {
    
    pop <- pop.raw
    
    ## A: Compute strata
    pop$strat <- ifelse(pop$y_obs <= median(pop$y_obs), 0, 1)
    
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(n/2, n/2), 
                         npop = c(table(pop$strat)))
    design$probs <- design$stratum_size/design$npop
    sampled_data <- optimall::sample_strata(pop, strata = "strat",
                                            id = "id", design_data = design,
                                            design_strata = "strat",
                                            probs = "probs",
                                            n_allocated = "stratum_size")
    ccs_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strat, probs = ~sampling_prob, 
                           data = ccs_data)
    ccsmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    ccs_B0[i] <- coef(ccsmodel)[1]
    ccs_B1[i] <- coef(ccsmodel)[2]
    ccs_B2[i] <- coef(ccsmodel)[3]
    ccs_B3[i] <- coef(ccsmodel)[4] 
    }
    
    ####
    #### Strategy 6: Stratified sampling using pilot
    if ("Stratified with Pilot" %in% methods) {
    
    pop <- pop.raw
    
    ## A: Compute influence functions from error-prone model
    inflB0_obs <- inf_fun_lm(popmodel_obs)[,1]
    pop$inflB0_obs <- inflB0_obs
    inflB1_obs <- inf_fun_lm(popmodel_obs)[,2]
    pop$inflB1_obs <- inflB1_obs
    inflB2_obs <- inf_fun_lm(popmodel_obs)[,3]
    pop$inflB2_obs <- inflB2_obs
    inflB3_obs <- inf_fun_lm(popmodel_obs)[,4]
    pop$inflB3_obs <- inflB3_obs
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- ifelse(pop$y_obs <= median(pop$y_obs), 0, 1)
    pop <- optimall::split_strata(pop,  strata = "strata",
                                  split_var = "inflB1_obs",
                                  type = "global quantile",
                                  split_at = split_at) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB2_obs",
                             type = "global quantile",
                             split_at = split_at) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB3_obs",
                             type = "global quantile",
                             split_at = split_at) %>%
      dplyr::rename(strata = new_strata)
    
    levels(pop$strata)<-1:length(table(pop$strata))
    
    ## Merge strata if any are size <= 1
    stratum_sizes <- table(pop$strata)
    small_strata <- names(stratum_sizes[stratum_sizes <= 1])
    
    if (length(small_strata) > 0) {
      # Get sizes of remaining strata (excluding size-1)
      larger_strata <- stratum_sizes[!(names(stratum_sizes) %in% small_strata)]
      # Identify the next smallest non-tiny stratum to merge into
      merge_target <- names(larger_strata)[which.min(larger_strata)]
      
      # Perform the merge: reassign tiny strata labels
      pop <- pop %>%
        dplyr::mutate(strata = ifelse(as.character(strata) %in% small_strata, merge_target,
                                      as.character(strata))) %>%
        dplyr::mutate(strata = droplevels(factor(strata)))  # drop unused levels
    }
    nstrata <- length(unique(pop$strata))
    design <- optimum_allocation(pop, strata = "strata", y = c("inflB0_obs", 
                                                               "inflB1_obs", 
                                                               "inflB2_obs", 
                                                               "inflB3_obs"),
                                 weights = c(0.25, 0.25, 0.25, 0.25), 
                                 nsample = r1 - 2*nstrata, 
                                 method = "Neyman", ndigits = 10)
    design$stratum_size <- pmin(design$npop, design$stratum_size + 2)
    
    design$probs1 <- design$stratum_size/design$npop
    
    ## Pilot
    sampled_data <- optimall::sample_strata(pop, strata = "strata",
                                            id = "id", design_data = design,
                                            probs = "probs1",
                                            n_allocated = "stratum_size")
    sampled_data$sampled_pilot <- sampled_data$sample_indicator
    strat_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strata, probs = ~sampling_prob, 
                           data = strat_data)
    stratmodel1 <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    
    ## Now use results from this to improve optimum allocation
    pilot_infs <- data.frame(inf_fun_lm(stratmodel1))
    names(pilot_infs) <- c("inflB0_pilot", "inflB1_pilot", "inflB2_pilot",
                           "inflB3_pilot")
    pilot_infs$id <- as.numeric(row.names(pilot_infs))
    pop_new <- dplyr::left_join(sampled_data, pilot_infs, by = "id")
    design2 <- optimall::allocate_wave(pop_new, strata = "strata", y = c("inflB0_pilot", 
                                                                         "inflB1_pilot", 
                                                                         "inflB2_pilot", 
                                                                         "inflB3_pilot"),
                                       weights = c(0.25, 0.25, 0.25, 0.25), nsample = n - r1, 
                                       allocation_method = "Neyman", 
                                       already_sampled = "sample_indicator")
    design2$n_to_sample <- pmin(design2$npop - design2$nsample_prior, design2$n_to_sample)
    design2$probs2 <- (design2$nsample_prior + design2$n_to_sample)/design$npop
    
    sampled_data2 <- optimall::sample_strata(sampled_data, strata = "strata",
                                             id = "id", design_data = design2,
                                             already_sampled = "sampled_pilot",
                                             n_allocated = "n_to_sample")
    strat_data2 <- sampled_data2[sampled_data2$sample_indicator == 1 |
                                   sampled_data2$sampled_pilot == 1,]
    strat_data2 <- dplyr::left_join(strat_data2, design2[,c("strata", "probs2")],
                                    by = "strata")
    svydesign <- svydesign(id = ~id, strata = ~strata, probs = ~probs2, 
                           data = strat_data2)
    stratmodel2 <- svyglm(y ~ X1 + X2 + X3, design = svydesign)
    
    strats_B0[i] <- coef(stratmodel2)[1]
    strats_B1[i] <- coef(stratmodel2)[2]
    strats_B2[i] <- coef(stratmodel2)[3]
    strats_B3[i] <- coef(stratmodel2)[4]
    }
    
    ####
    #### Strategy 7: Error-prone optimal individual probability
    if ("Error-prone Optimal Individual prob" %in% methods) {
    
    pop <- pop.raw
    
    ## A: Compute preliminaries from the error-prone full-data model.
    ## Match the full-data optimal Poisson score, replacing y with y_obs.
    Mx.obs.inv <- vcov(popmodel_obs)*N
    X.obs <- model.matrix(popmodel_obs)
    Mx_obs_inv_x <- X.obs %*% t(Mx.obs.inv)
    pop$Mx_obs_inv_x_norm <- sqrt(rowSums(Mx_obs_inv_x^2))
    pop$y_obs_hat <- fitted(popmodel_obs)
    pop$optp_raw <- with(pop, abs(y_obs - y_obs_hat) * Mx_obs_inv_x_norm)
    pop$optp <- n * pop$optp_raw / sum(pop$optp_raw)
    pop$optp <- pmin(pop$optp, 1) # cap at 1

    in_sample <- rbinom(N, 1, pop$optp)
    the_sample <- pop[in_sample == 1,]

    svydesign <- svydesign(id = ~id, probs = ~optp, data = the_sample)
    optmodel1 <- svyglm(y ~ X1 + X2 + X3, design = svydesign)

    opt1_B0[i] <- coef(optmodel1)[1]
    opt1_B1[i] <- coef(optmodel1)[2]
    opt1_B2[i] <- coef(optmodel1)[3]
    opt1_B3[i] <- coef(optmodel1)[4]
    }
      
  }
    
  # Gather results
  results <- data.frame(
    method = "SRS",
    mean_B0 = mean(srs_B0),
    mean_B1 = mean(srs_B1),
    mean_B2 = mean(srs_B2),
    mean_B3 = mean(srs_B3),
    var_B0 = var(srs_B0),
    var_B1 = var(srs_B1),
    var_B2 = var(srs_B2), 
    var_B3 = var(srs_B3), 
    MSE_B0 = mean((srs_B0 - B0)^2),
    MSE_B1 = mean((srs_B1 - B1)^2),
    MSE_B2 = mean((srs_B2 - B2)^2),
    MSE_B3 = mean((srs_B3 - B3)^2),
    MSE_B = mean((srs_B0 - B0)^2 + (srs_B1 - B1)^2 + (srs_B2 - B2)^2 + (srs_B3 - B3)^2)
  )
  results <- rbind(results, data.frame(
    method = "Case-control",
    mean_B0 = mean(cc_B0),
    mean_B1 = mean(cc_B1),
    mean_B2 = mean(cc_B2),
    mean_B3 = mean(cc_B3),
    var_B0 = var(cc_B0),
    var_B1 = var(cc_B1),
    var_B2 = var(cc_B2),
    var_B3 = var(cc_B3),
    MSE_B0 = mean((cc_B0 - B0)^2),
    MSE_B1 = mean((cc_B1 - B1)^2),
    MSE_B2 = mean((cc_B2 - B2)^2),
    MSE_B3 = mean((cc_B3 - B3)^2),
    MSE_B = mean((cc_B0 - B0)^2 + (cc_B1 - B1)^2 + (cc_B2 - B2)^2 + (cc_B3 - B3)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Stratified",
    mean_B0 = mean(strat_B0),
    mean_B1 = mean(strat_B1),
    mean_B2 = mean(strat_B2),
    mean_B3 = mean(strat_B3),
    var_B0 = var(strat_B0),
    var_B1 = var(strat_B1),
    var_B2 = var(strat_B2),
    var_B3 = var(strat_B3),
    MSE_B0 = mean((strat_B0 - B0)^2),
    MSE_B1 = mean((strat_B1 - B1)^2),
    MSE_B2 = mean((strat_B2 - B2)^2),
    MSE_B3 = mean((strat_B3 - B3)^2),
    MSE_B = mean((strat_B0 - B0)^2 + (strat_B1 - B1)^2 + (strat_B2 - B2)^2 + (strat_B3 - B3)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Optimal Poisson prob",
    mean_B0 = mean(optp_B0),
    mean_B1 = mean(optp_B1),
    mean_B2 = mean(optp_B2),
    mean_B3 = mean(optp_B3),
    var_B0 = var(optp_B0),
    var_B1 = var(optp_B1),
    var_B2 = var(optp_B2),
    var_B3 = var(optp_B3),
    MSE_B0 = mean((optp_B0 - B0)^2),
    MSE_B1 = mean((optp_B1 - B1)^2),
    MSE_B2 = mean((optp_B2 - B2)^2),
    MSE_B3 = mean((optp_B3 - B3)^2),
    MSE_B = mean((optp_B0 - B0)^2 + (optp_B1 - B1)^2 + (optp_B2 - B2)^2 + (optp_B3 - B3)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Case-control Surrogate",
    mean_B0 = mean(ccs_B0),
    mean_B1 = mean(ccs_B1),
    mean_B2 = mean(ccs_B2),
    mean_B3 = mean(ccs_B3),
    var_B0 = var(ccs_B0),
    var_B1 = var(ccs_B1),
    var_B2 = var(ccs_B2),
    var_B3 = var(ccs_B3),
    MSE_B0 = mean((ccs_B0 - B0)^2),
    MSE_B1 = mean((ccs_B1 - B1)^2),
    MSE_B2 = mean((ccs_B2 - B2)^2),
    MSE_B3 = mean((ccs_B3 - B3)^2),
    MSE_B = mean((ccs_B0 - B0)^2 + (ccs_B1 - B1)^2 + (ccs_B2 - B2)^2 + (ccs_B3 - B3)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Stratified with Pilot",
    mean_B0 = mean(strats_B0),
    mean_B1 = mean(strats_B1),
    mean_B2 = mean(strats_B2),
    mean_B3 = mean(strats_B3),
    var_B0 = var(strats_B0),
    var_B1 = var(strats_B1),
    var_B2 = var(strats_B2),
    var_B3 = var(strats_B3),
    MSE_B0 = mean((strats_B0 - B0)^2),
    MSE_B1 = mean((strats_B1 - B1)^2),
    MSE_B2 = mean((strats_B2 - B2)^2),
    MSE_B3 = mean((strats_B3 - B3)^2),
    MSE_B = mean((strats_B0 - B0)^2 + (strats_B1 - B1)^2 + (strats_B2 - B2)^2 + (strats_B3 - B3)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Error-prone Optimal Individual prob",
    mean_B0 = mean(opt1_B0),
    mean_B1 = mean(opt1_B1),
    mean_B2 = mean(opt1_B2),
    mean_B3 = mean(opt1_B3),
    var_B0 = var(opt1_B0),
    var_B1 = var(opt1_B1),
    var_B2 = var(opt1_B2),
    var_B3 = var(opt1_B3),
    MSE_B0 = mean((opt1_B0 - B0)^2),
    MSE_B1 = mean((opt1_B1 - B1)^2),
    MSE_B2 = mean((opt1_B2 - B2)^2),
    MSE_B3 = mean((opt1_B3 - B3)^2),
    MSE_B = mean((opt1_B0 - B0)^2 + (opt1_B1 - B1)^2 + (opt1_B2 - B2)^2 + (opt1_B3 - B3)^2)
  ))
  results$varsum <- results$var_B0 + results$var_B1 + results$var_B2 + results$var_B3
  results$MSEsum <- results$MSE_B0 + results$MSE_B1 + results$MSE_B2 + results$MSE_B3
  results <- results[results$method %in% methods, , drop = FALSE]
  return(results)
}
