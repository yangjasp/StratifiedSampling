######
######
###### Helper functions for simulations

## This version should be same as final, but has some extra strategies

######
###### Full data available

ThreeCovs <- function(N, n, nreps = 1000, scenario, r1, error = "low"){
  # Setup
  expit<-function(eta) exp(eta)/(1+exp(eta))
  logit<-function(p) log(p/(1-p))
  l2<-function(v) sqrt(sum(v*v))
  
  srs_B0 <- c()
  srs_B1 <- c()
  srs_B2 <- c()
  srs_B3 <- c()
  cc_B0 <- c()
  cc_B1 <- c()
  cc_B2 <- c()
  cc_B3 <- c()
  strat_B0 <- c()
  strat_B1 <- c()
  strat_B2 <- c()
  strat_B3 <- c()
  optp_B0 <- c()
  optp_B1 <- c()
  optp_B2 <- c()
  optp_B3 <- c()
  ccs_B0 <- c()
  ccs_B1 <- c()
  ccs_B2 <- c()
  ccs_B3 <- c()
  strats_B0 <- c()
  strats_B1 <- c()
  strats_B2 <- c()
  strats_B3 <- c()
  opt1_B0 <- c()
  opt1_B1 <- c()
  opt1_B2 <- c()
  opt1_B3 <- c()
  opt2_B0 <- c()
  opt2_B1 <- c()
  opt2_B2 <- c()
  opt2_B3 <- c()
  
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
    
    pop$mu <- expit(B0 + B1 * pop$X1 + B2 * pop$X2 + B3 * pop$X3)
    pop$y <- rbinom(N, 1, pop$mu)
    popmodel <- glm(y ~ X1 + X2 + X3, family=binomial, data=pop)
    
    ####
    #### Add measurement error
    if (!(error %in% c("low", "high"))){
      stop("error must be high or low")
    } else if(error == "low" & scenario != "catX"){
      pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.8
      pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.04, 0) + 0.95
    } else if(error == "high" & scenario != "catX"){
      pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.6
      pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.05, 0) + 0.9
    } else if(error == "low" & scenario == "catX"){
      pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.8
      pop$spec <- ifelse(pop$X1 == 1, 0.04, 0) + 0.95
    } else if(error == "high" & scenario == "catX"){
      pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.6
      pop$spec <- ifelse(pop$X1 == 1, 0.05, 0) + 0.9
    }
    pop$y_obs <- rbinom(N, 1, ifelse(pop$y == 1, pop$sens, 1 - pop$spec))
    
    pop.raw <- pop
    
    # Compute error-prone pop model
    popmodel_obs <- glm(y_obs ~ X1 + X2 + X3, family=binomial, data=pop)
    
    ###
    ### Strategy 1 : SRS
    srs <- sample(pop$id,size = n, replace = FALSE)
    srs_data <- pop[pop$id %in% srs,]
    srsmodel <- glm(y ~ X1 + X2 + X3, family=binomial, data=srs_data)
    srs_B0[i] <- coef(srsmodel)[1]
    srs_B1[i] <- coef(srsmodel)[2]
    srs_B2[i] <- coef(srsmodel)[3]
    srs_B3[i] <- coef(srsmodel)[4]
    
    ####
    #### Strategy 2: Case-control sample
    
    ## A: Compute strata
    pop$strat <- pop$y
    
    cases_ss <- min(n/2, sum(pop$y))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    cc_B0[i] <- coef(ccmodel)[1]
    cc_B1[i] <- coef(ccmodel)[2]
    cc_B2[i] <- coef(ccmodel)[3]
    cc_B3[i] <- coef(ccmodel)[4]
      
    ####
    #### Strategy 3: Stratified sample (A-optimal)
    
    pop <- pop.raw
    
    ## A: Compute influence functions
    inflB0 <- inf_fun_logit(popmodel)[,1]
    pop$inflB0 <- inflB0
    inflB1 <- inf_fun_logit(popmodel)[,2]
    pop$inflB1 <- inflB1
    inflB2 <- inf_fun_logit(popmodel)[,3]
    pop$inflB2 <- inflB2
    inflB3 <- inf_fun_logit(popmodel)[,4]
    pop$inflB3 <- inflB3
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y
    if(scenario != "catX"){
      pop <- optimall::split_strata(pop,  strata = "strata",
                                    split_var = "inflB1",
                                    type = "global quantile",
                                    split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB2",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB3",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
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
    stratmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    strat_B0[i] <- coef(stratmodel)[1]
    strat_B1[i] <- coef(stratmodel)[2]
    strat_B2[i] <- coef(stratmodel)[3]
    strat_B3[i] <- coef(stratmodel)[4]
    
    ####
    ### Strategy 4: Optimal Poisson sampling (uses full data MLE)
    
    pop <- pop.raw
    
    ## A: Compute preliminaries
    Mx.inv<- vcov(popmodel)*N
    X <- model.matrix(popmodel)
    Mx_inv_x <- X %*% t(Mx.inv)
    pop$Mx_inv_x_norm <- sqrt(rowSums(Mx_inv_x^2))
    # or equivalently
    # Ihat <- (t(dm) %*% (dm * popmodel$fitted.values * (1 - popmodel$fitted.values))) /nrow(dm)
    # Mx.in <- solve(Ihat)
    pop$p.hat<- fitted(popmodel)
    #pop$optp_raw <- with(pop, abs(y-p.hat)*l2(t(Mx.inv%*%rbind(1,pop$X1, pop$X2, pop$X3))))
    pop$optp_raw <- with(pop, abs(y-p.hat)*Mx_inv_x_norm)
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
    optmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    optp_B0[i] <- coef(optmodel)[1]
    optp_B1[i] <- coef(optmodel)[2]
    optp_B2[i] <- coef(optmodel)[3]
    optp_B3[i] <- coef(optmodel)[4]
     
    ####
    #### Strategy 5: Case-control sampling using surrogate 
    
    pop <- pop.raw
    
    ## A: Compute strata
    pop$strat <- pop$y_obs
    
    cases_ss <- min(n/2, sum(pop$y_obs))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccsmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    ccs_B0[i] <- coef(ccsmodel)[1]
    ccs_B1[i] <- coef(ccsmodel)[2]
    ccs_B2[i] <- coef(ccsmodel)[3]
    ccs_B3[i] <- coef(ccsmodel)[4] 
    
    ####
    #### Strategy 6: Stratified sampling using pilot
    
    pop <- pop.raw
    
    ## A: Compute influence functions from error-prone model
    inflB0_obs <- inf_fun_logit(popmodel_obs)[,1]
    pop$inflB0_obs <- inflB0_obs
    inflB1_obs <- inf_fun_logit(popmodel_obs)[,2]
    pop$inflB1_obs <- inflB1_obs
    inflB2_obs <- inf_fun_logit(popmodel_obs)[,3]
    pop$inflB2_obs <- inflB2_obs
    inflB3_obs <- inf_fun_logit(popmodel_obs)[,4]
    pop$inflB3_obs <- inflB3_obs
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y_obs
    pop <- optimall::split_strata(pop,  strata = "strata",
                                  split_var = "inflB1_obs",
                                  type = "global quantile",
                                  split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB2_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB3_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
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
    stratmodel1 <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    
    ## Now use results from this to improve optimum allocation
    pilot_infs <- data.frame(inf_fun_logit(stratmodel1))
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
    stratmodel2 <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    
    strats_B0[i] <- coef(stratmodel2)[1]
    strats_B1[i] <- coef(stratmodel2)[2]
    strats_B2[i] <- coef(stratmodel2)[3]
    strats_B3[i] <- coef(stratmodel2)[4]
    
    ####
    #### Strategy 7: Optimall surrogate sampling of Marks-Anglin et al (2025)
    
    pop <- pop.raw
    
    X <- as.matrix(cbind(1, pop[,c("X1", "X2", "X3")])) # intercept + covariates
    
    # generate true y and surrogate s
    y <- pop[,"y"]
    s <- pop[,"y_obs"]
    
    ## full sample MLE (benchmark)
    beta.model.y <- glm(y ~ 0 + X,
                        family=binomial(link="logit"))
    b.MLE = beta.model.y$coef
    # all(popmodel$coefficients == b.MLE) # Matches above
    hat.p.y <- 1/(1 + exp(-b.MLE %*%t(X)))
    w.y <- c(hat.p.y * (1-hat.p.y))
    Mx.y <- t(w.y*X) %*% X / n # Same as my Mx
    
    ## case-control sampling, pilot index 
    stage1.weights <- rep(NA, n)
    stage1.weights[s==1] <- 0.5*(1/length(s[s==1]))
    stage1.weights[s==0] <- 0.5*(1/length(s[s==0]))
    
    ## Sample:
    set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
    
    # table(y[set1], s[set1]) 
    converged <- 0
    # check convergence of logistic regression models
    while(min(table(y[set1], s[set1]))<3 | converged < 2){ # avoid numerical error
      set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
      tryCatch({
        ## pilot:  y ~ X to obtain approx p^(beta) 
        beta.model <- weighted.model.set1(set1, N, y, X, stage1.weights)
        converged <- sum(beta.model$converged)
        cat('.')
      }, error = function(e){
        converged <- 0
      })
    }
    
    ## pilot:  y ~ X to obtain approx p^(beta) 
    # beta.model <- weighted.model.set1(set1,n,y,X,stage1.weights)
    # fit = logistf(y~., data=data.frame(y,X[,-1])[set1,], weights = 1/stage1.weights[set1],family='binomial')
    hat.p <- 1/(1 + exp(-beta.model$coef %*%t(X)))
    w <- c(hat.p * (1-hat.p))
    Mx <- t(w*X) %*% X / N 
    
    ## full s ~ X to obtain p^*(gamma)
    gamma.model <- glm(s ~ 0+X, family=binomial(link="logit"))  
    hat.p.star <- 1/(1 + exp(-gamma.model$coef %*%t(X)))
    w.hat <- c(hat.p.star*(1-hat.p.star))
    Qx <- t(w.hat*X) %*% X / N   
    
    ## use Bayes rule to calculate p(y=1|s,X), instead of above pilot y~S+X
    # p(y=1|s,X) = p(s|y=1,X)p(y=1|X)/{p(s|y=1,X)p(y=1|X)+p(s|y=0,X)p(y=0|X)}
    hat.p.SX <- rep(NA, N)
    # estimated prevalence (p) ppv npv se sp 
    p = mean(s)
    ppv = mean(y[set1][s[set1]==1]==1)
    npv = mean(y[set1][s[set1]==0]==0)
    ppv <- ifelse(ppv == 1, 0.999, ppv)
    npv <- ifelse(npv == 1, 0.999, npv)
    se = (1-p-p*npv/(1-npv)) / (p*((1-ppv)/ppv-npv/(1-npv)))
    sp = (p-(1-p)*ppv/(1-ppv)) / ((1-p)*((1-npv)/npv-ppv/(1-ppv)))
    se = max(min(se, 0.99), 0.5) # truncate within 0.5-1 in extreme cases...
    sp = max(min(sp, 0.99), 0.5) 
    hat.p.SX[s==1] = se*hat.p[s==1] / (se*hat.p[s==1] + (1-sp)*(1-hat.p[s==1]))
    hat.p.SX[s==0] = (1-se)*hat.p[s==0] / ((1-se)*hat.p[s==0] + sp*(1-hat.p[s==0])) 
    
    ## sub-sampling prob
    numer01 <- c(hat.p.SX - 2*hat.p.SX*hat.p + hat.p^2) # OSSAT SSP (Proposition 1)
    ssp.OSSAT <- sqrt(numer01)*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)))
    ssp.OSSAT <- ssp.OSSAT / sum(ssp.OSSAT) 
    
    ## model 
    optmodel1 <- weighted.model.seq4(set1,N,n-r1,replace=TRUE,ssp=ssp.OSSAT,
                                     y,X,stage1.weights)
    
    opt1_B0[i] <- coef(optmodel1)[1]
    opt1_B1[i] <- coef(optmodel1)[2]
    opt1_B2[i] <- coef(optmodel1)[3]
    opt1_B3[i] <- coef(optmodel1)[4]
    
    ####
    ####
    #### Strategy 8 (Not used in paper): Marks-Anglin implementation of Wang et al approach
    hat.p <- 1/(1 + exp(-beta.model.y$coef %*%t(X))) # uses true beta model
    w <- c(hat.p * (1-hat.p))
    ssp.OSMAC <- abs(c(y-hat.p))*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)) )
    ssp.OSMAC <- ssp.OSMAC/sum(ssp.OSMAC)
    #optmodel3 <- weighted.model.seq4(set1,N,n-sum(in_sample),replace=TRUE,ssp=ssp.OSMAC,y,
    #                                 X,stage1.weights)
    #in_sample <- rbinom(N, 1, ssp.OSMAC)
    the_sample <- sample(seq(1,N,1), n, replace=TRUE, prob=ssp.OSMAC)
    optmodel2 <- weighted.model.set1(the_sample,N,y,X,ssp.OSMAC)
    opt2_B0[i] <- coef(optmodel2)[1]
    opt2_B1[i] <- coef(optmodel2)[2]
    opt2_B2[i] <- coef(optmodel2)[3]
    opt2_B3[i] <- coef(optmodel2)[4]
      
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
    method = "Surrogate Optimal Individual prob",
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
  results <- rbind(results, data.frame(
    method = "Optimal Poisson prob Wang M-A calculation",
    mean_B0 = mean(opt2_B0),
    mean_B1 = mean(opt2_B1),
    mean_B2 = mean(opt2_B2),
    mean_B3 = mean(opt2_B3),
    var_B0 = var(opt2_B0),
    var_B1 = var(opt2_B1),
    var_B2 = var(opt2_B2),
    var_B3 = var(opt2_B3),
    MSE_B0 = mean((opt2_B0 - B0)^2),
    MSE_B1 = mean((opt2_B1 - B1)^2),
    MSE_B2 = mean((opt2_B2 - B2)^2),
    MSE_B3 = mean((opt2_B3 - B3)^2),
    MSE_B = mean((opt2_B0 - B0)^2 + (opt2_B1 - B1)^2 + (opt2_B2 - B2)^2 + (opt2_B3 - B3)^2)
  ))
  results$varsum <- results$var_B0 + results$var_B1 + results$var_B2 + results$var_B3
  results$MSEsum <- results$MSE_B0 + results$MSE_B1 + results$MSE_B2 + results$MSE_B3
  return(results)
}

#####
##### Three covs with full data generated only once - only randomness is sampling

ThreeCovsDataOutside <- function(N, n, nreps = 1000, scenario, r1, error = "low"){
  # Setup
  expit<-function(eta) exp(eta)/(1+exp(eta))
  logit<-function(p) log(p/(1-p))
  l2<-function(v) sqrt(sum(v*v))
  
  srs_B0 <- c()
  srs_B1 <- c()
  srs_B2 <- c()
  srs_B3 <- c()
  cc_B0 <- c()
  cc_B1 <- c()
  cc_B2 <- c()
  cc_B3 <- c()
  strat_B0 <- c()
  strat_B1 <- c()
  strat_B2 <- c()
  strat_B3 <- c()
  optp_B0 <- c()
  optp_B1 <- c()
  optp_B2 <- c()
  optp_B3 <- c()
  ccs_B0 <- c()
  ccs_B1 <- c()
  ccs_B2 <- c()
  ccs_B3 <- c()
  strats_B0 <- c()
  strats_B1 <- c()
  strats_B2 <- c()
  strats_B3 <- c()
  opt1_B0 <- c()
  opt1_B1 <- c()
  opt1_B2 <- c()
  opt1_B3 <- c()
  opt2_B0 <- c()
  opt2_B1 <- c()
  opt2_B2 <- c()
  opt2_B3 <- c()
  
  # Setting
  B0 <- ifelse(scenario == "Exp", -0.5, 0.5)
  B1 <- 0.5
  B2 <- 0.5
  B3 <- 0.5
  
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
  
  pop$mu <- expit(B0 + B1 * pop$X1 + B2 * pop$X2 + B3 * pop$X3)
  pop$y <- rbinom(N, 1, pop$mu)
  popmodel <- glm(y ~ X1 + X2 + X3, family=binomial, data=pop)
  
  ####
  #### Add measurement error
  if (!(error %in% c("low", "high"))){
    stop("error must be high or low")
  } else if(error == "low" & scenario != "catX"){
    pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.8
    pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.04, 0) + 0.95
  } else if(error == "high" & scenario != "catX"){
    pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.6
    pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.05, 0) + 0.9
  } else if(error == "low" & scenario == "catX"){
    pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.8
    pop$spec <- ifelse(pop$X1 == 1, 0.04, 0) + 0.95
  } else if(error == "high" & scenario == "catX"){
    pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.6
    pop$spec <- ifelse(pop$X1 == 1, 0.05, 0) + 0.9
  }
  pop$y_obs <- rbinom(N, 1, ifelse(pop$y == 1, pop$sens, 1 - pop$spec))
  
  pop.raw <- pop
  
  # Compute error-prone pop model
  popmodel_obs <- glm(y_obs ~ X1 + X2 + X3, family=binomial, data=pop)
  
  # Loop over nreps (include data generation in loop)
  for(i in 1:nreps){
    # my_i <<- i
    set.seed(i)
    
    ###
    ### Strategy 1 : SRS
    srs <- sample(pop$id,size = n, replace = FALSE)
    srs_data <- pop[pop$id %in% srs,]
    srsmodel <- glm(y ~ X1 + X2 + X3, family=binomial, data=srs_data)
    srs_B0[i] <- coef(srsmodel)[1]
    srs_B1[i] <- coef(srsmodel)[2]
    srs_B2[i] <- coef(srsmodel)[3]
    srs_B3[i] <- coef(srsmodel)[4]
    
    ####
    #### Strategy 2: Case-control sample
    
    ## A: Compute strata
    pop$strat <- pop$y
    
    cases_ss <- min(n/2, sum(pop$y))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    cc_B0[i] <- coef(ccmodel)[1]
    cc_B1[i] <- coef(ccmodel)[2]
    cc_B2[i] <- coef(ccmodel)[3]
    cc_B3[i] <- coef(ccmodel)[4]
    
    ####
    #### Strategy 3: Stratified sample (A-optimal)
    
    pop <- pop.raw
    
    ## A: Compute influence functions
    inflB0 <- inf_fun_logit(popmodel)[,1]
    pop$inflB0 <- inflB0
    inflB1 <- inf_fun_logit(popmodel)[,2]
    pop$inflB1 <- inflB1
    inflB2 <- inf_fun_logit(popmodel)[,3]
    pop$inflB2 <- inflB2
    inflB3 <- inf_fun_logit(popmodel)[,4]
    pop$inflB3 <- inflB3
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y
    if(scenario != "catX"){
      pop <- optimall::split_strata(pop,  strata = "strata",
                                    split_var = "inflB1",
                                    type = "global quantile",
                                    split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB2",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB3",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
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
    
    # When scenario is catX, just do proportional allocation because
    # most strata have sd = 0.
    if(scenario == "catX"){
      design$prop <- design$npop /N
      design$stratum_size <- 2 + round(design$prop * (n-2*nstrata))
      design$stratum_size <- pmin(design$npop, design$stratum_size)
      design$probs <- design$stratum_size/design$npop
    }
    
    sampled_data <- optimall::sample_strata(pop, strata = "strata",
                                            id = "id", design_data = design,
                                            probs = "probs",
                                            n_allocated = "stratum_size")
    strat_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strata, probs = ~sampling_prob, 
                           data = strat_data)
    stratmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    strat_B0[i] <- coef(stratmodel)[1]
    strat_B1[i] <- coef(stratmodel)[2]
    strat_B2[i] <- coef(stratmodel)[3]
    strat_B3[i] <- coef(stratmodel)[4]
    
    ####
    ### Strategy 4: Optimal Poisson sampling (my version)
    
    pop <- pop.raw
    
    ## A: Compute preliminaries
    Mx.inv<- vcov(popmodel)*N
    X <- model.matrix(popmodel)
    Mx_inv_x <- X %*% t(Mx.inv)
    pop$Mx_inv_x_norm <- sqrt(rowSums(Mx_inv_x^2))
    # or equivalently
    # Ihat <- (t(dm) %*% (dm * popmodel$fitted.values * (1 - popmodel$fitted.values))) /nrow(dm)
    # Mx.in <- solve(Ihat)
    pop$p.hat<- fitted(popmodel)
    #pop$optp_raw <- with(pop, abs(y-p.hat)*l2(t(Mx.inv%*%rbind(1,pop$X1, pop$X2, pop$X3))))
    pop$optp_raw <- with(pop, abs(y-p.hat)*Mx_inv_x_norm)
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
    #svydesign <- svydesign(id = ~id, weights = ~wt, 
    #                       data = the_sample)
    optmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    optp_B0[i] <- coef(optmodel)[1]
    optp_B1[i] <- coef(optmodel)[2]
    optp_B2[i] <- coef(optmodel)[3]
    optp_B3[i] <- coef(optmodel)[4]
    
    ####
    #### Strategy 5: Case-control sampling using surrogate 
    
    pop <- pop.raw
    
    ## A: Compute strata
    pop$strat <- pop$y_obs
    
    cases_ss <- min(n/2, sum(pop$y_obs))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccsmodel <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    ccs_B0[i] <- coef(ccsmodel)[1]
    ccs_B1[i] <- coef(ccsmodel)[2]
    ccs_B2[i] <- coef(ccsmodel)[3]
    ccs_B3[i] <- coef(ccsmodel)[4] 
    
    ####
    #### Strategy 6: Stratified sampling using pilot
    
    pop <- pop.raw
    
    ## A: Compute influence functions from error-prone model
    inflB0_obs <- inf_fun_logit(popmodel_obs)[,1]
    pop$inflB0_obs <- inflB0_obs
    inflB1_obs <- inf_fun_logit(popmodel_obs)[,2]
    pop$inflB1_obs <- inflB1_obs
    inflB2_obs <- inf_fun_logit(popmodel_obs)[,3]
    pop$inflB2_obs <- inflB2_obs
    inflB3_obs <- inf_fun_logit(popmodel_obs)[,4]
    pop$inflB3_obs <- inflB3_obs
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y_obs
    pop <- optimall::split_strata(pop,  strata = "strata",
                                  split_var = "inflB1_obs",
                                  type = "global quantile",
                                  split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB2_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB3_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
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
    
    # If scenario is catX just do proportional allocation, because most
    # times sd of IFs will be 0 within strata
    if(scenario == "catX"){
      design$prop <- design$npop /N
      design$stratum_size <- 2 + round(design$prop * (n-2*nstrata))
      design$stratum_size <- pmin(design$npop, design$stratum_size)
      design$probs1 <- design$stratum_size/design$npop
    }
    
    ## Pilot
    sampled_data <- optimall::sample_strata(pop, strata = "strata",
                                            id = "id", design_data = design,
                                            probs = "probs1",
                                            n_allocated = "stratum_size")
    sampled_data$sampled_pilot <- sampled_data$sample_indicator
    strat_data <- sampled_data[sampled_data$sample_indicator == 1,]
    svydesign <- svydesign(id = ~id, strata = ~strata, probs = ~sampling_prob, 
                           data = strat_data)
    stratmodel1 <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    
    ## Now use results from this to improve optimum allocation
    pilot_infs <- data.frame(inf_fun_logit(stratmodel1))
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
    stratmodel2 <- svyglm(y ~ X1 + X2 + X3, design = svydesign, family=quasibinomial)
    
    strats_B0[i] <- coef(stratmodel2)[1]
    strats_B1[i] <- coef(stratmodel2)[2]
    strats_B2[i] <- coef(stratmodel2)[3]
    strats_B3[i] <- coef(stratmodel2)[4]
    
    ####
    #### Strategy 7: Optimall surrogate sampling of Marks-Anglin et al (2025)
    
    pop <- pop.raw
    
    X <- as.matrix(cbind(1, pop[,c("X1", "X2", "X3")])) # intercept + covariates
    
    # generate true y and surrogate s
    y <- pop[,"y"]
    s <- pop[,"y_obs"]
    
    ## full sample MLE (benchmark)
    beta.model.y <- glm(y ~ 0 + X,
                        family=binomial(link="logit"))
    b.MLE = beta.model.y$coef
    # all(popmodel$coefficients == b.MLE) # Matches above
    hat.p.y <- 1/(1 + exp(-b.MLE %*%t(X)))
    w.y <- c(hat.p.y * (1-hat.p.y))
    Mx.y <- t(w.y*X) %*% X / n # Same as my Mx
    
    ## case-control sampling, pilot index 
    stage1.weights <- rep(NA, n)
    stage1.weights[s==1] <- 0.5*(1/length(s[s==1]))
    stage1.weights[s==0] <- 0.5*(1/length(s[s==0]))
    
    ## Sample:
    set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
    
    # table(y[set1], s[set1]) 
    converged <- 0
    # check convergence of logistic regression models
    while(min(table(y[set1], s[set1]))<3 | converged < 2){ # avoid numerical error
      set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
      tryCatch({
        ## pilot:  y ~ X to obtain approx p^(beta) 
        beta.model <- weighted.model.set1(set1, N, y, X, stage1.weights)
        converged <- sum(beta.model$converged)
        cat('.')
      }, error = function(e){
        converged <- 0
      })
    }
    
    ## pilot:  y ~ X to obtain approx p^(beta) 
    # beta.model <- weighted.model.set1(set1,n,y,X,stage1.weights)
    # fit = logistf(y~., data=data.frame(y,X[,-1])[set1,], weights = 1/stage1.weights[set1],family='binomial')
    hat.p <- 1/(1 + exp(-beta.model$coef %*%t(X)))
    w <- c(hat.p * (1-hat.p))
    Mx <- t(w*X) %*% X / N 
    
    ## full s ~ X to obtain p^*(gamma)
    gamma.model <- glm(s ~ 0+X, family=binomial(link="logit"))  
    hat.p.star <- 1/(1 + exp(-gamma.model$coef %*%t(X)))
    w.hat <- c(hat.p.star*(1-hat.p.star))
    Qx <- t(w.hat*X) %*% X / N   
    
    ## use Bayes rule to calculate p(y=1|s,X), instead of above pilot y~S+X
    # p(y=1|s,X) = p(s|y=1,X)p(y=1|X)/{p(s|y=1,X)p(y=1|X)+p(s|y=0,X)p(y=0|X)}
    hat.p.SX <- rep(NA, N)
    # estimated prevalence (p) ppv npv se sp 
    p = mean(s)
    ppv = mean(y[set1][s[set1]==1]==1)
    npv = mean(y[set1][s[set1]==0]==0)
    ppv <- ifelse(ppv == 1, 0.999, ppv)
    npv <- ifelse(npv == 1, 0.999, npv)
    se = (1-p-p*npv/(1-npv)) / (p*((1-ppv)/ppv-npv/(1-npv)))
    sp = (p-(1-p)*ppv/(1-ppv)) / ((1-p)*((1-npv)/npv-ppv/(1-ppv)))
    se = max(min(se, 0.99), 0.5) # truncate within 0.5-1 in extreme cases...
    sp = max(min(sp, 0.99), 0.5) 
    hat.p.SX[s==1] = se*hat.p[s==1] / (se*hat.p[s==1] + (1-sp)*(1-hat.p[s==1]))
    hat.p.SX[s==0] = (1-se)*hat.p[s==0] / ((1-se)*hat.p[s==0] + sp*(1-hat.p[s==0])) 
    
    ## sub-sampling prob
    numer01 <- c(hat.p.SX - 2*hat.p.SX*hat.p + hat.p^2) # OSSAT SSP (Proposition 1)
    ssp.OSSAT <- sqrt(numer01)*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)))
    ssp.OSSAT <- ssp.OSSAT / sum(ssp.OSSAT) 
    
    ## model 
    optmodel1 <- weighted.model.seq4(set1,N,n-r1,replace=TRUE,ssp=ssp.OSSAT,
                                     y,X,stage1.weights)
    
    opt1_B0[i] <- coef(optmodel1)[1]
    opt1_B1[i] <- coef(optmodel1)[2]
    opt1_B2[i] <- coef(optmodel1)[3]
    opt1_B3[i] <- coef(optmodel1)[4]
    
    ####
    ####
    #### Strategy 8 (Not used in paper): Marks-Anglin implementation of Wang et al approach
    hat.p <- 1/(1 + exp(-beta.model.y$coef %*%t(X))) # uses true beta model
    w <- c(hat.p * (1-hat.p))
    ssp.OSMAC <- abs(c(y-hat.p))*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)) )
    ssp.OSMAC <- ssp.OSMAC/sum(ssp.OSMAC)
    #optmodel3 <- weighted.model.seq4(set1,N,n-sum(in_sample),replace=TRUE,ssp=ssp.OSMAC,y,
    #                                 X,stage1.weights)
    #in_sample <- rbinom(N, 1, ssp.OSMAC)
    the_sample <- sample(seq(1,N,1), n, replace=TRUE, prob=ssp.OSMAC)
    optmodel2 <- weighted.model.set1(the_sample,N,y,X,ssp.OSMAC)
    opt2_B0[i] <- coef(optmodel2)[1]
    opt2_B1[i] <- coef(optmodel2)[2]
    opt2_B2[i] <- coef(optmodel2)[3]
    opt2_B3[i] <- coef(optmodel2)[4]
    
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
    method = "Surrogate Optimal Individual prob",
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
  results <- rbind(results, data.frame(
    method = "Optimal Poisson prob Wang M-A calculation",
    mean_B0 = mean(opt2_B0),
    mean_B1 = mean(opt2_B1),
    mean_B2 = mean(opt2_B2),
    mean_B3 = mean(opt2_B3),
    var_B0 = var(opt2_B0),
    var_B1 = var(opt2_B1),
    var_B2 = var(opt2_B2),
    var_B3 = var(opt2_B3),
    MSE_B0 = mean((opt2_B0 - B0)^2),
    MSE_B1 = mean((opt2_B1 - B1)^2),
    MSE_B2 = mean((opt2_B2 - B2)^2),
    MSE_B3 = mean((opt2_B3 - B3)^2),
    MSE_B = mean((opt2_B0 - B0)^2 + (opt2_B1 - B1)^2 + (opt2_B2 - B2)^2 + (opt2_B3 - B3)^2)
  ))
  results$varsum <- results$var_B0 + results$var_B1 + results$var_B2 + results$var_B3
  results$MSEsum <- results$MSE_B0 + results$MSE_B1 + results$MSE_B2 + results$MSE_B3
  return(results)
}

SevenCovs <- function(N, n, nreps = 1000, scenario, r1, error = "low"){
  # Setup
  expit<-function(eta) exp(eta)/(1+exp(eta))
  logit<-function(p) log(p/(1-p))
  l2<-function(v) sqrt(sum(v*v))
  
  srs_B0 <- c()
  srs_B1 <- c()
  srs_B2 <- c()
  srs_B3 <- c()
  srs_B4 <- c()
  srs_B5 <- c()
  srs_B6 <- c()
  srs_B7 <- c()
  cc_B0 <- c()
  cc_B1 <- c()
  cc_B2 <- c()
  cc_B3 <- c()
  cc_B4 <- c()
  cc_B5 <- c()
  cc_B6 <- c()
  cc_B7 <- c()
  strat_B0 <- c()
  strat_B1 <- c()
  strat_B2 <- c()
  strat_B3 <- c()
  strat_B4 <- c()
  strat_B5 <- c()
  strat_B6 <- c()
  strat_B7 <- c()
  optp_B0 <- c()
  optp_B1 <- c()
  optp_B2 <- c()
  optp_B3 <- c()
  optp_B4 <- c()
  optp_B5 <- c()
  optp_B6 <- c()
  optp_B7 <- c()
  ccs_B0 <- c()
  ccs_B1 <- c()
  ccs_B2 <- c()
  ccs_B3 <- c()
  ccs_B4 <- c()
  ccs_B5 <- c()
  ccs_B6 <- c()
  ccs_B7 <- c()
  strats_B0 <- c()
  strats_B1 <- c()
  strats_B2 <- c()
  strats_B3 <- c()
  strats_B4 <- c()
  strats_B5 <- c()
  strats_B6 <- c()
  strats_B7 <- c()
  opt1_B0 <- c()
  opt1_B1 <- c()
  opt1_B2 <- c()
  opt1_B3 <- c()
  opt1_B4 <- c()
  opt1_B5 <- c()
  opt1_B6 <- c()
  opt1_B7 <- c()
  opt2_B0 <- c()
  opt2_B1 <- c()
  opt2_B2 <- c()
  opt2_B3 <- c()
  opt2_B4 <- c()
  opt2_B5 <- c()
  opt2_B6 <- c()
  opt2_B7 <- c()
  
  # Setting
  B0 <- ifelse(scenario == "Exp", -0.5, 0.5)
  B1 <- 0.5
  B2 <- 0.5
  B3 <- 0.5
  B4 <- 0.5
  B5 <- 0.5
  B6 <- 0.5
  B7 <- 0.5
  
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
                                          mu = c(0, 0, 0, 0, 0, 0, 0),
                                          Sigma = matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1), 
                                                         nrow = 7)))
    } else if(scenario == "unequalVar"){
      pop <- data.frame(id = 1:N, 
                        mvrnorm(n = N,
                                mu = c(0, 0, 0, 0, 0, 0, 0),
                                Sigma = matrix(c(1, 0.5/2, 0.5/3, 0.5/4, 0.5/5, 0.5/6, 0.5/7,
                                                 0.5/2, 1/4, 0.5/6, 0.5/8, 0.5/10, 0.5/12, 0.5/14,
                                                 0.5/3, 0.5/6, 1/9, 0.5/12, 0.5/15, 0.5/18, 0.5/21,
                                                 0.5/4, 0.5/8, 0.5/12, 1/16, 0.5/20, 0.5/24, 0.5/28,
                                                 0.5/5, 0.5/10, 0.5/15, 0.5/20, 1/25, 0.5/30, 0.5/35,
                                                 0.5/6, 0.5/12, 0.5/18, 0.5/24, 0.5/30, 1/36, 0.5/42,
                                                 0.5/7, 0.5/14, 0.5/21, 0.5/28, 0.5/35, 0.5/42, 1/49), 
                                               nrow = 7)))
    } else if(scenario == "rareEvent"){
      pop <- data.frame(id = 1:N, mvrnorm(n = N, 
                                          mu = c(-1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6),
                                          Sigma = matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                                                           0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1), 
                                                         nrow = 7)))
    } else if(scenario == "mixNormal"){
      pop <- 0.5*data.frame(id = 1:N, mvrnorm(n = N, 
                                              mu = c(1, 1, 1, 1, 1, 1, 1),
                                              Sigma = matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                               0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                               0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                                                               0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                                                               0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                                                               0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                                                               0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1),
                                                             nrow = 7))) + 
        0.5*data.frame(id = 1:N, mvrnorm(n = N, 
                                         mu = c(-1, -1, -1, -1, -1, -1, -1),
                                         Sigma = matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                          0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                          0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                                                          0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                                                          0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                                                          0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                                                          0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1),
                                                        nrow = 7)))
    } else if(scenario == "T3"){
      pop <- data.frame(mvtnorm::rmvt(n = N, df = 3,
                                      sigma = matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                       0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                                                       0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                                                       0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                                                       0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                                                       0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                                                       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1), nrow = 7)
      ))/10
      pop <- data.frame(id = 1:N, pop)
    } else if(scenario == "Exp"){
      pop <- data.frame(id = 1:N, X1 = rexp(N, rate = 2),
                        X2 = rexp(N, rate = 2),
                        X3 = rexp(N, rate = 2),
                        X4 = rexp(N, rate = 2),
                        X5 = rexp(N, rate = 2),
                        X6 = rexp(N, rate = 2),
                        X7 = rexp(N, rate = 2))
      B0 <- -0.5
    } 
    pop$mu <- expit(B0 + B1 * pop$X1 + B2 * pop$X2 + B3 * pop$X3 + 
                      B4 * pop$X4 + B5 * pop$X5 + B6 * pop$X6 + B7 * pop$X7)
    pop$y <- rbinom(N, 1, pop$mu)
    popmodel <- glm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, family=binomial, data=pop)
    
    ####
    #### Add measurement error
    if (!(error %in% c("low", "high"))){
      stop("error must be high or low")
    } else if(error == "low" & scenario != "catX"){
      pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.8
      pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.04, 0) + 0.95
    } else if(error == "high" & scenario != "catX"){
      pop$sens <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.1, 0) + 0.6
      pop$spec <- ifelse(pop$X1 < quantile(pop$X1, 0.3), 0.05, 0) + 0.9
    } else if(error == "low" & scenario == "catX"){
      pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.8
      pop$spec <- ifelse(pop$X1 == 1, 0.04, 0) + 0.95
    } else if(error == "high" & scenario == "catX"){
      pop$sens <- ifelse(pop$X1 == 0, 0.1, 0) + 0.6
      pop$spec <- ifelse(pop$X1 == 1, 0.05, 0) + 0.9
    }
    pop$y_obs <- rbinom(N, 1, ifelse(pop$y == 1, pop$sens, 1 - pop$spec))
    
    pop.raw <- pop
    
    # Compute error-prone pop model
    popmodel_obs <- glm(y_obs ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, family=binomial, data=pop)
    
    ###
    ### Strategy 1 : SRS
    srs <- sample(pop$id,size = n, replace = FALSE)
    srs_data <- pop[pop$id %in% srs,]
    srsmodel <- glm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, family=binomial, data=srs_data)
    srs_B0[i] <- coef(srsmodel)[1]
    srs_B1[i] <- coef(srsmodel)[2]
    srs_B2[i] <- coef(srsmodel)[3]
    srs_B3[i] <- coef(srsmodel)[4]
    srs_B4[i] <- coef(srsmodel)[5]
    srs_B5[i] <- coef(srsmodel)[6]
    srs_B6[i] <- coef(srsmodel)[7]
    srs_B7[i] <- coef(srsmodel)[8]
    
    ####
    #### Strategy 2: Case-control sample
    
    ## A: Compute strata
    pop$strat <- pop$y
    
    cases_ss <- min(n/2, sum(pop$y))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccmodel <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    cc_B0[i] <- coef(ccmodel)[1]
    cc_B1[i] <- coef(ccmodel)[2]
    cc_B2[i] <- coef(ccmodel)[3]
    cc_B3[i] <- coef(ccmodel)[4]
    cc_B4[i] <- coef(ccmodel)[5]
    cc_B5[i] <- coef(ccmodel)[6]
    cc_B6[i] <- coef(ccmodel)[7]
    cc_B7[i] <- coef(ccmodel)[8]
    
    ####
    #### Strategy 3: Stratified sample (A-optimal)
    
    pop <- pop.raw
    
    ## A: Compute influence functions
    inflB0 <- inf_fun_logit(popmodel)[,1]
    pop$inflB0 <- inflB0
    inflB1 <- inf_fun_logit(popmodel)[,2]
    pop$inflB1 <- inflB1
    inflB2 <- inf_fun_logit(popmodel)[,3]
    pop$inflB2 <- inflB2
    inflB3 <- inf_fun_logit(popmodel)[,4]
    pop$inflB3 <- inflB3
    inflB4 <- inf_fun_logit(popmodel)[,5]
    pop$inflB4 <- inflB4
    inflB5 <- inf_fun_logit(popmodel)[,6]
    pop$inflB5 <- inflB5
    inflB6 <- inf_fun_logit(popmodel)[,7]
    pop$inflB6 <- inflB6
    inflB7 <- inf_fun_logit(popmodel)[,8]
    pop$inflB7 <- inflB7
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y
    if(scenario != "catX"){
      pop <- optimall::split_strata(pop,  strata = "strata",
                                    split_var = "inflB1",
                                    type = "global quantile",
                                    split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB2",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
        dplyr::rename(strata = new_strata) %>%
        optimall::split_strata(strata = "strata",
                               split_var = "inflB3",
                               type = "global quantile",
                               split_at = c(0.2, 0.8)) %>%
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
                                 y = c("inflB0", "inflB1", "inflB2", "inflB3",
                                       "inflB4", "inflB5", "inflB6", "inflB7"),
                                 weights = rep(1/8, times = 8), nsample = n - 2*nstrata, 
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
    stratmodel <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    strat_B0[i] <- coef(stratmodel)[1]
    strat_B1[i] <- coef(stratmodel)[2]
    strat_B2[i] <- coef(stratmodel)[3]
    strat_B3[i] <- coef(stratmodel)[4]
    strat_B4[i] <- coef(stratmodel)[5]
    strat_B5[i] <- coef(stratmodel)[6]
    strat_B6[i] <- coef(stratmodel)[7]
    strat_B7[i] <- coef(stratmodel)[8]
    
    ####
    ### Strategy 4: Optimal Poisson sampling (my version)
    
    pop <- pop.raw
    
    ## A: Compute preliminaries
    Mx.inv<- vcov(popmodel)*N
    X <- model.matrix(popmodel)
    Mx_inv_x <- X %*% t(Mx.inv)
    pop$Mx_inv_x_norm <- sqrt(rowSums(Mx_inv_x^2))
    # or equivalently
    # Ihat <- (t(dm) %*% (dm * popmodel$fitted.values * (1 - popmodel$fitted.values))) /nrow(dm)
    # Mx.in <- solve(Ihat)
    pop$p.hat<- fitted(popmodel)
    #pop$optp_raw <- with(pop, abs(y-p.hat)*l2(t(Mx.inv%*%rbind(1,pop$X1, pop$X2, pop$X3))))
    pop$optp_raw <- with(pop, abs(y-p.hat)*Mx_inv_x_norm)
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
    #svydesign <- svydesign(id = ~id, weights = ~wt, 
    #                       data = the_sample)
    optmodel <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    optp_B0[i] <- coef(optmodel)[1]
    optp_B1[i] <- coef(optmodel)[2]
    optp_B2[i] <- coef(optmodel)[3]
    optp_B3[i] <- coef(optmodel)[4]
    optp_B4[i] <- coef(optmodel)[5]
    optp_B5[i] <- coef(optmodel)[6]
    optp_B6[i] <- coef(optmodel)[7]
    optp_B7[i] <- coef(optmodel)[8]
    
    ####
    #### Strategy 5: Case-control sampling using surrogate 
    
    pop <- pop.raw
    
    ## A: Compute strata
    pop$strat <- pop$y_obs
    
    cases_ss <- min(n/2, sum(pop$y_obs))
    controls_ss <- n - cases_ss
    design <- data.frame(strat = c(0,1), 
                         stratum_size = c(controls_ss, cases_ss), 
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
    ccsmodel <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    ccs_B0[i] <- coef(ccsmodel)[1]
    ccs_B1[i] <- coef(ccsmodel)[2]
    ccs_B2[i] <- coef(ccsmodel)[3]
    ccs_B3[i] <- coef(ccsmodel)[4] 
    ccs_B4[i] <- coef(ccsmodel)[5]
    ccs_B5[i] <- coef(ccsmodel)[6]
    ccs_B6[i] <- coef(ccsmodel)[7]
    ccs_B7[i] <- coef(ccsmodel)[8]
    
    ####
    #### Strategy 6: Stratified sampling using pilot
    
    pop <- pop.raw
    
    ## A: Compute influence functions from error-prone model
    inflB0_obs <- inf_fun_logit(popmodel_obs)[,1]
    pop$inflB0_obs <- inflB0_obs
    inflB1_obs <- inf_fun_logit(popmodel_obs)[,2]
    pop$inflB1_obs <- inflB1_obs
    inflB2_obs <- inf_fun_logit(popmodel_obs)[,3]
    pop$inflB2_obs <- inflB2_obs
    inflB3_obs <- inf_fun_logit(popmodel_obs)[,4]
    pop$inflB3_obs <- inflB3_obs
    inflB4_obs <- inf_fun_logit(popmodel_obs)[,5]
    pop$inflB4_obs <- inflB4_obs
    inflB5_obs <- inf_fun_logit(popmodel_obs)[,6]
    pop$inflB5_obs <- inflB5_obs
    inflB6_obs <- inf_fun_logit(popmodel_obs)[,7]
    pop$inflB6_obs <- inflB6_obs
    inflB7_obs <- inf_fun_logit(popmodel_obs)[,8]
    pop$inflB7_obs <- inflB7_obs
    
    ## B: Compute strata (on quantiles of inflB1, inflB2, inflB3)
    pop$strata <- pop$y_obs
    pop <- optimall::split_strata(pop,  strata = "strata",
                                  split_var = "inflB1_obs",
                                  type = "global quantile",
                                  split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB2_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
      dplyr::rename(strata = new_strata) %>%
      optimall::split_strata(strata = "strata",
                             split_var = "inflB3_obs",
                             type = "global quantile",
                             split_at = c(0.2, 0.8)) %>%
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
                                                               "inflB3_obs",
                                                               "inflB4_obs",
                                                               "inflB5_obs",
                                                               "inflB6_obs",
                                                               "inflB7_obs"),
                                 weights = rep(1/8, times = 8), 
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
    stratmodel1 <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    
    ## Now use results from this to improve optimum allocation
    pilot_infs <- data.frame(inf_fun_logit(stratmodel1))
    names(pilot_infs) <- c("inflB0_pilot", "inflB1_pilot", "inflB2_pilot",
                           "inflB3_pilot", "inflB4_pilot", "inflB5_pilot",
                           "inflB6_pilot", "inflB7_pilot")
    pilot_infs$id <- as.numeric(row.names(pilot_infs))
    pop_new <- dplyr::left_join(sampled_data, pilot_infs, by = "id")
    design2 <- optimall::allocate_wave(pop_new, strata = "strata", y = c("inflB0_pilot", 
                                                                         "inflB1_pilot", 
                                                                         "inflB2_pilot", 
                                                                         "inflB3_pilot",
                                                                         "inflB4_pilot", 
                                                                         "inflB5_pilot", 
                                                                         "inflB6_pilot", 
                                                                         "inflB7_pilot"),
                                       weights = rep(1/8, times = 8), nsample = n - r1, 
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
    stratmodel2 <- svyglm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, design = svydesign, family=quasibinomial)
    
    strats_B0[i] <- coef(stratmodel2)[1]
    strats_B1[i] <- coef(stratmodel2)[2]
    strats_B2[i] <- coef(stratmodel2)[3]
    strats_B3[i] <- coef(stratmodel2)[4]
    strats_B4[i] <- coef(stratmodel2)[5]
    strats_B5[i] <- coef(stratmodel2)[6]
    strats_B6[i] <- coef(stratmodel2)[7]
    strats_B7[i] <- coef(stratmodel2)[8]
    
    ####
    #### Strategy 7: Optimall surrogate sampling of Marks-Anglin et al (2025)
    
    pop <- pop.raw
    
    X <- as.matrix(cbind(1, pop[,c("X1", "X2", "X3", "X4", "X5",
                                   "X6", "X7")])) # intercept + covariates
    
    # generate true y and surrogate s
    y <- pop[,"y"]
    s <- pop[,"y_obs"]
    
    ## full sample MLE (benchmark)
    beta.model.y <- glm(y ~ 0 + X,
                        family=binomial(link="logit"))
    b.MLE = beta.model.y$coef
    # all(popmodel$coefficients == b.MLE) # Matches above
    hat.p.y <- 1/(1 + exp(-b.MLE %*%t(X)))
    w.y <- c(hat.p.y * (1-hat.p.y))
    Mx.y <- t(w.y*X) %*% X / n # Same as my Mx
    
    ## case-control sampling, pilot index 
    stage1.weights <- rep(NA, n)
    stage1.weights[s==1] <- 0.5*(1/length(s[s==1]))
    stage1.weights[s==0] <- 0.5*(1/length(s[s==0]))
    
    ## Sample:
    set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
    
    # table(y[set1], s[set1]) 
    converged <- 0
    # check convergence of logistic regression models
    while(min(table(y[set1], s[set1]))<3 | converged < 2){ # avoid numerical error
      set1 <- sample(seq(1,N,1), size=r1,replace=TRUE, prob=stage1.weights) 
      tryCatch({
        ## pilot:  y ~ X to obtain approx p^(beta) 
        beta.model <- weighted.model.set1(set1, N, y, X, stage1.weights)
        converged <- sum(beta.model$converged)
        cat('.')
      }, error = function(e){
        converged <- 0
      })
    }
    
    ## pilot:  y ~ X to obtain approx p^(beta) 
    # beta.model <- weighted.model.set1(set1,n,y,X,stage1.weights)
    # fit = logistf(y~., data=data.frame(y,X[,-1])[set1,], weights = 1/stage1.weights[set1],family='binomial')
    hat.p <- 1/(1 + exp(-beta.model$coef %*%t(X)))
    w <- c(hat.p * (1-hat.p))
    Mx <- t(w*X) %*% X / N 
    
    ## full s ~ X to obtain p^*(gamma)
    gamma.model <- glm(s ~ 0+X, family=binomial(link="logit"))  
    hat.p.star <- 1/(1 + exp(-gamma.model$coef %*%t(X)))
    w.hat <- c(hat.p.star*(1-hat.p.star))
    Qx <- t(w.hat*X) %*% X / N   
    
    ## use Bayes rule to calculate p(y=1|s,X), instead of above pilot y~S+X
    # p(y=1|s,X) = p(s|y=1,X)p(y=1|X)/{p(s|y=1,X)p(y=1|X)+p(s|y=0,X)p(y=0|X)}
    hat.p.SX <- rep(NA, N)
    # estimated prevalence (p) ppv npv se sp 
    p = mean(s)
    ppv = mean(y[set1][s[set1]==1]==1)
    npv = mean(y[set1][s[set1]==0]==0)
    ppv <- ifelse(ppv == 1, 0.999, ppv)
    npv <- ifelse(npv == 1, 0.999, npv)
    se = (1-p-p*npv/(1-npv)) / (p*((1-ppv)/ppv-npv/(1-npv)))
    sp = (p-(1-p)*ppv/(1-ppv)) / ((1-p)*((1-npv)/npv-ppv/(1-ppv)))
    se = max(min(se, 0.99), 0.5) # truncate within 0.5-1 in extreme cases...
    sp = max(min(sp, 0.99), 0.5) 
    hat.p.SX[s==1] = se*hat.p[s==1] / (se*hat.p[s==1] + (1-sp)*(1-hat.p[s==1]))
    hat.p.SX[s==0] = (1-se)*hat.p[s==0] / ((1-se)*hat.p[s==0] + sp*(1-hat.p[s==0])) 
    
    ## sub-sampling prob
    numer01 <- c(hat.p.SX - 2*hat.p.SX*hat.p + hat.p^2) # OSSAT SSP (Proposition 1)
    ssp.OSSAT <- sqrt(numer01)*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)))
    ssp.OSSAT <- ssp.OSSAT / sum(ssp.OSSAT) 
    
    ## model 
    optmodel1 <- weighted.model.seq4(set1,N,n-r1,replace=TRUE,ssp=ssp.OSSAT,
                                     y,X,stage1.weights)
    
    opt1_B0[i] <- coef(optmodel1)[1]
    opt1_B1[i] <- coef(optmodel1)[2]
    opt1_B2[i] <- coef(optmodel1)[3]
    opt1_B3[i] <- coef(optmodel1)[4]
    opt1_B4[i] <- coef(optmodel1)[5]
    opt1_B5[i] <- coef(optmodel1)[6]
    opt1_B6[i] <- coef(optmodel1)[7]
    opt1_B7[i] <- coef(optmodel1)[8]
    
    ####
    ####
    #### Strategy 8 (Not used in paper): Marks-Anglin implementation of Wang et al approach
    hat.p <- 1/(1 + exp(-beta.model.y$coef %*%t(X))) # uses true beta model
    w <- c(hat.p * (1-hat.p))
    ssp.OSMAC <- abs(c(y-hat.p))*apply(solve(Mx,t(X)), 2, function(a) sqrt(sum(a^2)) )
    ssp.OSMAC <- ssp.OSMAC/sum(ssp.OSMAC)
    #optmodel3 <- weighted.model.seq4(set1,N,n-sum(in_sample),replace=TRUE,ssp=ssp.OSMAC,y,
    #                                 X,stage1.weights)
    #in_sample <- rbinom(N, 1, ssp.OSMAC)
    the_sample <- sample(seq(1,N,1), n, replace=TRUE, prob=ssp.OSMAC)
    optmodel2 <- weighted.model.set1(the_sample,N,y,X,ssp.OSMAC)
    opt2_B0[i] <- coef(optmodel2)[1]
    opt2_B1[i] <- coef(optmodel2)[2]
    opt2_B2[i] <- coef(optmodel2)[3]
    opt2_B3[i] <- coef(optmodel2)[4]
    opt2_B4[i] <- coef(optmodel2)[5]
    opt2_B5[i] <- coef(optmodel2)[6]
    opt2_B6[i] <- coef(optmodel2)[7]
    opt2_B7[i] <- coef(optmodel2)[8]
  }
  
  # Gather results
  results <- data.frame(
    method = "SRS",
    mean_B0 = mean(srs_B0),
    mean_B1 = mean(srs_B1),
    mean_B2 = mean(srs_B2),
    mean_B3 = mean(srs_B3),
    mean_B4 = mean(srs_B4),
    mean_B5 = mean(srs_B5),
    mean_B6 = mean(srs_B6),
    mean_B7 = mean(srs_B7),
    var_B0 = var(srs_B0),
    var_B1 = var(srs_B1),
    var_B2 = var(srs_B2), 
    var_B3 = var(srs_B3),
    var_B4 = var(srs_B4),
    var_B5 = var(srs_B5),
    var_B6 = var(srs_B6),
    var_B7 = var(srs_B7),
    MSE_B0 = mean((srs_B0 - B0)^2),
    MSE_B1 = mean((srs_B1 - B1)^2),
    MSE_B2 = mean((srs_B2 - B2)^2),
    MSE_B3 = mean((srs_B3 - B3)^2),
    MSE_B4 = mean((srs_B4 - B4)^2),
    MSE_B5 = mean((srs_B5 - B5)^2),
    MSE_B6 = mean((srs_B6 - B6)^2),
    MSE_B7 = mean((srs_B7 - B7)^2),
    MSE_B = mean((srs_B0 - B0)^2 + (srs_B1 - B1)^2 + (srs_B2 - B2)^2 +
                     (srs_B3 - B3)^2 + (srs_B4 - B4)^2 + (srs_B5 - B5)^2 +
                     (srs_B6 - B6)^2 + (srs_B7 - B7)^2)
  )
  results <- rbind(results, data.frame(
    method = "Case-control",
    mean_B0 = mean(cc_B0),
    mean_B1 = mean(cc_B1),
    mean_B2 = mean(cc_B2),
    mean_B3 = mean(cc_B3),
    mean_B4 = mean(cc_B4),
    mean_B5 = mean(cc_B5),
    mean_B6 = mean(cc_B6),
    mean_B7 = mean(cc_B7),
    var_B0 = var(cc_B0),
    var_B1 = var(cc_B1),
    var_B2 = var(cc_B2),
    var_B3 = var(cc_B3),
    var_B4 = var(cc_B4),
    var_B5 = var(cc_B5),
    var_B6 = var(cc_B6),
    var_B7 = var(cc_B7),
    MSE_B0 = mean((cc_B0 - B0)^2),
    MSE_B1 = mean((cc_B1 - B1)^2),
    MSE_B2 = mean((cc_B2 - B2)^2),
    MSE_B3 = mean((cc_B3 - B3)^2),
    MSE_B4 = mean((cc_B4 - B4)^2),
    MSE_B5 = mean((cc_B5 - B5)^2),
    MSE_B6 = mean((cc_B6 - B6)^2),
    MSE_B7 = mean((cc_B7 - B7)^2),
    MSE_B = mean((cc_B0 - B0)^2 + (cc_B1 - B1)^2 + (cc_B2 - B2)^2 +
                     (cc_B3 - B3)^2 + (cc_B4 - B4)^2 + (cc_B5 - B5)^2 +
                     (cc_B6 - B6)^2 + (cc_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Stratified",
    mean_B0 = mean(strat_B0),
    mean_B1 = mean(strat_B1),
    mean_B2 = mean(strat_B2),
    mean_B3 = mean(strat_B3),
    mean_B4 = mean(strat_B4),
    mean_B5 = mean(strat_B5),
    mean_B6 = mean(strat_B6),
    mean_B7 = mean(strat_B7),
    var_B0 = var(strat_B0),
    var_B1 = var(strat_B1),
    var_B2 = var(strat_B2),
    var_B3 = var(strat_B3),
    var_B4 = var(strat_B4),
    var_B5 = var(strat_B5),
    var_B6 = var(strat_B6),
    var_B7 = var(strat_B7),
    MSE_B0 = mean((strat_B0 - B0)^2),
    MSE_B1 = mean((strat_B1 - B1)^2),
    MSE_B2 = mean((strat_B2 - B2)^2),
    MSE_B3 = mean((strat_B3 - B3)^2),
    MSE_B4 = mean((strat_B4 - B4)^2),
    MSE_B5 = mean((strat_B5 - B5)^2),
    MSE_B6 = mean((strat_B6 - B6)^2),
    MSE_B7 = mean((strat_B7 - B7)^2),
    MSE_B = mean((strat_B0 - B0)^2 + (strat_B1 - B1)^2 + (strat_B2 - B2)^2 +
                     (strat_B3 - B3)^2 + (strat_B4 - B4)^2 + (strat_B5 - B5)^2 +
                     (strat_B6 - B6)^2 + (strat_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Optimal Poisson prob",
    mean_B0 = mean(optp_B0),
    mean_B1 = mean(optp_B1),
    mean_B2 = mean(optp_B2),
    mean_B3 = mean(optp_B3),
    mean_B4 = mean(optp_B4),
    mean_B5 = mean(optp_B5),
    mean_B6 = mean(optp_B6),
    mean_B7 = mean(optp_B7),
    var_B0 = var(optp_B0),
    var_B1 = var(optp_B1),
    var_B2 = var(optp_B2),
    var_B3 = var(optp_B3),
    var_B4 = var(optp_B4),
    var_B5 = var(optp_B5),
    var_B6 = var(optp_B6),
    var_B7 = var(optp_B7),
    MSE_B0 = mean((optp_B0 - B0)^2),
    MSE_B1 = mean((optp_B1 - B1)^2),
    MSE_B2 = mean((optp_B2 - B2)^2),
    MSE_B3 = mean((optp_B3 - B3)^2),
    MSE_B4 = mean((optp_B4 - B4)^2),
    MSE_B5 = mean((optp_B5 - B5)^2),
    MSE_B6 = mean((optp_B6 - B6)^2),
    MSE_B7 = mean((optp_B7 - B7)^2),
    MSE_B = mean((optp_B0 - B0)^2 + (optp_B1 - B1)^2 + (optp_B2 - B2)^2 +
                     (optp_B3 - B3)^2 + (optp_B4 - B4)^2 + (optp_B5 - B5)^2 +
                     (optp_B6 - B6)^2 + (optp_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Case-control Surrogate",
    mean_B0 = mean(ccs_B0),
    mean_B1 = mean(ccs_B1),
    mean_B2 = mean(ccs_B2),
    mean_B3 = mean(ccs_B3),
    mean_B4 = mean(ccs_B4),
    mean_B5 = mean(ccs_B5),
    mean_B6 = mean(ccs_B6),
    mean_B7 = mean(ccs_B7),
    var_B0 = var(ccs_B0),
    var_B1 = var(ccs_B1),
    var_B2 = var(ccs_B2),
    var_B3 = var(ccs_B3),
    var_B4 = var(ccs_B4),
    var_B5 = var(ccs_B5),
    var_B6 = var(ccs_B6),
    var_B7 = var(ccs_B7),
    MSE_B0 = mean((ccs_B0 - B0)^2),
    MSE_B1 = mean((ccs_B1 - B1)^2),
    MSE_B2 = mean((ccs_B2 - B2)^2),
    MSE_B3 = mean((ccs_B3 - B3)^2),
    MSE_B4 = mean((ccs_B4 - B4)^2),
    MSE_B5 = mean((ccs_B5 - B5)^2),
    MSE_B6 = mean((ccs_B6 - B6)^2),
    MSE_B7 = mean((ccs_B7 - B7)^2),
    MSE_B = mean((ccs_B0 - B0)^2 + (ccs_B1 - B1)^2 + (ccs_B2 - B2)^2 +
                     (ccs_B3 - B3)^2 + (ccs_B4 - B4)^2 + (ccs_B5 - B5)^2 +
                     (ccs_B6 - B6)^2 + (ccs_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Stratified with Pilot",
    mean_B0 = mean(strats_B0),
    mean_B1 = mean(strats_B1),
    mean_B2 = mean(strats_B2),
    mean_B3 = mean(strats_B3),
    mean_B4 = mean(strats_B4),
    mean_B5 = mean(strats_B5),
    mean_B6 = mean(strats_B6),
    mean_B7 = mean(strats_B7),
    var_B0 = var(strats_B0),
    var_B1 = var(strats_B1),
    var_B2 = var(strats_B2),
    var_B3 = var(strats_B3),
    var_B4 = var(strats_B4),
    var_B5 = var(strats_B5),
    var_B6 = var(strats_B6),
    var_B7 = var(strats_B7),
    MSE_B0 = mean((strats_B0 - B0)^2),
    MSE_B1 = mean((strats_B1 - B1)^2),
    MSE_B2 = mean((strats_B2 - B2)^2),
    MSE_B3 = mean((strats_B3 - B3)^2),
    MSE_B4 = mean((strats_B4 - B4)^2),
    MSE_B5 = mean((strats_B5 - B5)^2),
    MSE_B6 = mean((strats_B6 - B6)^2),
    MSE_B7 = mean((strats_B7 - B7)^2),
    MSE_B = mean((strats_B0 - B0)^2 + (strats_B1 - B1)^2 + (strats_B2 - B2)^2 +
                     (strats_B3 - B3)^2 + (strats_B4 - B4)^2 + (strats_B5 - B5)^2 +
                     (strats_B6 - B6)^2 + (strats_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Surrogate Optimal Individual prob",
    mean_B0 = mean(opt1_B0),
    mean_B1 = mean(opt1_B1),
    mean_B2 = mean(opt1_B2),
    mean_B3 = mean(opt1_B3),
    mean_B4 = mean(opt1_B4),
    mean_B5 = mean(opt1_B5),
    mean_B6 = mean(opt1_B6),
    mean_B7 = mean(opt1_B7),
    var_B0 = var(opt1_B0),
    var_B1 = var(opt1_B1),
    var_B2 = var(opt1_B2),
    var_B3 = var(opt1_B3),
    var_B4 = var(opt1_B4),
    var_B5 = var(opt1_B5),
    var_B6 = var(opt1_B6),
    var_B7 = var(opt1_B7),
    MSE_B0 = mean((opt1_B0 - B0)^2),
    MSE_B1 = mean((opt1_B1 - B1)^2),
    MSE_B2 = mean((opt1_B2 - B2)^2),
    MSE_B3 = mean((opt1_B3 - B3)^2),
    MSE_B4 = mean((opt1_B4 - B4)^2),
    MSE_B5 = mean((opt1_B5 - B5)^2),
    MSE_B6 = mean((opt1_B6 - B6)^2),
    MSE_B7 = mean((opt1_B7 - B7)^2),
    MSE_B = mean((opt1_B0 - B0)^2 + (opt1_B1 - B1)^2 + (opt1_B2 - B2)^2 +
                     (opt1_B3 - B3)^2 + (opt1_B4 - B4)^2 + (opt1_B5 - B5)^2 +
                     (opt1_B6 - B6)^2 + (opt1_B7 - B7)^2)
  ))
  results <- rbind(results, data.frame(
    method = "Optimal Poisson prob Wang M-A calculation",
    mean_B0 = mean(opt2_B0),
    mean_B1 = mean(opt2_B1),
    mean_B2 = mean(opt2_B2),
    mean_B3 = mean(opt2_B3),
    mean_B4 = mean(opt2_B4),
    mean_B5 = mean(opt2_B5),
    mean_B6 = mean(opt2_B6),
    mean_B7 = mean(opt2_B7),
    var_B0 = var(opt2_B0),
    var_B1 = var(opt2_B1),
    var_B2 = var(opt2_B2),
    var_B3 = var(opt2_B3),
    var_B4 = var(opt2_B4),
    var_B5 = var(opt2_B5),
    var_B6 = var(opt2_B6),
    var_B7 = var(opt2_B7),
    MSE_B0 = mean((opt2_B0 - B0)^2),
    MSE_B1 = mean((opt2_B1 - B1)^2),
    MSE_B2 = mean((opt2_B2 - B2)^2),
    MSE_B3 = mean((opt2_B3 - B3)^2),
    MSE_B4 = mean((opt2_B4 - B4)^2),
    MSE_B5 = mean((opt2_B5 - B5)^2),
    MSE_B6 = mean((opt2_B6 - B6)^2),
    MSE_B7 = mean((opt2_B7 - B7)^2),
    MSE_B = mean((opt2_B0 - B0)^2 + (opt2_B1 - B1)^2 + (opt2_B2 - B2)^2 +
                     (opt2_B3 - B3)^2 + (opt2_B4 - B4)^2 + (opt2_B5 - B5)^2 +
                     (opt2_B6 - B6)^2 + (opt2_B7 - B7)^2)
  ))
  results$varsum <- results$var_B0 + results$var_B1 + results$var_B2 + results$var_B3 + 
    results$var_B4 + results$var_B5 + results$var_B6 + results$var_B7
  results$MSEsum <- results$MSE_B0 + results$MSE_B1 + results$MSE_B2 + results$MSE_B3 + 
    results$MSE_B4 + results$MSE_B5 + results$MSE_B6 + results$MSE_B7
  return(results)
}

