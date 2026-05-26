########
#### Supplementary stratification comparison analysis
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
library(ggplot2)
library(ggsci)

#####
## Source functions
#####
source("SimulationFunctions.R")

######
### Preliminaries
######

expit <- function(eta) exp(eta) / (1 + exp(eta))
logit <- function(p) log(p / (1 - p))
l2 <- function(v) sqrt(sum(v * v))
tr <- function(M) sum(diag(M))

inf_fun_logit <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  (dm * resid(fit, type = "response")) %*% solve(Ihat)
}

N <- 10000
scenario_levels <- c("Exp", "mixNormal", "rareEvent", "T3", "unequalVar", "zeroMean")

make_threecovs_with_split_points <- function(base_fun) {
  transformed_fun <- base_fun
  fun_text <- paste(deparse(body(transformed_fun)), collapse = "\n")
  fun_text <- gsub("split_at = c(0.2, 0.8)", "split_at = split_at", fun_text, fixed = TRUE)
  body(transformed_fun) <- parse(text = fun_text)[[1]]
  formals(transformed_fun)$split_at <- c(0.2, 0.8)
  transformed_fun
}

ThreeCovsSuppStrat <- make_threecovs_with_split_points(ThreeCovs)

run_threecovs_suppstrat <- function(n, split_at, strata_label, error = "low", r1 = 600, nreps = 1000) {
  scenarios <- c("zeroMean", "unequalVar", "rareEvent", "mixNormal", "T3", "Exp")
  methods <- c("Stratified", "Stratified with Pilot")

  # Force promises before parallel evaluation so workers get concrete values.
  split_at <- as.numeric(split_at)
  n <- as.integer(n)
  strata_label <- as.character(strata_label)
  error <- as.character(error)
  r1 <- as.integer(r1)
  nreps <- as.integer(nreps)

  detected_cores <- detectCores()
  ncores <- if (is.na(detected_cores) || detected_cores < 2) {
    1
  } else {
    min(length(scenarios), detected_cores - 1)
  }
  run_scenario <- function(scen) {
    data.frame(
      n = n,
      Scenario = scen,
      StrataDefinition = strata_label,
      ThreeCovsSuppStrat(
        N = N,
        n = n,
        nreps = nreps,
        r1 = r1,
        error = error,
        scenario = scen,
        split_at = split_at,
        methods = methods
      )
    )
  }

  cl <- tryCatch(makeCluster(ncores), error = function(e) NULL)
  if (is.null(cl)) {
    results_list <- lapply(scenarios, run_scenario)
  } else {
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(
      cl,
      c("ThreeCovsSuppStrat", "N", "expit", "l2", "inf_fun_logit", "logit", "tr", "methods")
    )
    parallel::clusterEvalQ(cl, {
      library(MASS)
      library(survey)
      library(dplyr)
      library(optimall)
      library(logistf)
    })
    results_list <- parLapply(cl, scenarios, run_scenario)
  }

  bind_rows(results_list)
}

######
### Supplementary Figure 1 data: low error, r1 = 600, ThreeCovs
######

results_dir <- file.path("results", "SuppStrat")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

split_rules <- list(
  "0.2 / 0.8" = c(0.2, 0.8),
  "0.33 / 0.67" = c(0.33, 0.67)
)

n_values <- c(800, 1200, 1600)

set.seed(101)
supp_results <- bind_rows(lapply(n_values, function(n) {
  bind_rows(lapply(names(split_rules), function(rule_name) {
    run_threecovs_suppstrat(
      n = n,
      split_at = split_rules[[rule_name]],
      strata_label = rule_name,
      error = "low",
      r1 = 600,
      nreps = 1000
    )
  }))
}))

saveRDS(supp_results, file.path(results_dir, "resultsThreeCovsSuppStratLowErrorR1600.rds"))

######
### Supplementary Figure S4: stratified approaches only
######

plot_data <- supp_results %>%
  filter(method %in% c("Stratified", "Stratified with Pilot")) %>%
  mutate(
    method = ifelse(method == "Stratified", "Optimal Stratified", "Optimal Stratified w/ pilot"),
    method = factor(method, levels = c("Optimal Stratified", "Optimal Stratified w/ pilot")),
    Scenario = factor(Scenario, levels = scenario_levels),
    StrataDefinition = factor(StrataDefinition, levels = c("0.2 / 0.8", "0.33 / 0.67"))
  )

figure1_scale_data <- bind_rows(
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n800.rds")) %>%
    mutate(n = 800),
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n1200.rds")) %>%
    mutate(n = 1200),
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n1600.rds")) %>%
    mutate(n = 1600)
) %>%
  filter(method != "SRS") %>%
  group_by(Scenario) %>%
  summarise(MSE_B = max(MSE_B, na.rm = TRUE), .groups = "drop") %>%
  tidyr::crossing(n = range(n_values))

p_supp_strat <- ggplot(
  plot_data,
  aes(x = n, y = MSE_B, color = method, shape = method, linetype = StrataDefinition)
) +
  geom_blank(
    data = figure1_scale_data,
    aes(x = n, y = MSE_B),
    inherit.aes = FALSE
  ) +
  geom_blank(
    data = tidyr::crossing(Scenario = unique(plot_data$Scenario), n = range(n_values), MSE_B = 0),
    aes(x = n, y = MSE_B),
    inherit.aes = FALSE
  ) +
  geom_point(size = 2) +
  geom_line(aes(group = interaction(method, StrataDefinition)), linewidth = 0.8) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Strata definition"
  ) +
  scale_color_manual(values = c(
    "Optimal Stratified" = "#0072B2",
    "Optimal Stratified w/ pilot" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "0.2 / 0.8" = "solid",
    "0.33 / 0.67" = "dashed"
  )) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = margin(t = 0, r = 4, b = 0, l = 4),
    plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 5.5),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave(
  file.path("Figures", "FigureSuppStrat1.png"),
  plot = p_supp_strat,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 300
)

######
### Supplementary Figure 2: Figure 1 plus alternative strata definition
######

figure1_results <- bind_rows(
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n800.rds")) %>%
    mutate(n = 800),
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n1200.rds")) %>%
    mutate(n = 1200),
  readRDS(file.path("results", "resultsThreeCovsPilotLowErrorR1600n1600.rds")) %>%
    mutate(n = 1600)
) %>%
  dplyr::select(n, Scenario, method, MSEsum, MSE_B) %>%
  filter(method != "SRS") %>%
  mutate(
    method = case_when(
      method == "Stratified" ~ "Optimal Stratified",
      method == "Optimal Poisson prob" ~ "Optimal Poisson",
      method == "Stratified with Pilot" ~ "Optimal Stratified w/ pilot",
      method == "Surrogate Optimal Poisson prob Wang M-A calculation" ~
        "Surrogate Optimal Individual prob",
      TRUE ~ method
    ),
    type = ifelse(
      method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"),
      "Optimal",
      "Surrogate-assisted"
    ),
    sampling_method = case_when(
      method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
      method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
      method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified (0.2 / 0.8)",
      TRUE ~ "Other"
    )
  )

alt_strata_results <- supp_results %>%
  filter(
    StrataDefinition == "0.33 / 0.67",
    method %in% c("Stratified", "Stratified with Pilot")
  ) %>%
  mutate(
    method = ifelse(method == "Stratified", "Optimal Stratified", "Optimal Stratified w/ pilot"),
    type = ifelse(method == "Optimal Stratified", "Optimal", "Surrogate-assisted"),
    sampling_method = "Stratified (0.33 / 0.67)"
  )

figure1_with_alt_strata <- bind_rows(figure1_results, alt_strata_results) %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = scenario_levels
    ),
    sampling_method = factor(
      sampling_method,
      levels = c(
        "Case-control",
        "Individual",
        "Stratified (0.2 / 0.8)",
        "Stratified (0.33 / 0.67)"
      )
    ),
    type = factor(type, levels = c("Optimal", "Surrogate-assisted"))
  )

p_supp_strat2 <- ggplot(
  figure1_with_alt_strata,
  aes(x = n, y = MSE_B, color = sampling_method, shape = sampling_method, linetype = type)
) +
  geom_point(size = 1.8) +
  geom_line(aes(group = interaction(method, sampling_method)), linewidth = 0.7) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  scale_linetype_manual(values = c(
    "Optimal" = "dashed",
    "Surrogate-assisted" = "solid"
  )) +
  scale_shape_manual(values = c(
    "Case-control" = 16,
    "Individual" = 17,
    "Stratified (0.2 / 0.8)" = 15,
    "Stratified (0.33 / 0.67)" = 1
  )) +
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = margin(t = 0, r = 4, b = 0, l = 4),
    plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 5.5),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

ggsave(
  file.path("Figures", "FigureS4.png"),
  plot = p_supp_strat2,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 300
)
