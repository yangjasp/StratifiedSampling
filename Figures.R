########
#### Stratified vs. Individualized sampling: Simulation Study and Data Example figures
######## 

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)

# Helpers
scenario_levels <- c("Exp", "mixNormal", "rareEvent", "T3", "unequalVar", "zeroMean")

order_scenario <- function(data) {
  if ("Scenario" %in% names(data)) {
    data$Scenario <- factor(data$Scenario, levels = scenario_levels)
  }
  data
}

figure_legend_scales <- list(
  scale_linetype_manual(values = c(
    "Optimal" = "dashed",
    "Surrogate-assisted" = "solid"
  )),
  scale_shape_manual(values = c(
    "Case-control" = 16,
    "Individual" = 17,
    "Stratified" = 15,
    "Other" = 4
  )),
  scale_color_jama(),
  theme_bw(),
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.spacing.y = grid::unit(0, "pt"),
    legend.margin = margin(t = 0, r = 4, b = 0, l = 4),
    plot.margin = margin(t = 5.5, r = 5.5, b = 10, l = 5.5),
    strip.text = element_text(face = "bold")
  ),
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
    shape = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )
)

save_figure <- function(filename, plot, height = 5.8) {
  ggsave(filename, plot = plot, width = 7.2, height = height, units = "in", dpi = 300)
}

###
###
### Figure 1: Low error, r1 = 600, Three covs


resultsThreeCovsKnownLowErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n800.rds")
resultsThreeCovsKnownLowErrorR1600n800$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n800$n <- 800

resultsThreeCovsKnownLowErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n1200.rds")
resultsThreeCovsKnownLowErrorR1600n1200$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n1200$n <- 1200

resultsThreeCovsKnownLowErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n1600.rds")
resultsThreeCovsKnownLowErrorR1600n1600$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n1600$n <- 1600

resultsThreeCovsKnown <- rbind(resultsThreeCovsKnownLowErrorR1600n800,
                               resultsThreeCovsKnownLowErrorR1600n1200,
                               resultsThreeCovsKnownLowErrorR1600n1600)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p1 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/Figure1.png", p1)


###
###
### Figure 2: Low error, r1 = 600, Seven covs


resultsSevenCovsKnownLowErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotLowErrorR1600n800.rds")
resultsSevenCovsKnownLowErrorR1600n800$Error <- "Low"
resultsSevenCovsKnownLowErrorR1600n800$n <- 800

resultsSevenCovsKnownLowErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotLowErrorR1600n1200.rds")
resultsSevenCovsKnownLowErrorR1600n1200$Error <- "Low"
resultsSevenCovsKnownLowErrorR1600n1200$n <- 1200

resultsSevenCovsKnownLowErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotLowErrorR1600n1600.rds")
resultsSevenCovsKnownLowErrorR1600n1600$Error <- "Low"
resultsSevenCovsKnownLowErrorR1600n1600$n <- 1600

resultsSevenCovsKnown <- rbind(resultsSevenCovsKnownLowErrorR1600n800,
                               resultsSevenCovsKnownLowErrorR1600n1200,
                               resultsSevenCovsKnownLowErrorR1600n1600)

#resultsSevenCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsSevenCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob", "Surrogate Optimal Poisson prob my calc", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p2 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
) & color != "Other"), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/Figure2.png", p2)

###
###
### Figure X: Categorical X Low error, r1 = 200, Three Covs


resultsThreeCovsKnowncatXLowError <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsCatXPilotLowErrorR1200.rds")
resultsThreeCovsKnowncatXLowError$Error <- "Low misclassification"
resultsThreeCovsKnowncatXHighError <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsCatXPilotHighErrorR1200.rds")
resultsThreeCovsKnowncatXHighError$Error <- "High misclassification"
resultsThreeCovsKnowncatX <- rbind(resultsThreeCovsKnowncatXLowError,
                                   resultsThreeCovsKnowncatXHighError)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnowncatX[,c("Error","n","method","varsum")]
results_summary$varsum <- round(results_summary$varsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob", "Surrogate Optimal Poisson prob my calc", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p3 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = n, y = varsum, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Error) +
  labs(
    #title = "Variance for each approach under discrete covariate scenario",
    x = "n (pilot study size = 200)",
    y = "Variance",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

# ggsave("Figures/Figure3.png", plot = p3, width = 7.2, height = 5, units = "in", dpi = 300)

###
###
### Figure 3: Categorical X, low error, r1 = 600


resultsThreeCovsKnowncatXLowError <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsCatXPilotLowErrorR1600.rds")
resultsThreeCovsKnowncatXLowError$Error <- "Low misclassification"
resultsThreeCovsKnowncatXHighError <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsCatXPilotHighErrorR1600.rds")
resultsThreeCovsKnowncatXHighError$Error <- "High misclassification"
resultsThreeCovsKnowncatX <- rbind(resultsThreeCovsKnowncatXLowError,
                                   resultsThreeCovsKnowncatXHighError)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnowncatX[,c("Error","n","method","varsum")]
results_summary$varsum <- round(results_summary$varsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob", "Surrogate Optimal Poisson prob my calc", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p3 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = n, y = varsum, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Error) +
  labs(
    title = "Variance for each approach under discrete covariate scenario",
    x = "n (pilot study size = 600)",
    y = "Variance",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales +
  scale_y_continuous(
    breaks = c(0, 0.05, 0.10),
    labels = c("0.00", "0.05", "0.10")
  )

save_figure("Figures/Figure3.png", p3)

####
####
#### Data Example Figure
dataExampleResultsr1150 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1150.csv")
dataExampleResultsr1150$r1 <- 150
dataExampleResultsr1125 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1125.csv")
dataExampleResultsr1125$r1 <- 125
dataExampleResultsr1100 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1100.csv")
dataExampleResultsr1100$r1 <- 100
dataExampleResultsr175 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr175.csv")
dataExampleResultsr175$r1 <- 75
dataExampleResultsr1175 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1175.csv")
dataExampleResultsr1175$r1 <- 175
dataExampleResultsr10 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr10.csv")
dataExampleResultsr10$r1 <- 200

resultsDE <- rbind(dataExampleResultsr10, dataExampleResultsr175,
                   dataExampleResultsr1100, dataExampleResultsr1125, dataExampleResultsr1150, 
                   dataExampleResultsr1175)


# Add sum columns to results
results_summary <- resultsDE[,c("r1","method","MSEsum")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p4 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
) & color != "Other"), aes(x = r1, y = MSEsum, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "pilot study size (n = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales +
  scale_x_continuous(breaks = c(75, 100, 125, 150, 175, 200)) +
  scale_y_continuous(limits = c(0, 1))

save_figure("Figures/Figure4.png", p4)

#######
####### Figure S2: Low error, R1 = 200, Three covs

resultsThreeCovsKnownLowErrorR1200n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1200n800.rds")
resultsThreeCovsKnownLowErrorR1200n800$Error <- "Low"
resultsThreeCovsKnownLowErrorR1200n800$n <- 800

resultsThreeCovsKnownLowErrorR1200n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1200n1200.rds")
resultsThreeCovsKnownLowErrorR1200n1200$Error <- "Low"
resultsThreeCovsKnownLowErrorR1200n1200$n <- 1200

resultsThreeCovsKnownLowErrorR1200n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1200n1600.rds")
resultsThreeCovsKnownLowErrorR1200n1600$Error <- "Low"
resultsThreeCovsKnownLowErrorR1200n1600$n <- 1600

resultsThreeCovsKnown <- rbind(resultsThreeCovsKnownLowErrorR1200n800,
                               resultsThreeCovsKnownLowErrorR1200n1200,
                               resultsThreeCovsKnownLowErrorR1200n1600)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")


# Create plot
pS2 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c(
  "SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/FigureS2.png", pS2)

#####
#####
##### Figure S3: Three covariates + intercept with a pilot sample, replicate Marks-Anglin et al (2018) under high error, r1 = 200


resultsThreeCovsKnownHighErrorR1200n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1200n800.rds")
resultsThreeCovsKnownHighErrorR1200n800$Error <- "High"
resultsThreeCovsKnownHighErrorR1200n800$n <- 800

resultsThreeCovsKnownHighErrorR1200n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1200n1200.rds")
resultsThreeCovsKnownHighErrorR1200n1200$Error <- "High"
resultsThreeCovsKnownHighErrorR1200n1200$n <- 1200

resultsThreeCovsKnownHighErrorR1200n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1200n1600.rds")
resultsThreeCovsKnownHighErrorR1200n1600$Error <- "High"
resultsThreeCovsKnownHighErrorR1200n1600$n <- 1600

resultsThreeCovsKnown <- rbind(resultsThreeCovsKnownHighErrorR1200n800,
                               resultsThreeCovsKnownHighErrorR1200n1200,
                               resultsThreeCovsKnownHighErrorR1200n1600)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
pS3 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/FigureS3.png", pS3)

### Figure S1: high error, r1 = 600, Three covs

resultsThreeCovsKnownHighErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n800.rds")
resultsThreeCovsKnownHighErrorR1600n800$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n800$n <- 800

resultsThreeCovsKnownHighErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n1200.rds")
resultsThreeCovsKnownHighErrorR1600n1200$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n1200$n <- 1200

resultsThreeCovsKnownHighErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n1600.rds")
resultsThreeCovsKnownHighErrorR1600n1600$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n1600$n <- 1600

resultsThreeCovsKnown <- rbind(resultsThreeCovsKnownHighErrorR1600n800,
                               resultsThreeCovsKnownHighErrorR1600n1200,
                               resultsThreeCovsKnownHighErrorR1600n1600)

#resultsThreeCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsThreeCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
pS1 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/FigureS1.png", pS1)

### Figure S7: high error, r1 = 600, Seven Covs

resultsSevenCovsKnownHighErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotHighErrorR1600n800.rds")
resultsSevenCovsKnownHighErrorR1600n800$Error <- "High"
resultsSevenCovsKnownHighErrorR1600n800$n <- 800

resultsSevenCovsKnownHighErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotHighErrorR1600n1200.rds")
resultsSevenCovsKnownHighErrorR1600n1200$Error <- "High"
resultsSevenCovsKnownHighErrorR1600n1200$n <- 1200

resultsSevenCovsKnownHighErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsSevenCovsPilotHighErrorR1600n1600.rds")
resultsSevenCovsKnownHighErrorR1600n1600$Error <- "High"
resultsSevenCovsKnownHighErrorR1600n1600$n <- 1600

resultsSevenCovsKnown <- rbind(resultsSevenCovsKnownHighErrorR1600n800,
                               resultsSevenCovsKnownHighErrorR1600n1200,
                               resultsSevenCovsKnownHighErrorR1600n1600)

#resultsSevenCovsKnown$method[seq(7,126, by = 7)] <- "Surrogate Optimal Poisson Prob"

# Add sum columns to results
results_summary <- resultsSevenCovsKnown[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

# Pivot wider
#results_summary <- results_summary %>%
#  tidyr::pivot_wider(
#    names_from = c(method),
#    values_from = MSEsum
# ) 
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
pS7 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
) & color != "Other"), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  figure_legend_scales

save_figure("Figures/FigureS7.png", pS7)

### Figure S5: Coverage of 95% confidence intervals, Three covs

resultsThreeCovsKnownLowErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n800.rds")
resultsThreeCovsKnownLowErrorR1600n800$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n800$r1 <- 600
resultsThreeCovsKnownLowErrorR1600n800$n <- 800

resultsThreeCovsKnownLowErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n1200.rds")
resultsThreeCovsKnownLowErrorR1600n1200$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n1200$r1 <- 600
resultsThreeCovsKnownLowErrorR1600n1200$n <- 1200

resultsThreeCovsKnownLowErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotLowErrorR1600n1600.rds")
resultsThreeCovsKnownLowErrorR1600n1600$Error <- "Low"
resultsThreeCovsKnownLowErrorR1600n1600$r1 <- 600
resultsThreeCovsKnownLowErrorR1600n1600$n <- 1600

resultsThreeCovsKnownHighErrorR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n800.rds")
resultsThreeCovsKnownHighErrorR1600n800$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n800$r1 <- 600
resultsThreeCovsKnownHighErrorR1600n800$n <- 800

resultsThreeCovsKnownHighErrorR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n1200.rds")
resultsThreeCovsKnownHighErrorR1600n1200$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n1200$r1 <- 600
resultsThreeCovsKnownHighErrorR1600n1200$n <- 1200

resultsThreeCovsKnownHighErrorR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsPilotHighErrorR1600n1600.rds")
resultsThreeCovsKnownHighErrorR1600n1600$Error <- "High"
resultsThreeCovsKnownHighErrorR1600n1600$r1 <- 600
resultsThreeCovsKnownHighErrorR1600n1600$n <- 1600

resultsThreeCovsKnown <- rbind(resultsThreeCovsKnownLowErrorR1600n800,
                               resultsThreeCovsKnownLowErrorR1600n1200,
                               resultsThreeCovsKnownLowErrorR1600n1600,
                               resultsThreeCovsKnownHighErrorR1600n800,
                               resultsThreeCovsKnownHighErrorR1600n1200,
                               resultsThreeCovsKnownHighErrorR1600n1600)

results_summary <- resultsThreeCovsKnown[,c("Error", "n", "Scenario", "method",
                                            "coverage_B1", "coverage_B2",
                                            "coverage_B3")]
results_summary <- results_summary %>%
  tidyr::pivot_longer(cols = c("coverage_B1", "coverage_B2", "coverage_B3"),
                      names_to = "coef", values_to = "coverage")
results_summary$coef <- gsub("coverage_", "", results_summary$coef)
results_summary <- dplyr::filter(results_summary, coef != "B0")
results_summary$method <- ifelse(results_summary$method == "Stratified", "Optimal Stratified", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Optimal Poisson prob", "Optimal Poisson", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Stratified with Pilot", "Optimal Stratified w/ pilot", results_summary$method)
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob Wang M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")
results_summary$shape_key <- paste0("n=", results_summary$n, ", ", results_summary$type)
results_summary$shape_key <- factor(results_summary$shape_key, 
                                    levels = c("n=800, Optimal", "n=800, Surrogate-assisted",
                                               "n=1200, Optimal", "n=1200, Surrogate-assisted",
                                               "n=1600, Optimal", "n=1600, Surrogate-assisted"))

pS5 <- ggplot(dplyr::filter(order_scenario(results_summary), !(method %in% c("SRS")
)), aes(x = Scenario, y = coverage, color = color, shape = shape_key, fill = type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.45)) +
  facet_grid(Error ~ coef, labeller = labeller(
    Error = c("Low" = "Low misclassification", "High" = "High misclassification"),
    coef = as_labeller(c("B1" = "beta[1]", "B2" = "beta[2]", "B3" = "beta[3]"), label_parsed)
  )) +
  labs(
    x = "Scenario (pilot study size = 600)",
    y = "Coverage of 95% confidence interval",
    color = "Sampling Method",
    shape = "n",
    fill = "Method Type"
  ) +
  scale_shape_manual(values = c(
    "n=800, Optimal" = 17,
    "n=800, Surrogate-assisted" = 2,
    "n=1200, Optimal" = 16,
    "n=1200, Surrogate-assisted" = 1,
    "n=1600, Optimal" = 15,
    "n=1600, Surrogate-assisted" = 0
  )) +
  scale_fill_manual(values = c("Optimal" = "black", "Surrogate-assisted" = NA)) +
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.spacing = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  guides(
    color = guide_legend(
      order = 1,
      ncol = 3,
      keyheight = unit(0.2, "cm"),
      override.aes = list(shape = 18)
    ),
    fill = "none",
    shape = guide_legend(
      order = 2,
      ncol = 2,
      keyheight = unit(0.2, "cm"),
      breaks = c("n=800, Optimal", "n=800, Surrogate-assisted", 
                 "n=1200, Optimal", "n=1200, Surrogate-assisted",
                 "n=1600, Optimal", "n=1600, Surrogate-assisted"),
      labels = c("n=800", "", "n=1200", "", "n=1600", "")
    )
  )

save_figure("Figures/FigureS5.png", pS5, height = 7.2)

### Figure S6: Linear model (Three covariates), pilot-assisted

# Read results and keep the r1 = 600 family used in Figure 1
resultsThreeCovsLinearR1600n800 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsLinearPilotR1600n800.rds")
resultsThreeCovsLinearR1600n800$n <- 800

resultsThreeCovsLinearR1600n1200 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsLinearPilotR1600n1200.rds")
resultsThreeCovsLinearR1600n1200$n <- 1200

resultsThreeCovsLinearR1600n1600 <- readRDS("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/resultsThreeCovsLinearPilotR1600n1600.rds")
resultsThreeCovsLinearR1600n1600$n <- 1600

resultsThreeCovsLinear <- rbind(resultsThreeCovsLinearR1600n800,
                                resultsThreeCovsLinearR1600n1200,
                                resultsThreeCovsLinearR1600n1600)

results_summary <- resultsThreeCovsLinear[,c("n","Scenario","method","MSEsum", "MSE_B")]
results_summary$MSEsum <- round(results_summary$MSEsum, 4)

results_summary$method_plot <- dplyr::case_when(
  results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Outcome-Dependent",
  results_summary$method %in% c("Optimal Poisson prob", "Surrogate Optimal Individual prob", "Error-prone Optimal Individual prob") ~ "Individual",
  results_summary$method %in% c("Stratified", "Optimal Stratified", "Optimal Stratified w/ pilot", "Stratified with Pilot") ~ "Stratified",
  results_summary$method == "SRS" ~ "SRS"
)
results_summary$method_plot <- factor(results_summary$method_plot,
                                     levels = c("Outcome-Dependent", "Individual", "Stratified", "SRS"))
results_summary$type <- dplyr::case_when(
  results_summary$method == "SRS" ~ "SRS",
  results_summary$method %in% c("Case-control", "Optimal Poisson prob", "Stratified", "Optimal Stratified w/ pilot") ~ "Optimal",
  TRUE ~ "Surrogate-assisted"
)
results_summary$type <- factor(results_summary$type, levels = c("Optimal", "Surrogate-assisted", "SRS"))
results_summary <- dplyr::filter(results_summary, method_plot != "Outcome-Dependent")

pS6 <- ggplot(dplyr::filter(order_scenario(results_summary)), aes(x = n, y = MSE_B, color = method_plot, shape = method_plot, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Sampling method",
    shape = "Sampling method",
    linetype = "Method type"
  ) +
  scale_color_manual(values = c(
    "SRS" = "#374e55",
    "Individual" = "#df8f44",
    "Stratified" = "#00a1d5"
  )) +
  scale_shape_manual(values = c(
    "SRS" = 16,
    "Individual" = 17,
    "Stratified" = 15
  )) +
  scale_linetype_manual(values = c(
    "Optimal" = "dashed",
    "Surrogate-assisted" = "solid",
    "SRS" = "dotted"
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

save_figure("Figures/FigureS6.png", pS6)
