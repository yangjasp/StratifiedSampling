########
#### Stratified vs. Individualized sampling: Simulation Study and Data Example figures
######## 

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)

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
p1 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/Figure1.png", plot = p1, width = 7.2, height = 5, units = "in", dpi = 300)


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
p2 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/Figure2.png", plot = p2, width = 7.2, height = 5, units = "in", dpi = 300)

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
p3 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = varsum, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Error) +
  labs(
    #title = "Variance for each approach under discrete covariate scenario",
    x = "n (pilot study size = 200)",
    y = "Variance",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

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
p3 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = varsum, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Error) +
  labs(
    title = "Variance for each approach under discrete covariate scenario",
    x = "n (pilot study size = 600)",
    y = "Variance",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/Figure3.png", plot = p3, width = 7.2, height = 5, units = "in", dpi = 300)

####
####
#### Data Example Figure
dataExampleResultsr1125 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1125.csv")
dataExampleResultsr1125$r1 <- 125
dataExampleResultsr1100 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr1100.csv")
dataExampleResultsr1100$r1 <- 100
dataExampleResultsr175 <- read.csv("/Users/jasper/Library/Mobile Documents/com~apple~CloudDocs/Desktop/UW Biostatistics 2023-24/Stats:HIV RA/StratifiedSampling/SimulationStudy/results/data_example_resultsr175.csv")
dataExampleResultsr175$r1 <- 75

resultsDE <- rbind(dataExampleResultsr175, dataExampleResultsr1100, dataExampleResultsr1125)

# Add sum columns to results
results_summary <- resultsDE[,c("r1","method","MSEsum", "MSE_B")]
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
results_summary$method <- ifelse(results_summary$method == "Surrogate Optimal Poisson prob M-A calculation", "Surrogate Optimal Individual prob", results_summary$method)
results_summary$type <- ifelse(results_summary$method %in% c("Case-control", "Optimal Poisson", "Optimal Stratified"), "Optimal", "Surrogate-assisted")
results_summary$color <- case_when(results_summary$method %in% c("Case-control", "Case-control Surrogate") ~ "Case-control",
                                   results_summary$method %in% c("Optimal Poisson", "Surrogate Optimal Individual prob") ~ "Individual",
                                   results_summary$method %in% c("Optimal Stratified", "Optimal Stratified w/ pilot") ~ "Stratified",
                                   TRUE ~ "Other")

# Create plot
p4 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = r1, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "r1 (subsample size = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  scale_x_continuous(breaks = c(75, 100, 125)) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/Figure4.png", plot = p4, width = 7.2, height = 5, units = "in", dpi = 300)

#######
####### Figure S1: Low error, R1 = 200, Three covs

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
pS1 <- ggplot(dplyr::filter(results_summary, !(method %in% c(
  "SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Scenario",
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  #scale_color_manual(
  #  values = c("Case-control" = "#374e55",
  #               "Case-control Surrogate" = "#374e55",
  #               "Optimal Poisson" = "#df8f44",
  #               "Surrogate Optimal Individual prob" = "#df8f44",
  #               "Optimal Stratified" = "#00a1d5",
  #               "Optimal Stratified w/ pilot" = "#00a1d5",
  #               "Other" = "grey50")) +
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) # +
#guides(
#  color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
#)

ggsave("Figures/FigureS1.png", plot = pS1, width = 7.2, height = 5, units = "in")

#####
#####
##### Figure S2: Three covariates + intercept with a pilot sample, replicate Marks-Anglin et al (2018) under high error, r1 = 200


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
pS2 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 200)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/FigureS2.png", plot = pS2, width = 7.2, height = 5, units = "in", dpi = 300)

### Figure S3: high error, r1 = 600, Three covs

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
pS3 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/FigureS3.png", plot = pS3, width = 7.2, height = 5, units = "in", dpi = 300)

### Figure S4: high error, r1 = 600, Seven Covs

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
pS4 <- ggplot(dplyr::filter(results_summary, !(method %in% c("SRS")
)), aes(x = n, y = MSE_B, color = color, shape = color, linetype = type)) +
  geom_point() +
  geom_line(aes(group = method)) +
  facet_wrap(~ Scenario, scales = "free_y") +
  labs(
    #title = "MSE for each approach under different data generating scenarios ",
    x = "n (pilot study size = 600)",
    y = "Mean Squared Error (MSE)",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  scale_linetype_manual(values = c("Optimal" = "dashed", "Surrogate-assisted" = "solid"))+
  scale_color_jama() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16))  # e.g. solid dots + lines
  )

ggsave("Figures/FigureS4.png", plot = pS4, width = 7.2, height = 5, units = "in", dpi = 300)
