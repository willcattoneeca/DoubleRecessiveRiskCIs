# Load required libraries
library(dplyr)
library(PropCIs)
library(tidyverse)

# Load the data
data <- read.csv("20240830_forAussieMit_NuVs.csv")

# Select relevant columns
data <- data %>%
  select(gene, variant, RSID,
         Homozygotes.Count, Heterozygotes.Count, Allele.Count, Allele.Frequency,
         hrcAF, gnomadAF, consequences, clinvar, verdict)

# Group by gene and summarize
data <- data %>%
  group_by(gene) %>%
  summarise(gene_variant_count = n(),
            gene_allele_count = sum(Allele.Count),
            genomes = 5690)

# Calculate confidence intervals and risk
data <- data %>%
  rowwise() %>%
  mutate(
    p = gene_allele_count / genomes,
    ci = list(PropCIs::exactci(gene_allele_count, genomes, 0.95)$conf.int),
    lower_p = ci[[1]],
    upper_p = ci[[2]],
    risk = p^2,
    lower_risk = lower_p^2,
    upper_risk = upper_p^2
  ) %>%
  ungroup() %>%
  select(gene, gene_variant_count, gene_allele_count, genomes, lower_p, p, upper_p, lower_risk, risk, upper_risk)

# Save the results to a CSV file
write.csv(data, 'results_trial_pathogenic_tocalculateLR_bygene.csv', row.names = FALSE)

# Calculate the overall probability of being homozygous at least at one locus
P <- 1 - prod(1 - data$risk, na.rm = TRUE)

# Delta method: calculate partial derivatives and variances
data <- data %>%
  rowwise() %>%
  mutate(
    partial = 2 * p * prod(1 - risk[-cur_group_id()], na.rm = TRUE),
    var_p = p * (1 - p) / (genomes - 1)
  ) %>%
  ungroup()

# Variance of P
var_P <- sum(data$partial^2 * data$var_p, na.rm = TRUE)

# Standard deviation and confidence interval for P
sigma_P <- sqrt(var_P)
CI_lower <- P - 1.96 * sigma_P
CI_upper <- P + 1.96 * sigma_P

# Print the results
cat("Probability of being homozygous at least at one locus:", P, "\n")
cat("95% Confidence Interval for P:", CI_lower, "to", CI_upper, "\n")
