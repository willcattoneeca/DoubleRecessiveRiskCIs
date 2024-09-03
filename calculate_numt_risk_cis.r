# Load required libraries
library(dplyr)
library(PropCIs)
library(tidyverse)

# Load the data
trial_pathogenic_tocalculateLR <- read.csv("20240830_forAussieMit_NuVs.csv")

# Select relevant columns
trial_pathogenic_tocalculateLR <- trial_pathogenic_tocalculateLR %>%
  select(gene, variant, RSID,
         Homozygotes.Count, Heterozygotes.Count, Allele.Count, Allele.Frequency,
         hrcAF, gnomadAF, consequences, clinvar, verdict)

# Group by gene and summarize
trial_pathogenic_tocalculateLR_bygene <- trial_pathogenic_tocalculateLR %>%
  group_by(gene) %>%
  summarise(gene_variant_count = n(),
            gene_allele_count = sum(Allele.Count)) %>%
  mutate(genomes = 5690)

# Prepare the data for confidence interval calculation
data <- trial_pathogenic_tocalculateLR_bygene
data$lower <- NA
data$upper <- NA

# Loop through each row and calculate confidence intervals for p and p^2
for (row in 1:nrow(data)) {
  x = as.numeric(data[row, "gene_allele_count"])
  n = as.numeric(data[row, "genomes"])
  # Check for non-zero values to avoid issues with the exactci function
  if (n > 0 & x >= 0) {
    c = PropCIs::exactci(x, n, 0.95)
    p_hat = x / n
    lower_p = c$conf.int[1]
    upper_p = c$conf.int[2]
    # Store p and the confidence intervals
    data[row, "p"] = p_hat
    data[row, "lower_p"] = lower_p
    data[row, "upper_p"] = upper_p
    # Calculate risk as p^2
    data[row, "risk"] = p_hat^2
    data[row, "lower_risk"] = lower_p^2
    data[row, "upper_risk"] = upper_p^2
  } else {
    data[row, "p"] = NA
    data[row, "lower"] = NA
    data[row, "upper"] = NA
    data[row, "risk"] = NA
    data[row, "lower_risk"] = NA
    data[row, "upper_risk"] = NA
  }
}

# Reorder and rename the columns
data <- data %>%
  select(gene, gene_variant_count, gene_allele_count, genomes, lower_p, p, upper_p,
         lower_risk, risk, upper_risk)

# Save the results to a CSV file
write.csv(data, 'results_trial_pathogenic_tocalculateLR_bygene.csv')

# Calculate the overall probability of being homozygous at least at one locus
P = 1 - prod(1 - data$risk, na.rm = TRUE)

# Delta method: calculate partial derivatives and variances
# Loop through each gene to exclude the i-th term from the product
data$partial <- NA
for (i in 1:nrow(data)) {
  # Exclude the i-th term from the product
  other_risks = data$risk[-i]
  # Calculate the partial derivative
  data$partial[i] = 2 * data$p[i] * prod(1 - other_risks, na.rm = TRUE)
}

data$var_p = data$p * (1 - data$p) / (data$genomes - 1)

# Variance of P
var_P = sum(data$partial^2 * data$var_p, na.rm = TRUE)

# Standard deviation and confidence interval for P
sigma_P = sqrt(var_P)
CI_lower = P - 1.96 * sigma_P
CI_upper = P + 1.96 * sigma_P

# Print the results
cat("Probability of being homozygous at least at one locus:", P, "\n")
cat("95% Confidence Interval for P:", CI_lower, "to", CI_upper, "\n")
