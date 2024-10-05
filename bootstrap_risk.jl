# Import necessary packages
using CSV
using DataFrames
using Random
using Statistics
using Distributions
using ProgressMeter

# Set the seed for reproducibility
Random.seed!(12345)

# Read the CSV file into a DataFrame
data = CSV.read("results_trial_pathogenic_tocalculateLR_bygene.csv", DataFrame)

# Calculate the total number of gene copies examined (n_i) for each gene
data.n_i = data.genomes

# For each gene, s_i is the gene_allele_count (number of non-wild-type alleles observed)
data.s_i = data.gene_allele_count

# List of genes
genes = data.gene

# Number of bootstrap iterations
num_bootstrap = 10000

# Initialize an array to store bootstrap estimates of P
P_bootstrap = zeros(num_bootstrap)

# Progress meter to track bootstrap iterations
@showprogress 1 "Performing bootstrap..." for b in 1:num_bootstrap
    # Initialize an array to store p_i^(b) for each gene in this bootstrap iteration
    p_b = zeros(length(genes))
    
    # For each gene, perform bootstrap resampling
    for (index, row) in enumerate(eachrow(data))
        n_i = row.n_i
        s_i = row.s_i
        
        # Create the original dataset for this gene
        # 1 represents a non-wild-type allele (success), 0 represents a wild-type allele (failure)
        original_data = vcat(ones(Int64, s_i), zeros(Int64, n_i - s_i))
        
        # Bootstrap resampling with replacement
        bootstrap_sample = sample(original_data, n_i, replace = true)
        
        # Calculate the bootstrap proportion p_i^(b)
        p_i_b = mean(bootstrap_sample)
        
        # Store p_i_b in the array
        p_b[index] = p_i_b
    end
    
    # Calculate P^(b) for this bootstrap iteration
    # P^(b) = 1 - product over i of (1 - (p_i^(b))^2)
    P_b = 1.0 - prod(1 .- (p_b .^ 2))
    
    # Store P_b in the bootstrap estimates array
    P_bootstrap[b] = P_b
end

# Calculate the 95% confidence interval for P
lower_P = quantile(P_bootstrap, 0.025)
upper_P = quantile(P_bootstrap, 0.975)
mean_P = mean(P_bootstrap)

# Print the results
println("Estimated Risk (P): $mean_P")
println("95% Confidence Interval: [$lower_P, $upper_P]")
