library(dplyr)
library(data.table)
library(TwoSampleMR)

# Path to your GWAS summary files and SNP list
path_smoking_gwas <- "smokergwas.txt"
path_chd_gwas <- "coronarygwas.tsv"
path_snp_list <- "exclusive_smoking_snps.csv"

# Read the GWAS summary statistics for smoking, with progress
smoking_gwas <- fread(path_smoking_gwas, showProgress = TRUE)
# Read the CHD GWAS data, with progress
chd_gwas <- fread(path_chd_gwas, showProgress = TRUE)
# Define your list of SNPs
exclusive_snps <- fread(path_snp_list, showProgress = TRUE)

# Smoking - exposure object
smoking_gwas <- read_exposure_data(
  filename = path_smoking_gwas,
  snp_col = "MarkerName",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF_A1",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "Pval",
  sep = "\t"
)

# CHD - outcome object
chd_gwas <- read_outcome_data(
  filename = path_chd_gwas,
  snp_col = "variant_id",  # Updated to 'variant_id'
  beta_col = "beta",  # Updated to 'beta' (use 'hm_beta' if you want the harmonized beta)
  se_col = "standard_error",  # Updated to 'standard_error' (use 'hm_standard_error' if available and you want harmonized values)
  eaf_col = "effect_allele_frequency",  # Updated to 'effect_allele_frequency' (use 'hm_effect_allele_frequency' for harmonized values)
  effect_allele_col = "effect_allele",  # Updated to 'effect_allele' (use 'hm_effect_allele' for harmonized values)
  other_allele_col = "other_allele",  # Updated to 'other_allele' (use 'hm_other_allele' for harmonized values)
  pval_col = "p-value",  # Updated to 'p-value'
  sep = "\t"  # The delimiter used in your file (update if different)
)

# Filtering out only the smoking SNPs
smoking_gwas_filtered <- smoking_gwas[smoking_gwas$SNP %in% exclusive_snps$riskAllele, ]

# Harmonizing data
harmonized_data <- harmonise_data(smoking_gwas_filtered, chd_gwas)

# IVW
mr_result <- mr_ivw(
  b_exp = harmonized_data$beta.exposure, 
  se_exp = harmonized_data$se.exposure, 
  b_out = harmonized_data$beta.outcome, 
  se_out = harmonized_data$se.outcome
)

# View the results of the MR analysis
print(mr_result)

# IVW, MR EGGER and Weighted Median
res <- mr(harmonized_data, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res)















# Perform MR-Egger regression
mr_egger_result <- mr(harmonized_data, method = "MR-Egger")

# View the results
print(mr_egger_result)



smoking_gwas_filtered <- smoking_gwas %>% filter(MarkerName %in% exclusive_snps$rs12130857)
chd_gwas_filtered <- chd_gwas %>% filter(SNP %in% exclusive_snps$rs12130857)

# Align alleles and flip effects if necessary
smoking_gwas_filtered <- smoking_gwas_filtered %>%
  mutate(
    # Find the corresponding SNP in the CHD dataset
    chd_ref_allele = chd_gwas_filtered$reference_allele[match(MarkerName, chd_gwas_filtered$SNP)],
    
    # Determine if the effect should be flipped
    flipped = ifelse(A1 != chd_ref_allele, TRUE, FALSE),
    
    # Flip the Beta if necessary
    Beta = ifelse(flipped, -Beta, Beta),
    
    # Swap A1 and A2 if flipped
    A1 = ifelse(flipped, A2, A1),
    A2 = ifelse(flipped, chd_ref_allele, A2)
  ) %>%
  # Optionally, remove the temporary column chd_ref_allele if it's no longer needed
  select(-chd_ref_allele)

# Merge datasets
merged_data <- merge(smoking_gwas_filtered, chd_gwas_filtered, by.x = "MarkerName", by.y = "SNP")

# Quality control
pvalue_threshold <- 0.05
merged_data_filtered <- merged_data %>%
  filter(Pval < pvalue_threshold, pvalue < pvalue_threshold)

# The merged_data_filtered dataframe now contains harmonized data ready for MR analysis

# Assuming merged_data_filtered is your final harmonized dataset
# Format the data for TwoSampleMR
# Assuming merged_data_filtered is your harmonized dataset
# Format the exposure data for TwoSampleMR
exposure_data <- merged_data %>%
  select(SNP = MarkerName, 
         id.exposure = MarkerName, # Use SNP as the identifier if no separate id is available
         beta.exposure = Beta, 
         se.exposure = SE, 
         effect_allele.exposure = A1, 
         other_allele.exposure = A2) %>%
  mutate(exposure = "smoking") # Assign a constant value to the exposure column

# Format the outcome data for TwoSampleMR
outcome_data <- merged_data %>%
  select(SNP = MarkerName, 
         id.outcome = MarkerName, # Use SNP as the identifier if no separate id is available
         beta.outcome = log_odds, 
         se.outcome = log_odds_se, 
         effect_allele.outcome = reference_allele, 
         other_allele.outcome = other_allele) %>%
  mutate(outcome = "CHD") # Assign a constant value to the outcome column

# Harmonize datasets
harmonized_data <- harmonise_data(exposure_data, outcome_data)

# Extract necessary components from the harmonized data
b_exposure <- harmonized_data$beta.exposure
se_exposure <- harmonized_data$se.exposure
b_outcome <- harmonized_data$beta.outcome
se_outcome <- harmonized_data$se.outcome

# Perform MR analysis using the IVW method
mr_result <- mr_ivw(b_exposure, se_exposure, b_outcome, se_outcome)

# View results
print(mr_result)

# Perform MR-Egger regression
mr_egger_result <- mr(harmonized_data, method = "Egger")

# View the results
print(mr_egger_result)




























