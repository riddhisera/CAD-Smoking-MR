#####################################################
#List of potentially important SNPs
install.packages("readr")
library(readr)

smoke_snp_list <- read_tsv("smokingsnps.tsv")
cad_snp_list <- read_tsv("cadsnps.tsv")

# Selecting only those with a love p-value
smoke_snp_list <- smoke_snp_list[smoke_snp_list$pValue < 10^-7, ]

# Strip allele information from the SNP identifiers in the snp_list
smoke_snp_list$riskAlleleless <- gsub("-.*$", "", smoke_snp_list$riskAllele)
cad_snp_list$riskAlleleless <- gsub("-.*$", "", cad_snp_list$riskAllele)

# Find common SNPs
common_snps <- smoke_snp_list$riskAlleleless %in% cad_snp_list$riskAlleleless 
smoke_snp_list[common_snps, ]

# Remove common SNPs from smoke_snp_list
smoke_snp_list <- smoke_snp_list[!common_snps, ]

# Some of these SNPs come from combined studies and might possess 
# confounding factors such as BMI as indicated by the traitName column
# On visual observation, row number until 230 has only smoking behaviour
# as its traitName and hence we will select only those SNPs

# Selecting the first 230 rows
smoke_snp_list <- smoke_snp_list[1:230, ]

# Exporting the data
write.csv(smoke_snp_list$riskAlleleless, "CADfreeSNPs.csv")

#####################################################

library(R.utils)
aa_100022 <- gunzip("/Users/riddhisera/gwas.tsv.bgz", "/Users/riddhisera/done.tsv")

library(Rsamtools)

# Load the TSV file
data <- read.table("/Users/riddhisera/done.tsv", header = TRUE, sep = "\t")
head(data)

# Sort by chromosome and start position (assuming these are in the first two columns)
# Add two new columns for chromosome and position
data$chromosome <- sapply(strsplit(as.character(data$variant), ":"), function(x) x[1])
data$position <- as.numeric(sapply(strsplit(as.character(data$variant), ":"), function(x) x[2]))
# Sort by chromosome and position
data <- data[order(data$chromosome, data$position),]
# Remove the temporary columns if they are no longer needed
data$chromosome <- NULL
data$position <- NULL

write.table(data, "sorted.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Compress the file
bgzip("sorted.tsv")

# Index the file
tabix("done.tsv.gz", preset="vcf")
tbx <- TabixFile("done.tsv.gz")


# Assuming 'data' is your data frame
library(dplyr)
library(tidyr)

# Splitting the variant column into multiple columns
data <- data %>% 
  separate(variant, into = c("chromosome", "position", "ref", "alt"), sep = ":") %>%
  mutate(position = as.numeric(position))

# Checking the modified data
head(data)

# Write the modified data frame to a TSV file
write.table(data, "modified_data.tsv", sep="\t", quote = FALSE, row.names = FALSE)




