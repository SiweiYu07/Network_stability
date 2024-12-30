
library(betapart)  # Load the 'betapart' package for beta diversity calculations.
library(vegan)     # Load the 'vegan' package for ecological analysis and diversity calculations.

# 1. Abundance-based pair-wise dissimilarities (Bray-Curtis)
abund_phyto <- read.table("CGP_abund.txt", header = TRUE, row.names = 1)  
# Read abundance data from a text file where the first row contains OTU names, 
# and the first column contains sample names. Assign OTU names as row names.


abund_total <- beta.multi.abund(abund_phyto, index.family = "bray")  
# Calculate multiple-site abundance-based beta diversity metrics using Bray-Curtis index.
# Outputs include total dissimilarity (β_total), turnover (β_turn), and nestedness (β_nest).

abund_total  # Display the total beta diversity results.

abund_phyto_detail <- beta.pair.abund(abund_phyto, index.family = "bray")  
# Calculate pairwise abundance-based beta diversity metrics using Bray-Curtis index.
# Outputs are pairwise matrices for turnover, nestedness, and total beta diversity.

abund_phyto_detail  # Display pairwise beta diversity results.

write.table(abund_phyto_detail, 'CGP_abund.csv')  
# Save the pairwise beta diversity results to a CSV file.

# 2. Data transformation for incidence-based analysis
library(dplyr)  # Load the 'dplyr' package for data manipulation.

incid_phyto %>% mutate(across(where(is.numeric), ~ +as.logical(.x)))  
# Transform the abundance dataset into a presence-absence dataset using `mutate`.
# Convert numeric values to logical (TRUE/FALSE) and then to binary (1/0).

incid_phyto %>% mutate_if(is.numeric, ~1 * (. != 0))  
# Alternative way to transform numeric abundance values to binary presence-absence (1/0).

write.table(incid_phyto, 'CGP_incid.csv')  
# Save the incidence-based dataset as a CSV file.

# 3. Incidence-based pair-wise dissimilarities (Jaccard)
incid_phyto <- read.table("CGP_incid.txt", header = TRUE, row.names = 1)  
# Read incidence (presence-absence) data from a text file.
# First row contains OTU names, and the first column contains sample names.

incid_total <- beta.multi(incid_phyto, index.family = "jaccard")  
# Calculate multiple-site incidence-based beta diversity metrics using the Jaccard index.
# Outputs include total dissimilarity (β_total), turnover (β_turn), and nestedness (β_nest).

incid_total  # Display the total incidence-based beta diversity results.

incid_phyto <- beta.pair(incid_phyto, index.family = "jac")  
# Calculate pairwise incidence-based beta diversity metrics using the Jaccard index.
# Outputs are pairwise matrices for turnover, nestedness, and total beta diversity.

incid_phyto  # Display the pairwise incidence-based beta diversity results.

write.table(incid_phyto_detail, 'CGP_incid.csv')  
# Save the pairwise beta diversity results for incidence-based data to a CSV file.
