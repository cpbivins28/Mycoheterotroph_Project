---
title: "Mycoheterotroph Project NEATER"
output: html_notebook
---

# 1. Install and load packages
## Packages Installation
```{r}
# Packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")


install.packages("phyloseq")
install.packages("file2meco")
install.packages("microeco")
install.packages("ggplot2")
install.packages("vegan")
install.packages("tidyr")
install.packages("dplyr")
install.packages("MASS")
install.packages("cowplot")
install.packages("data.table")
install.packages("tidyverse")
install.packages("boot")
install.packages("dunn.test")
install.packages("agricolae")
install.packages("lawstat")
install.packages("patchwork")
install.packages("lmerTest")
install.packages("lme4")
```

### Load packages:
```{r}
library(biomformat)
library(phyloseq)
library(file2meco)
library(microeco)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(MASS)
library(cowplot)
library(data.table)
library(tidyverse)
library(boot)
library(forcats)
library(purrr)
library(broom)
library(dunn.test) 
library(stringr)
library(stats)
library(agricolae) 
library(lawstat)
library(patchwork)
```




## Project Description:
#### This project investigates the symbiotic relationship between two mycoheterotrophic plants (*Sarcodes* *sanguinea* and *Corallorhiza* *striata*) and their fungal hosts (*Rhizopogon* *ellenae* and *Tomentella* *fuscocinerea*, respectively). Specifically, this project utilizes fungal DNA metabarcoding to investigate the relative abundance of these fungi and how this abundance changes based on proximity to their respective mycoheterotroph associates. Bidartondo et al. (2000) found that the absolute abundance of *R. ellenae*, as measured on ectomycorrhizal (ECM) fir root tips is extremely high within the rootball of *S. sanguinea*, and immediately decreases at further distances away from the plant. Additionally, other ECM fungi are far less abundant or entirely absent from fir root tips within the rootball of *S. sanguinea*, but become more abundant and heterogenous at further distances away from the plant. This study aims to replicate Bidartondo et al.'s (2000) results with a fungal DNA metabarcoding approach from rhizosphere soil adhering to fir roots. Additionally, this study aims to replicate Bidartondo et al.'s (2000) study on another mycoheterotroph, *C. striata* to see if its host fungus (*T. fuscocinerea*) also follows a similar pattern. For this aim, we sequenced fungal DNA from rhizosphere soil adhering to fir roots with the same sampling design as Bidartondo et al. (2000), in addition to collecting root tips at each distance class and pooling them, and collectively sequencing them with fungal DNA metabarcoding. Using a phylogenetic approach, we obtained (from GenBank) the sequences of each respective fungus from the papers that identified them as the host fungi of these two mycoheterotrophic plants. We built Maximum Likelihood trees with all of the OTUs identified to Rhizopogonaceae and Thelephoraceaecombined them with the GenBank sequences of *R. ellenae* and *T. fuscocinerea*, respectively. This led us to identify OTU511 as *R. ellenae* and OTU8 as *T. fuscocinerea*.

### 2. Import data, clean up, subset, and assign functional guilds

#### Import biom file with OTU and taxonomy data from amptk and filter out negative control reads. 
```{r}

# Read the .biom file and convert to a phyloseq object
physeq <- import_biom("/Users/cbivins/Desktop/AMPTK_Outputs/AMPTK_Outputs_Mycoheterotroph_Project/taxonomyITS.biom")
# Convert the OTU table to a data.frame
otu_df <- as.data.frame(otu_table(physeq))
# Identify the negative control samples
neg_samples <- grep("NEG", colnames(otu_df), value = TRUE)
# Compute the total count of each OTU across the negative control samples
neg_counts <- rowSums(otu_df[neg_samples])
# Iterate over each OTU found in the negative controls
for (otu in names(neg_counts)) {
  # Check if the OTU is in the sample data
  if (otu %in% rownames(otu_df)) {
    # Subtract the negative control counts from the sample counts, ensuring no negative values
    otu_df[otu,] <- pmax(otu_df[otu,] - neg_counts[otu], 0)
  }
}
# Convert the adjusted DataFrame back to a phyloseq OTU table
otu_table_adj <- otu_table(otu_df, taxa_are_rows = TRUE)
# Create the new adjusted phyloseq object without phylogenetic tree
physeq_adj <- phyloseq(otu_table_adj, sample_data(physeq), tax_table(physeq))
# Check that negative control filtering worked: 
# Extract count for OTU8 in CS1-10N-S50 from physeq
otu_count_physeq <- otu_table(physeq)["OTU8", "CS1-10N-S50"]
# Print the count
print(paste("Count in physeq: ", otu_count_physeq)) 
# Should be 1007
# Extract count for OTU8 in CS1-10N-S50 from physeq_adj
otu_count_physeq_adj <- otu_table(physeq_adj)["OTU8", "CS1-10N-S50"]
# Print the count
print(paste("Count in physeq_adj: ", otu_count_physeq_adj))
# Should be 1001
# Change the name of physeq_adj to amptk_filtered_nc for downstream analysis
amptk_filtered_nc <- physeq_adj
# Remove Negative Control Samples: 
# Define a vector of sample names to be removed
samples_to_remove <- c("NEG1-S121", "NEG2-S122", "NEG3-S123", "NEG4-S124")
# Reinitialize phyloseq object
amptk_filtered_nc <- subset_samples(amptk_filtered_nc, !phinchID %in% samples_to_remove)
## Within the biom_data file, the "tax_table" that contains all the taxonomic information is incorrectly formated. The column names for different taxonomic levels are "Rank1", "Rank2", "Rank3", etc., when they really should be "kingdom", "phylum", "class", etc... I need to rename these columns appropriately so that the guild assignment step recognizes the required taxonomic headers. 
rank_names(amptk_filtered_nc)
colnames(tax_table(amptk_filtered_nc)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(amptk_filtered_nc)
# Clean up the names in the taxonomy table - get rid of all that junk!
tax_table(amptk_filtered_nc)
tax_table(amptk_filtered_nc)[, colnames(tax_table(amptk_filtered_nc))] <- gsub(tax_table(amptk_filtered_nc)[, colnames(tax_table(amptk_filtered_nc))], pattern = "[a-z]__", replacement = "")
tax_table(amptk_filtered_nc)
```


```{r}
# Load necessary library
library(phyloseq)

# Specify the output directory and file name
output_directory <- "/Users/cbivins/Desktop/Mycoheterotroph_Project_NEATER/Raw_Abundance"
output_file <- "physeq_raw.csv"
output_path <- file.path(output_directory, output_file)

# Extract the OTU table from the phyloseq object
otu_table_raw <- as(otu_table(physeq), "matrix")  # Convert OTU table to a matrix
otu_table_raw_df <- as.data.frame(otu_table_raw)  # Convert to a data frame

# Write the OTU table to a CSV file
write.csv(otu_table_raw_df, file = output_path, row.names = TRUE)

# Confirmation message
cat("OTU table exported to:", output_path, "\n")

```


#### Assign functional guilds
```{r}
# Guild assignment function 
assign_guild <- function(object, database) { 
  ### This function assigns trophic modes (and other traits) to OTUs based on taxonomy
  # returning a table with taxonomy and guild/traits for each OTU 
  ### Arguments:
  # object: a phyloseq object, for example imported from a .biom file
  # database: Reference used, for example FungalTraits (Polme et al. 2020, Fungal Diversity 105, 1-16).
  ### Function:
  # load required packages 
  require(file2meco)
  require(microeco)
  # convert phyloseq object into a "microtable"
  meco <- phyloseq2meco(object)
  # verify that OTUs and samples information is consistent across files
  meco$tidy_dataset()
  # assign guilds
  t1 <- trans_func$new(dataset = meco)
  t1$for_what <- "fungi"
  t1$cal_spe_func(fungi_database = database)
  # create a dataframe with taxonomy and guild/traits information
  as.data.frame(t1$res_spe_func_raw_FungalTraits)
}
# Assign guilds to OTUs
# I'm getting an error message saying "Error in if (any(apply(otu_table, 1, sum) == 0)) { :
# I need to first investigate what values have either NA or NaN before I remove them (or change them to zeros), as will likely need to be done in order to assign guilds 
# Check for NA values in the OTU table
any(is.na(otu_table(amptk_filtered_nc)))
# Convert the OTU table to a matrix
otu_matrix <- as.matrix(otu_table(amptk_filtered_nc))
# Find the indices of NA values
na_indices <- which(is.na(otu_matrix), arr.ind = TRUE)
# Print the indices
print(na_indices)
# Ah, the Negative control samples (NEG1-4) have all been changed to NaN - we just need to convert these to zeroes
# Get the OTU table
otu_table <- otu_table(amptk_filtered_nc)
# Replace NaN values with 0
otu_table[is.nan(otu_table)] <- 0
# Assign the modified otu_table back to the phyloseq object
otu_table(amptk_filtered_nc) <- otu_table
# Assign guilds to the entire dataset
guild_table = assign_guild(object = amptk_filtered_nc, database = "FungalTraits")
```

### Create a new phyloseq object with OTU counts transformed into relative abundances
```{r}
# Transform the OTU table in amptk_filtered_nc into relative abundances
MH_phyloseq_RA <- transform_sample_counts(
  amptk_filtered_nc,
  function(x) x / sum(x)
)

# Quick check of the transformation:
# (1) Are the sums of OTU counts now 1 for each sample?
sample_sums(MH_phyloseq_RA)

```

### Export guild table, relative abundance otu table, and sample table as CSV files
```{r}
# Define output directory
out_dir <- "/Users/christopherbivins/Desktop/Mycoheterotroph_Project_Neater/Relative_Abundance"

# 1) Export the guild_table data frame
tryCatch({
  write.csv(guild_table,
            file = file.path(out_dir, "guild_table.csv"),
            row.names = TRUE)
  message("Successfully saved guild_table.csv")
}, error = function(e) {
  message("ERROR saving guild_table.csv: ", e$message)
})

# 2) Export the OTU table from MH_phyloseq_RA
#    Convert the phyloseq OTU table to a data frame for writing to CSV.
tryCatch({
  write.csv(as.data.frame(otu_table(MH_phyloseq_RA)),
            file = file.path(out_dir, "MH_phyloseq_RA_otu_table.csv"),
            row.names = TRUE)
  message("Successfully saved MH_phyloseq_RA_otu_table.csv")
}, error = function(e) {
  message("ERROR saving MH_phyloseq_RA_otu_table.csv: ", e$message)
})

# 3) Sample data
# Forcefully convert to data frame
sample_data_df <- as.data.frame(as.matrix(sample_data(MH_phyloseq_RA)))

# Specify the output path
output_file <- "/Users/christopherbivins/Desktop/Mycoheterotroph_Project_Neater/Relative_Abundance/sample_data_MH_phyloseq_RA.csv"

# Export the data frame
write.csv(sample_data_df, file = output_file, row.names = TRUE)

# Confirmation message
cat("Sample data has been exported to:", output_file)


```
