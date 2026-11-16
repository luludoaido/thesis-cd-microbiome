# --- Data Cleaning: Franzosa ---

# What this script does
# Reads the Supplementary Datasets of Franzosa 2019 and turns them into clean tables with one
# row per sample: metadata, species, genus, enzymes, metabolites, butyrate and propionate.
# The species profiles were generated with MetaPhlAn 2 (Truong 2015).

# Decisions
# 1. CD and controls only, because the thesis compares these two groups. UC samples are removed.
# 2. Each sample is one person, so sample_id is also used as subject_id.
# 3. The sub-cohort (PRISM or Validation) is used as the centre.
# 4. A genus is the sum of its species.
# 5. Butyrate and isobutyrate cannot be separated in this dataset, so they stay together as "butyrate".

library(readxl)
library(tidyverse)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "franzosa") 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# function to turn a block of feature rows into one row per sample
features_by_sample <- function(block) {
  df <- as.data.frame(t(block[, -1]))
  colnames(df) <- block[[1]]
  df <- as.data.frame(lapply(df, as.numeric))
  data.frame(sample_id = sample_ids, df, check.names = FALSE)
}

# --- Metadata (rows 3-10) ---
raw_species <- read_excel(file.path(dir_raw, "franzosa/SupplDataset4_species.xlsx"),
                          col_names = FALSE)

sample_ids <- as.character(raw_species[2, -1])

meta_block <- raw_species[3:10, ]
meta_raw <- as.data.frame(t(meta_block[,-1]))
colnames(meta_raw) <- meta_block[[1]]
meta_raw$sample_id <- sample_ids

# CD and controls only. Franzosa has one sample per person and no read counts,
# so time and reads stay empty and the read limit from 00_setup.R does not apply.
franzosa_meta <- meta_raw %>%
  filter(Diagnosis != "UC") %>%
  transmute(sample_id, subject_id = sample_id, diagnosis = Diagnosis, time = NA_real_,
            age = as.numeric(Age), antibiotic, site = sub("[|].*", "", sample_id),
            reads = NA_real_,
            calprotectin = as.numeric(Fecal.Calprotectin), immunosuppressant, mesalamine, steroids, SRA_metagenome_name) %>%
  finish_meta()

write_csv(franzosa_meta, file.path(out_dir, "franzosa_meta.csv"))

# --- Species (from row 11 of Supplementary Dataset 4) ---
species <- features_by_sample(raw_species[11:nrow(raw_species), ])

# --- Genus ---
# sum the species of each genus (genus = everything before the first "_")
genus <- species %>%
  pivot_longer(-sample_id, names_to = "species", values_to = "abundance") %>%
  mutate(genus = sub("_.*", "", species)) %>%
  group_by(sample_id, genus) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = abundance)
save_table(genus, franzosa_meta, "franzosa", "genus")

colnames(species)[-1] <- clean_species(colnames(species)[-1])
save_table(species, franzosa_meta, "franzosa", "species")

# --- Enzymes (from row 11 of Supplementary Dataset 6) ---
raw_enzymes <- read_excel(file.path(dir_raw, "franzosa/SupplDataset6_enzymes.xlsx"), col_names = FALSE)
enzymes <- features_by_sample(raw_enzymes[11:nrow(raw_enzymes), ])
save_table(enzymes, franzosa_meta, "franzosa", "enzymes")

# --- Metabolites (from row 10 of Supplementary Dataset 2) ---
raw_metab <- read_excel(file.path(dir_raw, "franzosa/SupplDataset2_metabolites.xlsx"),
                        col_names =  FALSE)
metabolites <- features_by_sample(raw_metab[10:nrow(raw_metab), ])
save_table(metabolites, franzosa_meta, "franzosa", "metabolites")

# --- Butyrate and Propionate ---
# Dataset 1 annotates the clusters:
# Cluster_0013 = butyrate/isobutyrate (not separated)
# Cluster_0003 = propionate
# The values are peak intensities, not concentrations.

metab_info <- read_excel(file.path(dir_raw, "franzosa/SupplDataset1_metabolite_info.xlsx"), skip = 1)
colnames(metab_info) <- c("feature", "retention_time", "mz", "cluster_da", "chemical_class", "match_to_standard", "adduct")
scfa_lookup <- metab_info[metab_info$feature %in% c("HILIC-neg_Cluster_0013", "HILIC-neg_Cluster_0003"), ]

write.csv(scfa_lookup, file.path(dir_processed, "franzosa", "franzosa_scfa_lookup.csv"), row.names = FALSE)

scfa <- metabolites[, c("sample_id", "HILIC.neg_Cluster_0013", "HILIC.neg_Cluster_0003")]
colnames(scfa) <- c("sample_id", "butyrate", "propionate")
save_table(scfa, franzosa_meta, "franzosa", "scfa")

# --- Literature ---
# Franzosa EA et al. (2019). Gut microbiome structure and metabolic activity in inflammatory bowel
#   disease. Nat Microbiol 4(2):293-305. doi:10.1038/s41564-018-0306-4
# Truong DT et al. (2015). MetaPhlAn2 for enhanced metagenomic taxonomic profiling.
#   Nat Methods 12(10):902-903. doi:10.1038/nmeth.3589
