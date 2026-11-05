# --- Data Cleaning: Franzosa ---
library(readxl)
library(tidyverse)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "franzosa") # Output folder
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#function to turn a block of feature rows into one row per sample
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

# CD and controls only, Franzosa has one sample per person and no read counts
franzosa_meta <- meta_raw %>%
  filter(Diagnosis != "UC") %>%
  transmute(sample_id, subject_id = sample_id, diagnosis = Diagnosis, time = NA_real_,
            age = as.numeric(Age), antibiotic, site = sub("[|].*", "", sample_id),
            reads = NA_real_, 
            calprotectin = as.numeric(Fecal.Calprotectin), immunosuppressant, mesalamine, steroids, SRA_metagenome_name) %>%
  finish_meta()

write_csv(franzosa_meta, file.path(out_dir, "franzosa_meta.csv"))

# --- Species (from row 11 on Supplementary Dataset 4) ---
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

# --- Enzymes (from row 11 of Supplementary Data set 6) ---
raw_enzymes <- read_excel(file.path(dir_raw, "franzosa/SupplDataset6_enzymes.xlsx"), col_names = FALSE)
enzymes <- features_by_sample(raw_enzymes[11:nrow(raw_enzymes), ])
save_table(enzymes, franzosa_meta, "franzosa", "enzymes")

# --- Metabolites (from row 10 Supplementary Data set 2) ---
raw_metab <- read_excel(file.path(dir_raw, "franzosa/SupplDataset2_metabolites.xlsx"),
                        col_names =  FALSE)
metabolites <- features_by_sample(raw_metab[10:nrow(raw_metab), ])
save_table(metabolites, franzosa_meta, "franzosa", "metabolites")

# --- Butyrate and Propionate ---
# Data set 1 annotates the clusters: 
# Cluster_0013 = butyrate/isobutyrate (not separated)
# Cluster_0003 = propionate

metab_info <- read_excel(file.path(dir_raw, "franzosa/SupplDataset1_metabolite_info.xlsx"), skip = 1)
colnames(metab_info) <- c("feature", "retention_time", "mz", "cluster_da", "chemical_class", "match_to_standard", "adduct")
scfa_lookup <- metab_info[metab_info$feature %in% c("HILIC-neg_Cluster_0013", "HILIC-neg_Cluster_0003"), ]

write.csv(scfa_lookup, file.path(dir_processed, "franzosa", "franzosa_scfa_lookup.csv"), row.names = FALSE)

scfa <- metabolites[, c("sample_id", "HILIC.neg_Cluster_0013", "HILIC.neg_Cluster_0003")]
colnames(scfa) <- c("sample_id", "butyrate", "propionate")
save_table(scfa, franzosa_meta, "franzosa", "scfa")
