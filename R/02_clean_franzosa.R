# --- Data Cleaning: Franzosa ---
library(readxl)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "franzosa") # Output folder

#function to turn a block of feature rows into one row per sample
features_by_sample <- function(block) {
  df <- as.data.frame(t(block[, -1]))
  colnames(df) <- block[[1]]
  df <- as.data.frame(lapply(df, as.numeric))
  data.frame(sample_id = sample_ids, df, check.names = FALSE)
}

# --- Metadata (rows 3-10) ---
raw_species <- read_excel(file.path(dir_raw, "franzosa/supplDataset4_species.xlsx"),
                          col_names = FALSE)

sample_ids <- as.character(raw_species[2, -1])

meta_block <- raw_species[3:10, ]
franzosa_meta <- as.data.frame(t(meta_block[,-1]))
colnames(franzosa_meta) <- meta_block[[1]]
# sample Id as first column
franzosa_meta <- data.frame(sample_id = sample_ids, franzosa_meta, check.names = FALSE)
franzosa_meta$Age <- as.numeric(franzosa_meta$Age)
franzosa_meta$Fecal.Calprotectin <- as.numeric(franzosa_meta$Fecal.Calprotectin)

#CD and Control only
franzosa_meta <- franzosa_meta[franzosa_meta$Diagnosis != "UC", ]
write.csv(franzosa_meta, file.path(out_dir, "franzosa_meta.csv"), row.names = FALSE)

# --- Species (from row 11 on Supplementary Dataset 4) ---
species <- features_by_sample(raw_species[11:nrow(raw_species), ])
write.csv(merge(franzosa_meta, species, by = "sample_id"), 
          file.path(out_dir, "franzosa_features_merged.csv"), 
          row.names = FALSE)

# --- Enzymes (from row 11 of Supplementary Data set 6) ---
raw_enzymes <- read_excel(file.path(dir_raw, "franzosa/SupplDataset6_enzymes.xlsx"), col_names = FALSE)
enzymes <- features_by_sample(raw_enzymes[11:nrow(raw_enzymes), ])
write.csv(merge(franzosa_meta, enzymes, by = "sample_id"),
          file.path(out_dir, "franzosa_enz_merged.csv"), row.names = FALSE)

# --- Metabolites (from row 10 Supplementary Data set 2) ---
raw_metab <- read_excel(file.path(dir_raw, "franzosa/SupplDataset2_metabolites.xlsx"),
                        col_names =  FALSE)
metabolites <- features_by_sample(raw_metab[10:nrow(raw_metab), ])

write.csv(merge(franzosa_meta, metabolites, by = "sample_id"),
          file.path(out_dir,
                    "franzosa_metab_merged.csv"), row.names = FALSE)

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

write.csv(merge(franzosa_meta, scfa, by = "sample_id"),
          file.path(out_dir, "franzosa_h3a_merged.csv"), row.names = FALSE)
