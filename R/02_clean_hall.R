# --- 02 Clean Hall 2017 ---
# Turns the raw Hall files from 01_download_cmd.R into clean tables with one row per sample.

# Results: metadata of CD and controls; species, genus and taxonomy per species

# Methods: docs/methods.md, 1.1 and section 3

library(tidyverse)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "hall")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hall_meta_raw <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_metadata.csv"))
hall_taxa <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_taxa.csv"), check.names = FALSE)
hall_taxonomy <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_taxonomy.csv"))

# --- Metadata ---
# controls get "control" (see methods 3.2), centre (see methods 3.3), time in weeks (see methods 3.4)
hall_meta <- hall_meta_raw %>%
  transmute(sample_id = X, subject_id, diagnosis = ifelse(study_condition == "control",
                                                          "control", disease_subtype),
            time = days_from_first_collection / 7,
            age = NA_real_, antibiotic = NA_character_, site = ifelse(grepl("^SKST", subject_id),
                                                                      "Atlanta", "Boston"),
            reads = number_reads,
            age_category, country, calprotectin, HBI, days_from_first_collection, visit_number) %>%
  filter(diagnosis %in% c("CD", "control"))

selected_samples <- hall_meta$sample_id

hall_meta <- finish_meta(hall_meta)
write_csv(hall_meta, file.path(out_dir, "hall_meta.csv"))

# --- Species ---
species_rows <- hall_taxa[grepl("^species:", hall_taxa[[1]]), ]
species_values <- as.matrix(as.data.frame(species_rows)[, -1])  # species in rows, samples in columns

species <- as.data.frame(t(species_values))
colnames(species) <- clean_species(gsub("^species:", "", species_rows[[1]]))

# rename to the names used in HMP2 and Franzosa (see methods 3.5)
species_synonyms <- c(
  Phocaeicola_vulgatus     = "Bacteroides_vulgatus",
  Phocaeicola_dorei        = "Bacteroides_dorei",
  Phocaeicola_massiliensis = "Bacteroides_massiliensis",
  Phocaeicola_coprocola    = "Bacteroides_coprocola",
  Phocaeicola_plebeius     = "Bacteroides_plebeius",
  Phocaeicola_coprophilus  = "Bacteroides_coprophilus",
  Lachnospira_eligens      = "Eubacterium_eligens",
  Segatella_copri          = "Prevotella_copri"
)

colnames(species) <- ifelse(colnames(species) %in%
                              names(species_synonyms),
                            species_synonyms[colnames(species)],
                            colnames(species))

# order, family and genus of each species
taxonomy_rows <-
  hall_taxonomy[match(species_rows[[1]],
                      hall_taxonomy[[1]]), ]
hall_taxonomy_clean <- data.frame(species = colnames(species), genus = taxonomy_rows$genus,
                                  family = taxonomy_rows$family, order = taxonomy_rows$order)
write_csv(distinct(hall_taxonomy_clean, species, .keep_all = TRUE), file.path(out_dir, "hall_taxonomy.csv"))

# columns with the same name after renaming are added up
species <- data.frame(sample_id = rownames(species),
                      species, check.names = FALSE, row.names = NULL) %>%
  pivot_longer(-sample_id, names_to = "taxon", values_to = "abundance") %>%
  group_by(sample_id, taxon) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>% pivot_wider(names_from = taxon, values_from = abundance)

save_table(species, hall_meta, "hall", "species")

# --- Genus ---
# sums the species of each genus; species without a genus are dropped here (see methods 3.6)

genus_lookup <- hall_taxonomy$genus[match(species_rows[[1]],
                                          hall_taxonomy[[1]])]

has_genus <- !is.na(genus_lookup)
genus_sums <- rowsum(species_values[has_genus, ],
                     group = genus_lookup[has_genus], na.rm = TRUE)

# 3 rows are classified only to genus level (see methods 3.7)
genus_only <- hall_taxa[grepl("^genus:", hall_taxa[[1]]), ]
genus_only_values <- as.matrix(as.data.frame(genus_only)[, -1])
rownames(genus_only_values) <- gsub("^genus:", "", genus_only[[1]])

for (g in rownames(genus_only_values)) {
  if (g %in% rownames(genus_sums)) {
    genus_sums[g, ] <- genus_sums[g, ] + genus_only_values[g, ]
  } else {
    genus_sums <- rbind(genus_sums, genus_only_values[g, , drop = FALSE])
  }
}

genus <- as.data.frame(t(genus_sums))
genus <- data.frame(sample_id = rownames(genus), genus, check.names = FALSE, row.names = NULL)
save_table(genus, hall_meta, "hall", "genus")

# all CD and control samples before the exclusions, for R/checks/check_hall_samples.R
write_csv(genus[genus$sample_id %in% selected_samples, ], file.path(out_dir, "hall_genus_before_exclusions.csv"))
