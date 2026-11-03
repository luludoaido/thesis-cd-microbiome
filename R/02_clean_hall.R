# --- Data Cleaning: Hall 2017 ---
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "hall")

hall_meta_raw <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_metadata.csv"))
hall_taxa <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_taxa.csv"), check.names = FALSE)
hall_taxonomy <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_taxonomy.csv"))

# --- Metadata ---
# In CuratedMetagenomicData controls have no disease_subtype, so they get "control"
hall_meta <- data.frame(
  ID           = hall_meta_raw$X,
  subject_id   = hall_meta_raw$subject_id,
  diagnosis    = ifelse(hall_meta_raw$study_condition == "control", "control", hall_meta_raw$disease_subtype),
  age_category = hall_meta_raw$age_category,
  country      = hall_meta_raw$country,
  calprotectin = hall_meta_raw$calprotectin,
  HBI          = hall_meta_raw$HBI,
  days_from_first_collection = hall_meta_raw$days_from_first_collection,
  visit_number = hall_meta_raw$visit_number
)

#CD and controls only
hall_meta <- hall_meta[hall_meta$diagnosis %in% c("CD", "control"), ]
write.csv(hall_meta, file.path(out_dir, "hall_meta.csv"), row.names = FALSE)

# rename column name "ID" to "sample_id" for merge()
meta_for_merge <- hall_meta
names(meta_for_merge)[names(meta_for_merge) == "ID"] <- "sample_id"

# --- Species ---
species_rows <- hall_taxa[grepl("^species:", hall_taxa[[1]]), ]
species_values <- as.matrix(as.data.frame(species_rows)[, -1])
#species x samples

species <- as.data.frame(t(species_values))
colnames(species) <- gsub("^species:", "", species_rows[[1]])
species <- data.frame(sample_id = rownames(species), species, check.names = FALSE, row.names = NULL)

write.csv(merge(meta_for_merge, species, by = "sample_id"), 
          file.path(out_dir, "hall_species_merged.csv"),
          row.names = FALSE)

# --- Genus ---
# Sum the species of each genus - there are 23 species that have no genus in the taxonomy
# those species are dropped here and kept at the species level.

genus_lookup <- hall_taxonomy$genus[match(species_rows[[1]],
                                          hall_taxonomy[[1]])]

has_genus <- !is.na(genus_lookup)
genus_sums <- rowsum(species_values[has_genus, ],
                     group = genus_lookup[has_genus], na.rm = TRUE)

# there are 3 rows that are classified only to genus level
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

write.csv(merge(meta_for_merge, genus, by = "sample_id"), 
          file.path(out_dir, "hall_genus_merged.csv"),
          row.names = FALSE)

