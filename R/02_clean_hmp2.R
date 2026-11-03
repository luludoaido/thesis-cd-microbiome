# --- Data cleaning: HMP2 ---
library(data.table)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "hmp2")

# Function to turn a table with features in rows and samples in columns into one row per sample
features_by_sample <- function(tab, suffix) {
  mat <- t(as.data.frame(tab)[, -1])
  colnames(mat) <- tab[[1]]
  data.frame(sample_id = sub(suffix, "", rownames(mat)), mat, 
             check.names = FALSE, row.names = NULL)
}
# --- Metadata ---
hmp2_meta_raw <- read.csv(file.path(dir_raw, "hmp2/hmp2_metadata_2018-08-20.csv"), sep = ";", check.names = FALSE)

# Meta-genome samples
hmp2_meta <- hmp2_meta_raw[hmp2_meta_raw$data_type == "metagenomics" &
                             hmp2_meta_raw$diagnosis %in% c("CD", "nonIBD") &
                             !is.na(hmp2_meta_raw$consent_age)& 
                             hmp2_meta_raw$consent_age >= hmp2_min_age, 
                                    c("Project", "External ID", "Participant ID",
                                      "site_sub_coll", "week_num", "interval_days",
                                      "visit_num", "consent_age", "diagnosis")]

hmp2_meta$diagnosis[hmp2_meta$diagnosis == "nonIBD"] <- "control"

write.csv(hmp2_meta, file.path(out_dir, "hmp2_meta.csv"), row.names = FALSE)

# Rename column "External ID" to "sample_id" for merge
# MSM6J2MF has no EC data. so exclusion
meta_for_merge <- hmp2_meta
colnames(meta_for_merge)[colnames(meta_for_merge) == "External ID"] <- "sample_id"
meta_for_merge <- meta_for_merge[meta_for_merge$sample_id != "MSM6J2MF", ]


# --- Species and genus ---
# Species rows contain "s__" but no strain ("t__"); genus row contain "g__" but no "s__"
taxa <- fread(file.path(dir_raw, "hmp2/taxonomic_profiles_3.tsv"))
species <- features_by_sample(taxa[grepl("s__", taxa[[1]]) & !grepl("t__", taxa[[1]]), ], "_profile$")
genus <- features_by_sample(taxa[grepl("g__", taxa[[1]]) & !grepl("s__", taxa[[1]]), ], "_profile$")

write.csv(merge(meta_for_merge, species, by = "sample_id"), file.path(out_dir, "hmp2_species_merged.csv"),
          row.names = FALSE)
write.csv(merge(meta_for_merge, genus, by = "sample_id"), file.path(out_dir, "hmp2_genus_merged.csv"),
          row.names = FALSE)


# --- Butyrate pathways and enzymes ---
butanoate_pathways <- c("CENTFERM-PWY: pyruvate fermentation to butanoate", 
                        "PWY-5676: acetyl-CoA fermentation to butanoate II",
                        "PWY-5677: succinate fermentation to butanoate")

butyrate_ecs <- c("2.3.1.19: Phosphate butyryltransferase",
                  "2.7.2.7: Butyrate kinase",
                  "2.8.3.8: Acetate CoA-transferase")

pathways <- fread(file.path(dir_raw, "hmp2/pathabundances_3.tsv"), header = TRUE)
pathways <- features_by_sample(pathways[pathways[[1]] %in% butanoate_pathways], "_pathabundance_cpm$")

ecs <- fread(file.path(dir_raw, "hmp2/ecs_3.tsv"), header = TRUE)
ecs <- ecs[ecs[[1]] %in% butyrate_ecs]
ecs[, (2:ncol(ecs)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(ecs)]
ecs <- features_by_sample(ecs, "_level4ec$")

write.csv(merge(meta_for_merge, pathways, by = "sample_id"), 
          file.path(out_dir,"hmp2_pathway_butanoate_merged.csv"),
          row.names = FALSE)
write.csv(merge(meta_for_merge, ecs, by = "sample_id"),
          file.path(out_dir, "hmp2_ec_butyrate_merged.csv"),
          row.names = FALSE)


