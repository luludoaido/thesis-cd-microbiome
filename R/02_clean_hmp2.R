# --- 02 Clean HMP2 ---
# Turns the HMP2 metadata and bioBakery 3 profiles into clean tables with one row per sample.

# Results: metadata of adult CD and controls with HBI, calprotectin and antibiotics;
# species, genus, taxonomy, butyrate pathways and enzymes

# Methods: docs/methods.md, 1.1 and section 4

library(data.table)
library(tidyverse)
source("R/00_setup.R")

out_dir <- file.path(dir_processed, "hmp2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Function to turn a table with features in rows and samples in columns into one row per sample
features_by_sample <- function(tab, suffix) {
  mat <- t(as.data.frame(tab)[, -1])
  colnames(mat) <- tab[[1]]
  data.frame(sample_id = sub(suffix, "", rownames(mat)), mat,
             check.names = FALSE, row.names = NULL)
}
# --- Metadata ---
hmp2_meta_raw <- read.csv(file.path(dir_raw, "hmp2/hmp2_metadata_2018-08-20.csv"), sep = ";", check.names = FALSE)

# metagenome samples, CD and controls, adults with known age only (see methods 4.2, 4.3)
hmp2_meta <- hmp2_meta_raw[hmp2_meta_raw$data_type == "metagenomics" &
                             hmp2_meta_raw$diagnosis %in% c("CD", "nonIBD") &
                             !is.na(hmp2_meta_raw$consent_age) &
                             hmp2_meta_raw$consent_age >= hmp2_min_age,
                           c("Project", "External ID", "Participant ID", "site_sub_coll", "week_num",
                             "interval_days", "visit_num", "consent_age", "diagnosis", "reads_filtered")]
hmp2_meta$diagnosis[hmp2_meta$diagnosis == "nonIBD"] <- "control"

# --- Clinical values ---
# same sample, else the nearest visit within max_gap weeks (see methods 4.4)
max_gap <- 2

hmp2_visits <- hmp2_meta_raw %>% transmute(visit_id = `External ID`, `Participant ID`, visit_week = week_num,
                                           hbi = suppressWarnings(as.numeric(hbi)),
                                           calpro = suppressWarnings(as.numeric(fecalcal)),
                                           abx = na_if(trimws(`Has the subject used any antibiotics since the last visit?`), ""))

# finds the best value of one variable for each sample and adds three columns:
# the value, where it comes from (_how) and the distance in weeks (_gap)
match_var <- function(var) {
  hmp2_meta %>%
    select(`External ID`, `Participant ID`, week_num) %>%
    inner_join(filter(hmp2_visits, !is.na(.data[[var]])), by = "Participant ID", relationship = "many-to-many") %>%
    mutate(same_sample = visit_id == `External ID`, gap = abs(week_num - visit_week)) %>%
    filter(same_sample | gap <= max_gap) %>%
    arrange(`External ID`, desc(same_sample), gap) %>%
    distinct(`External ID`, .keep_all = TRUE) %>%
    transmute(`External ID`, value = .data[[var]], how = case_when(same_sample ~"same sample",
                                                                   gap == 0 ~ "same visit",
                                                                   TRUE ~"nearest visit"),
              gap = if_else(same_sample, 0, gap)) %>%
    set_names("External ID", var, paste0(var, "_how"),
              paste0(var, "_gap"))
}

hmp2_meta <- hmp2_meta %>%
  left_join(match_var("hbi"), by = "External ID") %>%
  left_join(match_var("calpro"), by = "External ID") %>%
  left_join(match_var("abx"), by = "External ID")

# MSM6J2MF has no enzyme data (see methods 4.5)
hmp2_meta <- hmp2_meta %>%
  filter(`External ID` != "MSM6J2MF") %>%
  transmute(sample_id = `External ID`, subject_id = `Participant ID`, diagnosis, time = week_num,
            age = consent_age, antibiotic = NA_character_, site = substr(`Participant ID`, 1, 1),
            reads = reads_filtered, Project, site_sub_coll, interval_days, visit_num, hbi, hbi_how, hbi_gap, calpro, calpro_how, calpro_gap, abx, abx_how, abx_gap) %>%
  finish_meta()

write_csv(hmp2_meta, file.path(out_dir, "hmp2_meta.csv"))

# --- Species and genus (MetaPhlAn 3) ---
# Each row is one taxon with its full lineage.
# Species rows contain "s__" but no strain ("t__"); genus rows contain "g__" but no "s__".
taxa <- fread(file.path(dir_raw, "hmp2/taxonomic_profiles_3.tsv"))

species <- features_by_sample(taxa[grepl("s__", taxa[[1]]) & !grepl("t__", taxa[[1]]), ], "_profile$")

# order, family and genus of each species
lineage <- colnames(species)[-1]
hmp2_taxonomy <- data.frame(species = clean_species(sub(".*s__", "", lineage)),
                            genus   = str_match(lineage, "g__([^|]+)")[, 2],
                            family  = str_match(lineage, "f__([^|]+)")[, 2],
                            order   = str_match(lineage, "o__([^|]+)")[, 2])

write_csv(hmp2_taxonomy, file.path(out_dir, "hmp2_taxonomy.csv"))

colnames(species)[-1] <- hmp2_taxonomy$species
save_table(species, hmp2_meta, "hmp2", "species")

genus <- features_by_sample(taxa[grepl("g__", taxa[[1]]) & !grepl("s__", taxa[[1]]), ], "_profile$")
colnames(genus)[-1] <- sub(".*g__", "", colnames(genus)[-1])
save_table(genus, hmp2_meta, "hmp2", "genus")

# --- Butyrate pathways and enzymes (HUMAnN 3) ---
# total rows only (see methods 4.8)
butanoate_pathways <- c("CENTFERM-PWY: pyruvate fermentation to butanoate",
                        "PWY-5676: acetyl-CoA fermentation to butanoate II",
                        "PWY-5677: succinate fermentation to butanoate")

butyrate_ecs <- c("2.3.1.19: Phosphate butyryltransferase",
                  "2.7.2.7: Butyrate kinase",
                  "2.8.3.8: Acetate CoA-transferase")

pathways <- fread(file.path(dir_raw, "hmp2/pathabundances_3.tsv"), header = TRUE)
pathways <- features_by_sample(pathways[pathways[[1]] %in% butanoate_pathways], "_pathabundance_cpm$")
save_table(pathways, hmp2_meta, "hmp2", "pathways")

ecs <- fread(file.path(dir_raw, "hmp2/ecs_3.tsv"), header = TRUE)   # 700 MB, takes a moment
ecs <- ecs[ecs[[1]] %in% butyrate_ecs]
ecs[, (2:ncol(ecs)) := lapply(.SD, as.numeric), .SDcols = 2:ncol(ecs)]
ecs <- features_by_sample(ecs, "_level4ec$")
save_table(ecs, hmp2_meta, "hmp2", "enzymes")
