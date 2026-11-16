# --- Data cleaning: HMP2 ---

# What this script does
# Reads the HMP2 metadata (Lloyd-Price 2019) and the bioBakery 3 profiles (Beghini 2021) and turns
# them into clean tables with one row per sample: metadata with clinical values, species, genus,
# taxonomy, and the butyrate pathways and enzymes.

# Decisions
# 1. Metagenome samples of CD and non-IBD controls only. nonIBD is called control.
# 2. Adults with a known age only, for comparability with the other cohorts (own decision).
# 3. HBI, calprotectin and antibiotics come from the same sample, or else from the nearest visit
#    of the same person at most 2 weeks away (own decision).
# 4. Sample MSM6J2MF is removed, because it has no enzyme data and all tables should contain the same samples.
# 5. The centre is the first letter of the Participant ID.
# 6. Butyrate is described by three pathways and three enzymes, using only the total rows.

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

# metagenome samples, CD and controls, adults with known age only
hmp2_meta <- hmp2_meta_raw[hmp2_meta_raw$data_type == "metagenomics" &
                             hmp2_meta_raw$diagnosis %in% c("CD", "nonIBD") &
                             !is.na(hmp2_meta_raw$consent_age) &
                             hmp2_meta_raw$consent_age >= hmp2_min_age,
                           c("Project", "External ID", "Participant ID", "site_sub_coll", "week_num",
                             "interval_days", "visit_num", "consent_age", "diagnosis", "reads_filtered")]
hmp2_meta$diagnosis[hmp2_meta$diagnosis == "nonIBD"] <- "control"

# --- Clinical values ---
# HBI, calprotectin and antibiotics are often recorded at a different visit than the stool sample.
# For each sample, the value comes from the same sample if there is one, otherwise from the nearest
# visit of the same person within max_gap weeks. The 2 weeks are an own decision;
# 03_eda_clinical.R repeats the comparison with 1 and 0 weeks.
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

# MSM6J2MF has no enzyme data and is removed, so that all tables contain the same samples
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
# Butyrate production is commonly studied through the genes of its pathways. Single genes such as but
# and buk can be misannotated, so whole pathways are more robust (Vital 2014).
# The selection of these pathways and enzymes is an own decision; succinate fermentation is not covered by Vital 2014.
# Only the total rows are used; rows split by species are not.
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

# --- Literature ---
# Beghini F et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse
#   microbial communities with bioBakery 3. eLife 10:e65088. doi:10.7554/eLife.65088
# Lloyd-Price J et al. (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel
#   diseases. Nature 569(7758):655-662. doi:10.1038/s41586-019-1237-9
# Vital M et al. (2014). Revealing the bacterial butyrate synthesis pathways by analyzing
#   (meta)genomic data. mBio 5(2):e00889-14. doi:10.1128/mbio.00889-14
