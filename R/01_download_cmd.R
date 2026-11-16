# --- 01_download_cmd ---

# What this script does
# Downloads the Hall 2017 data from curatedMetagenomicData (Pasolli 2017) and saves them
# unchanged, so that the original data stay available even if the database changes.
# HMP2 and Franzosa are downloaded by hand.

# Decisions
# 1. A fixed database version (2021-03-31), so that a later download gives the same data.
# 2. Relative abundances instead of estimated counts, as for the other cohorts.

library(curatedMetagenomicData)   # Bioconductor
source("R/00_setup.R")

# fixed database version
hall <- curatedMetagenomicData("2021-03-31.HallAB_2017.relative_abundance",
                               dryrun = FALSE, # to load the dataset
                               counts = FALSE, # relative abundances, not counts
                               rownames = "short")[[1]] # species names without the full lineage

out <- file.path(dir_raw, "cmd")

write.csv(assay(hall), file.path(out, "HallAB_2017_taxa.csv"))
write.csv(as.data.frame(colData(hall)), file.path(out, "HallAB_2017_metadata.csv"))
write.csv(as.data.frame(rowData(hall)), file.path(out, "HallAB_2017_taxonomy.csv"))

# --- Literature ---
# Hall AB et al. (2017). A novel Ruminococcus gnavus clade enriched in inflammatory bowel disease
#   patients. Genome Med 9:103. doi:10.1186/s13073-017-0490-5
# Pasolli E et al. (2017). Accessible, curated metagenomic data through ExperimentHub.
#   Nat Methods 14(11):1023-1024. doi:10.1038/nmeth.4468
