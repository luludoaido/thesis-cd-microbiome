# --- 01 Download Hall 2017 ---
# Downloads Hall 2017 from curatedMetagenomicData and saves it unchanged.
# HMP2 and Franzosa are downloaded by hand.

# Results: raw Hall tables (taxa, metadata, taxonomy), fixed database version

# Methods: docs/methods.md, 3.1

library(curatedMetagenomicData)   # Bioconductor
source("R/00_setup.R")

# fixed database version (see methods 3.1)
hall <- curatedMetagenomicData("2021-03-31.HallAB_2017.relative_abundance",
                               dryrun = FALSE, # to load the dataset
                               counts = FALSE, # relative abundances, not counts
                               rownames = "short")[[1]] # species names without the full lineage

out <- file.path(dir_raw, "cmd")

write.csv(assay(hall), file.path(out, "HallAB_2017_taxa.csv"))
write.csv(as.data.frame(colData(hall)), file.path(out, "HallAB_2017_metadata.csv"))
write.csv(as.data.frame(rowData(hall)), file.path(out, "HallAB_2017_taxonomy.csv"))
