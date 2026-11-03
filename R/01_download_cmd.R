# --- 01_download_cmd ---
# Download Hall 2017 from curatedMetagenomic Data and save it unchanged as 
# three CSV files in data/raw/cmd/
# HMP2 and Franzosa are downloaded by hand

library(curatedMetagenomicData)   # Bioconductor
source("R/00_setup.R")

#fixed database version
hall <- curatedMetagenomicData("2021-03-31.HallAB_2017.relative_abundance",
                               dryrun = FALSE, #to load the dataset
                               counts = FALSE, #bc we want relative abundance
                               rownames = "short")[[1]] # species names without the full lineage

out <- file.path(dir_raw, "cmd")

write.csv(assay(hall), file.path(out, "HallAB_2017_taxa.csv"))
write.csv(as.data.frame(colData(hall)), file.path(out, "HallAB_2017_metadata.csv"))
write.csv(as.data.frame(rowData(hall)), file.path(out, "HallAB_2017_taxonomy.csv"))
