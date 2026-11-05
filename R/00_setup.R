# --- 00_Set up ---

# Shared settings for all scripts
# Every script starts with source("R/00_setup.R"), so paths, exclusion rules,
# seed and colors are the same everywhere.

# Run the scripts from the project folder (open thesis-cd-microbiome.Rproj in Rstudio)

library(ggplot2)

if (!file.exists("thesis-cd-microbiome.Rproj")) {
  stop("Working directory is not the project folder (", getwd(), ").\n",
       "Open thesis-cd-microbiome.Rproj in RStudio or run Rscript from the project folder.")
}

# --- Paths ---
dir_raw       <- "data/raw"         # original data, never changed
dir_processed <- "data/processed"   # cleaned data 
dir_results   <- "results"          # plots and tables

# --- Exclusion rules ---
# The reason for each rule are in docs/decisions.md

# HMP2: only adults with a known age (Franzosa and Hall already all adults)
hmp2_min_age <- 18

# Hall sample with 100% Streptococcus and Shannon = 0, not plausible
exclude_samples <- c("SKST006_10_G102994")

# HMP2 samples without a single count (row sum 0) -> cannot be normalised
hmp2_empty_samples <- c("CSM7KOOR")

#Samples with fewer than 1 million reads are excluded (Franzosa has no read counts,
# so this only applies to Hall and HMP2)
min_reads <- 1e6

# --- Seed ---
# Same results on every run.
seed <- 42
set.seed(seed)


# --- Color ---
# Okabe-Ito palette, readable with color blindness

okabe_ito <- c(orange = "#E69F00", sky_blue = "#56B4E9", green = "#009E73", yellow = "#F0E442",
               blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999",
               black = "#000000")

col_diag <- c(control = okabe_ito[["blue"]], cd = okabe_ito[["vermillion"]])

col_severity <- unname(okabe_ito[c("green", "yellow", "orange", "vermillion")])
#green = normal ; vermillion = severe

#Order in which ggplot assigns the colors to groups. A plot colored by diagnosis
# takes the first two, as in col_diag
palette_order <- unname(okabe_ito[c("blue", "vermillion", "green", "orange", "sky_blue",
                                    "purple", "yellow", "black")])

options(ggplot2.discrete.colour = palette_order, 
        ggplot2.discrete.fill = palette_order)

# --- Plot theme ---
theme_plot <- theme_minimal(base_size = 13) +
  theme(text = 
          element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = 
          element_line(color = "grey92"),
        legend.position = "bottom",
        strip.text = 
          element_text(face = "bold", colour = "black"),
        plot.title = 
          element_text(face= "bold", size = 14))

# --- Helper Functions for cleaning scripts ---
# harmonize species names, so that they are comparable between cohorts
# example : "[Eubacterium] rectale" -> Eubacterium_rectale
# "Roseburia sp. CAG:309" -> "Roseburia_sp_CAG_309"

clean_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  sub("_+$", "", x)
}

# Last cleaning step for all cohorts:
# writed diagnosis in lower case (cd or control), removes faulty samples and samples
# with too few reads, sets the time so that each person's firts remaining sample has time = 0
# and marks that samples in the new column baseline. 
# this table needs columns samples_id, subject_id, diagnosis, time and reads.
finish_meta <- function(meta) {
  #lower-case diagnosis
  meta$diagnosis <- tolower(meta$diagnosis)
  
  #remove faulty samples and samples with too few reads
  low_reads <- !is.na(meta$reads) & meta$reads < min_reads
  meta <- meta [!meta$sample_id %in% c(exclude_samples, hmp2_empty_samples) & !low_reads, ]
  
  #time since the person's first sample
  meta <- meta %>%
    group_by(subject_id) %>%
    mutate(time = if (all(is.na(time))) time else time - min(time, na.rm = TRUE)) %>%
    ungroup()
  
  #baseline = the earliest sample of each person
  base_ids <- meta %>%
    arrange(subject_id, time, sample_id) %>%
    group_by(subject_id) %>%
    slice(1)%>%
    ungroup() %>%
    pull(sample_id)
  meta$baseline <- meta$sample_id %in% base_ids
  
  meta
}

# save one table of a cohort twice: only the abundance, and merged with the metadata
# keeps only the good samples

save_table <- function(table, meta, cohort, name) {
  table <- table[table$sample_id %in% meta$sample_id, ]
  folder <- file.path(dir_processed, cohort)
  write_csv(table, file.path(folder, paste0(cohort, "_", name, ".csv")))
  write_csv(merge(meta, table, by = "sample_id"),
            file.path(folder, paste0(cohort, "_", name, "_merged.csv")))
}
