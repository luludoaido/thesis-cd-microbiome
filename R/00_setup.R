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

# reads per sample from the original metadata
read_depth <- function() {
  hall <- read.csv(file.path(dir_raw, "cmd/HallAB_2017_metadata.csv"))
  hmp2 <- read.csv(file.path(dir_raw, "hmp2/hmp2_metadata_2018-08-20.csv"), sep = ";", check.names = FALSE)
  # with check.name = FALSE, the names stay exactly as in the file
  hmp2 <- hmp2[hmp2$data_type == "metagenomics", ]
  rbind(data.frame(cohort = "hall", sample_id = hall$X, reads = hall$number_reads),
        data.frame(cohort = "hmp2", sample_id = hmp2$`External ID`, reads = hmp2$reads_filtered))
}

# IDs of all samples below min_reads
low_read_samples <- function(cohort = c("hall", "hmp2")) {
  d <- read_depth()
  d$sample_id[d$cohort %in% cohort & !is.na(d$reads) & d$reads < min_reads]
}

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


