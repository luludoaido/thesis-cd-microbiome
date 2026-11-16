# --- 00_Set up ---

# What this script does
# Holds the settings that all scripts share. Every script starts with source("R/00_setup.R"),
# so paths, exclusion rules, the seed, colours and helper functions are the same everywhere.
# The scripts have to run from the project folder (open thesis-cd-microbiome.Rproj in RStudio).

# Decisions
# 1. HMP2 only contains adults with a known age, for comparability with the other cohorts.
# 2. One implausible Hall sample and one empty HMP2 sample are excluded.
# 3. Samples with fewer than 1 million reads are excluded.
# 4. The first remaining sample of each person is the baseline, because samples of the same
#    person are not independent.
# 5. The prevalence filter works at genus level and species are not filtered further.
# 6. The seed is fixed, so that every run gives the same results.
# More detail on the exclusions and their effect is in docs/decisions.md.

library(ggplot2)

if (!file.exists("thesis-cd-microbiome.Rproj")) {
  stop("Working directory is not the project folder (", getwd(), ").\n",
       "Open thesis-cd-microbiome.Rproj in RStudio or run Rscript from the project folder.")
}

# --- Paths ---
dir_raw       <- "data/raw"         
dir_processed <- "data/processed"   
dir_results   <- "results"          

# --- Exclusion rules ---
# HMP2: only adults with a known age 
hmp2_min_age <- 18

# Hall sample with 100 % Streptococcus and Shannon = 0 despite 7.3 million reads, not plausible
exclude_samples <- c("SKST006_10_G102994")

# HMP2 samples without a single count (row sum 0) cannot be normalised
hmp2_empty_samples <- c("CSM7KOOR")

# Samples with very few reads show only the most abundant taxa, which distorts richness and Shannon.
# The limit of 1 million reads is an own decision. Franzosa has no read counts, so this only applies
# to Hall and HMP2.
min_reads <- 1e6

# --- Seed ---
seed <- 42
set.seed(seed)


# --- Color ---
# Okabe-Ito palette, readable with color blindness (Wong 2011)

okabe_ito <- c(orange = "#E69F00", sky_blue = "#56B4E9", green = "#009E73", yellow = "#F0E442",
               blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", grey = "#999999",
               black = "#000000")

col_diag <- c(control = okabe_ito[["blue"]], cd = okabe_ito[["vermillion"]])

# green = normal, vermillion = severe
col_severity <- unname(okabe_ito[c("green", "yellow", "orange", "vermillion")])

# Order in which ggplot assigns the colors to groups. A plot colored by diagnosis
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

# --- Helper functions for the cleaning scripts ---
# Harmonizes species names, so that they are comparable between cohorts
# example: "[Eubacterium] rectale" -> "Eubacterium_rectale"
#          "Roseburia sp. CAG:309" -> "Roseburia_sp_CAG_309"

clean_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  sub("_+$", "", x)
}

# Last cleaning step for all cohorts:
# writes the diagnosis in lower case (cd or control), removes faulty samples and samples
# with too few reads, sets the time so that each person's first sample has time = 0
# and marks that sample in the new column baseline.
# The table needs the columns sample_id, subject_id, diagnosis, time and reads.
finish_meta <- function(meta) {
  meta$diagnosis <- tolower(meta$diagnosis)

  low_reads <- !is.na(meta$reads) & meta$reads < min_reads
  meta <- meta [!meta$sample_id %in% c(exclude_samples, hmp2_empty_samples) & !low_reads, ]

  meta <- meta %>%
    group_by(subject_id) %>%
    mutate(time = if (all(is.na(time))) time else time - min(time, na.rm = TRUE)) %>%
    ungroup()

  base_ids <- meta %>%
    arrange(subject_id, time, sample_id) %>%
    group_by(subject_id) %>%
    slice(1)%>%
    ungroup() %>%
    pull(sample_id)
  meta$baseline <- meta$sample_id %in% base_ids

  meta
}

# Saves one table of a cohort twice: only the abundance, and merged with the metadata.
# Keeps only the samples that passed finish_meta.

save_table <- function(table, meta, cohort, name) {
  table <- table[table$sample_id %in% meta$sample_id, ]
  folder <- file.path(dir_processed, cohort)
  write_csv(table, file.path(folder, paste0(cohort, "_", name, ".csv")))
  write_csv(merge(meta, table, by = "sample_id"),
            file.path(folder, paste0(cohort, "_", name, "_merged.csv")))
}


# --- Prevalence filter ---
# Own decision: the filter always works at genus level. At species level, a species stays if
# its genus passes, with no further filter on the species. This keeps species that are rare
# themselves but belong to a common genus. Species without a genus cannot be filtered and stay.
# X: relative abundances, samples in rows. 

passes_filter <- function(X, cohort, level, thr){
  folder <- file.path(dir_processed, cohort)
  if (level == "genus"){
    genus_prev <- colMeans(X > 0)
    taxon_genus <- colnames(X)
  } else {
    genus <- read_csv(file.path(folder, paste0(cohort, "_genus.csv")), show_col_types = FALSE)
    genus <- as.matrix(genus[match(rownames(X), genus$sample_id), -1])
    genus_prev <- colMeans(genus > 0)
    taxonomy <- file.path(folder, paste0(cohort, "_taxonomy.csv"))
    taxon_genus <- if (file.exists(taxonomy)) {
      lookup <- read_csv(taxonomy, show_col_types = FALSE)
      lookup$genus[match(colnames(X), lookup$species)]
    } else {
      sub("_.*", "", colnames(X)) #Franzosa: the genus is the first part of the name
    }
  }
  passing <- names(genus_prev)[if (thr == 0)
    genus_prev > 0 else genus_prev >= thr]
  colSums(X > 0) > 0 & (is.na(taxon_genus)| taxon_genus %in% passing)
}

# --- Literature ---
# Hillmann B et al. (2018). Evaluating the information content of shallow shotgun metagenomics.
#   mSystems 3(6):e00069-18. doi:10.1128/msystems.00069-18
# Wong B (2011). Points of view: Color blindness. Nat Methods 8(6):441. doi:10.1038/nmeth.1618
