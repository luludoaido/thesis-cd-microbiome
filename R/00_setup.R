# --- 00 Setup ---
# Settings shared by all scripts: paths, exclusion rules, seed, colours, plot theme and helper
# functions. Every script starts with source("R/00_setup.R"). Run from the project folder.

# Results: none, no files are written.

# Methods: docs/methods.md, section 1 (and 3.8, 4.3, 4.6)

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
# HMP2: adults with a known age only (see methods 4.3)
hmp2_min_age <- 18

# implausible Hall sample (see methods 3.8)
exclude_samples <- c("SKST006_10_G102994")

# HMP2 sample without a single count (see methods 4.6)
hmp2_empty_samples <- c("CSM7KOOR")

# minimum sequencing depth (see methods 1.2)
min_reads <- 1e6

# --- Seed (see methods 1.7) ---
seed <- 42
set.seed(seed)


# --- Color ---
# Okabe-Ito palette (see methods 1.9)

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
# and marks that sample in the new column baseline (see methods 1.2, 1.5).
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


# --- Confidence interval of the CD effect ---
# 95 % CI of the coefficient diagnosiscd in lm or lmer: t quantile with the model's degrees of
# freedom (Satterthwaite for lmer, see methods 1.11). Returns c(low, high).
ci_cd <- function(fit) {
  s <- coef(summary(fit))["diagnosiscd", ]
  df <- if ("df" %in% names(s)) s[["df"]] else df.residual(fit)
  s[["Estimate"]] + c(-1, 1) * qt(0.975, df) * s[["Std. Error"]]
}

# --- Prevalence filter ---
# Prevalence filter at genus level; a species stays if its genus stays (see methods 1.3).
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
