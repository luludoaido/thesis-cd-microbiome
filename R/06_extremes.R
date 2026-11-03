library(tidyverse)
library(readxl)
library(vegan)
source("R/00_setup.R")

# =====================================================================
# Extreme genera and samples per cohort
#   1. the 10 most and the 10 least abundant genera (mean relative abundance)
#   2. the 10 samples with the lowest and the highest Shannon diversity
# Data preparation as in R/06_diversity.R:
#   same exclusions, TSS, prevalence filter (main analysis 10 %), TSS again.
# Output: results/extremes/ (tables as CSV + plots), each row with column "cohort".
# =====================================================================

path   <- dir_processed
level  <- "genus"     # "genus" (main analysis) or "species" (as in R/06_diversity.R)
unit   <- if (level == "genus") "genera" else "species"
outdir <- if (level == "genus") file.path(dir_results, "extremes") else
  file.path(dir_results, "extremes_species")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prev_thr <- 0.10   # 0.10 = main analysis; 0 = no filter (then the rarest genera are almost all single occurrences)
n_top    <- 10

# exclusions as in 06_diversity.R (see R/00_setup.R): faulty samples + too few reads
exclude_samples <- union(exclude_samples, low_read_samples())

# colours col_diag and theme_plot come from R/00_setup.R

# ---- 0. Load data (same as R/06_diversity.R) ---------------------------
# returns: list(abund = sample_id + genera, meta = sample_id, subject_id, diagnosis)
# harmonise species names (same as in R/06_diversity.R)
clean_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  sub("_+$", "", x)
}
species_synonyms <- c(
  Phocaeicola_vulgatus = "Bacteroides_vulgatus", Phocaeicola_dorei = "Bacteroides_dorei",
  Phocaeicola_massiliensis = "Bacteroides_massiliensis", Phocaeicola_coprocola = "Bacteroides_coprocola",
  Phocaeicola_plebeius = "Bacteroides_plebeius", Phocaeicola_coprophilus = "Bacteroides_coprophilus",
  Lachnospira_eligens = "Eubacterium_eligens", Segatella_copri = "Prevotella_copri")
sum_duplicate_cols <- function(abund) {
  abund %>% pivot_longer(-sample_id, names_to = "taxon", values_to = "abundance") %>%
    group_by(sample_id, taxon) %>% summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abundance)
}

load_cohort <- function(cohort) {
  if (cohort == "hall") {
    meta_raw   <- read.csv(file.path(path, "hall/hall_meta.csv"))
    if (level == "genus") {
      hall_genus <- read_csv(file.path(path, "hall/hall_genus_merged.csv"), show_col_types = FALSE)
      first_bac  <- which(colnames(hall_genus) == "Abiotrophia")
      abund <- hall_genus %>% select(sample_id, all_of(first_bac:ncol(hall_genus)))
    } else {
      hall_species <- read_csv(file.path(path, "hall/hall_species_merged.csv"), show_col_types = FALSE)
      first_bac <- which(colnames(hall_species) == "visit_number") + 1
      abund <- hall_species %>% select(sample_id, all_of(first_bac:ncol(hall_species)))
      nn <- clean_species(colnames(abund)[-1])
      colnames(abund)[-1] <- ifelse(nn %in% names(species_synonyms), species_synonyms[nn], nn)
      abund <- sum_duplicate_cols(abund)
    }
    meta  <- meta_raw %>% transmute(sample_id = ID, subject_id = as.character(subject_id), diagnosis = diagnosis)
  }
  if (cohort == "hmp2") {
    meta_raw <- read.csv(file.path(path, "hmp2/hmp2_meta.csv"), check.names = FALSE)
    pat <- if (level == "genus") "g__" else "s__"
    hmp2_tab <- read_csv(file.path(path, paste0("hmp2/hmp2_", level, "_merged.csv")), show_col_types = FALSE)
    tax_cols <- grep(pat, colnames(hmp2_tab), value = TRUE)
    abund <- hmp2_tab %>% select(sample_id, all_of(tax_cols)) %>%
      rename_with(~ clean_species(sub(paste0(".*", pat), "", .x)), all_of(tax_cols))
    meta <- meta_raw %>% transmute(sample_id = `External ID`, subject_id = as.character(`Participant ID`), diagnosis = diagnosis)
  }
  if (cohort == "franzosa") {
    meta_raw <- read.csv(file.path(path, "franzosa/franzosa_meta.csv"))
    franzosa_species <- read_csv(file.path(path, "franzosa/franzosa_features_merged.csv"), show_col_types = FALSE)
    first_bac <- which(colnames(franzosa_species) == "Methanobrevibacter_smithii")
    abund <- franzosa_species %>% select(sample_id, all_of(first_bac:ncol(franzosa_species)))
    if (level == "genus") {
      # species -> genus: everything before the first "_"
      abund <- abund %>%
        pivot_longer(-sample_id, names_to = "species", values_to = "abundance") %>%
        mutate(genus = sub("_.*", "", species)) %>%
        group_by(sample_id, genus) %>% summarise(abundance = sum(abundance), .groups = "drop") %>%
        pivot_wider(names_from = genus, values_from = abundance)
    } else {
      colnames(abund)[-1] <- clean_species(colnames(abund)[-1])
    }
    meta <- meta_raw %>% transmute(sample_id = sample_id, subject_id = as.character(sample_id), diagnosis = Diagnosis)
  }
  meta <- meta %>%
    mutate(diagnosis = factor(tolower(diagnosis), levels = c("control", "cd"))) %>%
    filter(!sample_id %in% exclude_samples)
  list(abund = abund, meta = meta)
}

# ---- 1. TSS + prevalence filter + TSS again ------------------------------
prepare <- function(d) {
  X <- d$abund %>%
    filter(sample_id %in% d$meta$sample_id) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  X <- X[rowSums(X) > 0, , drop = FALSE]          # HMP2: CSM7KOOR has row sum 0
  X <- X / rowSums(X)
  prev <- colMeans(X > 0)
  keep <- if (prev_thr == 0) prev > 0 else prev >= prev_thr
  X <- X[, keep, drop = FALSE]
  X <- X[rowSums(X) > 0, , drop = FALSE]
  X / rowSums(X)
}

cohorts <- c(franzosa = "Franzosa", hmp2 = "HMP2", hall = "Hall")
data <- imap(cohorts, function(label, key) {
  d <- load_cohort(key)
  X <- prepare(d)
  list(X = X, meta = d$meta %>% filter(sample_id %in% rownames(X)))
})

# ---- 2. Most and least abundant genera ------------------------------------
# mean relative abundance over all samples (in %), plus prevalence and mean per diagnosis
genera <- imap_dfr(data, function(d, key) {
  X <- d$X
  diag <- d$meta$diagnosis[match(rownames(X), d$meta$sample_id)]
  tibble(cohort      = cohorts[[key]],
         genus       = colnames(X),
         mean_pct    = colMeans(X) * 100,
         prevalence  = colMeans(X > 0) * 100,
         mean_cd     = colMeans(X[diag == "cd", , drop = FALSE]) * 100,
         mean_control = colMeans(X[diag == "control", , drop = FALSE]) * 100)
})

top_genera <- bind_rows(
  genera %>% group_by(cohort) %>% slice_max(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "most abundant"),
  genera %>% group_by(cohort) %>% slice_min(mean_pct, n = n_top, with_ties = FALSE) %>% mutate(group = "least abundant")
) %>% ungroup() %>% arrange(cohort, group, desc(mean_pct))

write_csv(top_genera, file.path(outdir, paste0(unit, "_top_bottom.csv")))
print(top_genera, n = Inf)

# plot: one panel per cohort, bars = mean per diagnosis, sorted by overall mean
plot_genera <- function(which, file, log_x = FALSE) {
  d <- top_genera %>% filter(group == which) %>%
    pivot_longer(c(mean_cd, mean_control), names_to = "diagnosis", values_to = "pct") %>%
    mutate(diagnosis = factor(sub("mean_", "", diagnosis), levels = c("control", "cd")),
           cohort = factor(cohort, levels = cohorts),
           genus_lab = tidytext_reorder(genus, mean_pct, cohort))
  p <- ggplot(d, aes(pct, genus_lab, fill = diagnosis)) +
    geom_col(position = position_dodge(width = .75), width = .7) +
    facet_wrap(~ cohort, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_fill_manual(values = col_diag) +
    labs(x = "mean relative abundance (%)", y = NULL, fill = NULL,
         title = sprintf("%d %s %s per cohort (%s prevalence filter)", n_top, which, unit,
                         if (prev_thr == 0) "no" else paste0(prev_thr * 100, " %"))) +
    theme_plot
  if (log_x) p <- p + scale_x_log10()
  ggsave(file.path(outdir, file), p, width = 12, height = 5, dpi = 200)
  p
}
# helper: sort genera separately within each cohort (like tidytext::reorder_within)
tidytext_reorder <- function(x, by, within) {
  new_x <- paste(x, within, sep = "___")
  reorder(new_x, by)
}
plot_genera("most abundant",  paste0(unit, "_most_abundant.png"))
plot_genera("least abundant", paste0(unit, "_least_abundant.png"), log_x = TRUE)

# ---- 3. Samples with lowest / highest Shannon diversity -------------------
samples <- imap_dfr(data, function(d, key) {
  X <- d$X
  tibble(cohort    = cohorts[[key]],
         sample_id = rownames(X),
         Shannon   = diversity(X, index = "shannon"),
         Richness  = rowSums(X > 0)) %>%
    left_join(d$meta, by = "sample_id")
})

extremes <- bind_rows(
  samples %>% group_by(cohort) %>% slice_min(Shannon, n = n_top, with_ties = FALSE) %>% mutate(group = "lowest diversity"),
  samples %>% group_by(cohort) %>% slice_max(Shannon, n = n_top, with_ties = FALSE) %>% mutate(group = "highest diversity")
) %>% ungroup() %>% arrange(cohort, group, Shannon)

write_csv(extremes, file.path(outdir, "samples_lowest_highest_shannon.csv"))
print(extremes %>% select(cohort, group, sample_id, subject_id, diagnosis, Shannon, Richness), n = Inf)

# how many different persons are behind the extremes? (longitudinal: often the same person)
extremes %>% count(cohort, group, name = "samples") %>%
  left_join(extremes %>% group_by(cohort, group) %>% summarise(people = n_distinct(subject_id), .groups = "drop"),
            by = c("cohort", "group")) %>% print()

# plot: point = sample, label = sample (person), rows = lowest / highest, columns = cohort
d <- extremes %>%
  mutate(cohort = factor(cohort, levels = cohorts),
         group  = factor(group, levels = c("lowest diversity", "highest diversity")),
         label  = ifelse(sample_id == subject_id, sample_id, paste0(sample_id, " (", subject_id, ")")),
         label  = tidytext_reorder(label, Shannon, paste(cohort, group)))
p <- ggplot(d, aes(Shannon, label, colour = diagnosis)) +
  geom_segment(aes(x = 0, xend = Shannon, yend = label), colour = "grey85", linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(group ~ cohort, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_colour_manual(values = col_diag) +
  labs(x = "Shannon diversity", y = NULL, colour = NULL,
       title = sprintf("%d samples with the lowest and highest diversity per cohort", n_top)) +
  theme_plot + theme(axis.text.y = element_text(size = 8))
ggsave(file.path(outdir, "samples_lowest_highest_shannon.png"), p, width = 14, height = 8, dpi = 200)

# the same plots separately (lowest and highest), with larger labels
for (g in levels(d$group)) {
  p <- ggplot(filter(d, group == g), aes(Shannon, label, colour = diagnosis)) +
    geom_segment(aes(x = 0, xend = Shannon, yend = label), colour = "grey85", linewidth = 1) +
    geom_point(size = 3) +
    facet_wrap(~ cohort, scales = "free_y", ncol = 3) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(limits = c(0, ceiling(max(d$Shannon)))) +
    scale_colour_manual(values = col_diag) +
    labs(x = "Shannon diversity", y = NULL, colour = NULL,
         title = sprintf("%d samples with the %s per cohort", n_top, g)) +
    theme_plot + theme(axis.text.y = element_text(size = 10))
  ggsave(file.path(outdir, paste0("samples_", sub(" diversity", "", g), "_shannon.png")), p, width = 14, height = 5, dpi = 200)
}
