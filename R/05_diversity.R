library(tidyverse)
library(readxl)
library(vegan)
library(uwot)      # UMAP
library(lme4)      # Mixed Models
library(lmerTest)  # p-values for mixed models
library(limma)     # Bioconductor: removeBatchEffect
library(sva)       # Bioconductor: ComBat
library(MMUPHin)   # Bioconductor: batch correction for microbiome data
source("R/00_setup.R")

# =====================================================================
# Abundances: sparsity, few genera, alpha diversity, PCA/UMAP
# Workflow: choose the cohort at the top, then run the whole script.
# For the next cohort only change `cohort` and run again.
# Prevalence filter with 0 / 5 / 10 % threshold; longitudinal for Hall + HMP2.
# `level` selects the taxonomic level: "genus" (main analysis)
# or "species". Variables and plots are still called "genera" / "n_genera" –
# with level = "species" they refer to species.
# =====================================================================

cohort <- "hall"      # "hall", "hmp2" or "franzosa"
level  <- "genus"     # "genus" (main analysis) or "species"
# can be overridden in the terminal: Rscript R/06_diversity.R hmp2 genus
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) cohort <- args[1]
if (length(args) >= 2) level  <- args[2]

# sensitivity analysis: exclude persons, e.g. c("C3009") for HMP2
# (results then go to a separate folder, e.g. results/hmp2_without_C3009)
exclude_subjects <- character(0)
# or all persons whose ID matches a pattern, e.g. "^SKST" for Hall
# (SKST = Atlanta + different DNA extraction kit; all CD persons are Boston + Qiagen)
exclude_pattern <- NA_character_

# exclusions that always apply (exclude_samples, min_reads) are in R/00_setup.R
raw_path  <- dir_raw

path   <- dir_processed
excl_label <- c(exclude_subjects, if (!is.na(exclude_pattern)) gsub("[^A-Za-z0-9]", "", exclude_pattern))
run_name <- if (length(excl_label) == 0) cohort else
  paste0(cohort, "_without_", paste(excl_label, collapse = "_"))
# genus: results/<run_name>; species: results/species/<run_name>
outdir <- if (level == "genus") file.path(dir_results, run_name) else
  file.path(file.path(dir_results, "species"), run_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# harmonise species names so that they are comparable between cohorts:
# "[Eubacterium] rectale" -> "Eubacterium_rectale", "Roseburia sp. CAG:309" -> "Roseburia_sp_CAG_309"
clean_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  sub("_+$", "", x)
}
# Hall uses a newer taxonomy -> map back to the MetaPhlAn2 names of HMP2/Franzosa
# (only the cases relevant for the hypotheses or frequent ones; extend if needed)
species_synonyms <- c(
  Phocaeicola_vulgatus    = "Bacteroides_vulgatus",
  Phocaeicola_dorei       = "Bacteroides_dorei",
  Phocaeicola_massiliensis = "Bacteroides_massiliensis",
  Phocaeicola_coprocola   = "Bacteroides_coprocola",
  Phocaeicola_plebeius    = "Bacteroides_plebeius",
  Phocaeicola_coprophilus = "Bacteroides_coprophilus",
  Lachnospira_eligens     = "Eubacterium_eligens",
  Segatella_copri         = "Prevotella_copri"
)
# sum columns with the same name (after cleaning/mapping)
sum_duplicate_cols <- function(abund) {
  abund %>%
    pivot_longer(-sample_id, names_to = "taxon", values_to = "abundance") %>%
    group_by(sample_id, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abundance)
}

# ---- 0. Data loading --------------------------------------------------
# goal of this block: two objects
#   abund: table, first column sample_id, rest = genera
#   meta : table with sample_id, subject_id, diagnosis, time (weeks since first sample),
#          age, antibiotic (NA if the cohort has no data)

if (cohort == "hall") {
  meta_raw   <- read.csv(file.path(path, "hall/hall_meta.csv"))
  if (level == "genus") {
    hall_genus <- read_csv(file.path(path, "hall/hall_genus_merged.csv"))
    first_bac <- which(colnames(hall_genus) == "Abiotrophia")
    abund <- hall_genus %>% select(sample_id, all_of(first_bac:ncol(hall_genus)))
  } else {
    # species table: the first columns are metadata, bacteria start at "visit_number" + 1
    hall_species <- read_csv(file.path(path, "hall/hall_species_merged.csv"))
    first_bac <- which(colnames(hall_species) == "visit_number") + 1
    abund <- hall_species %>% select(sample_id, all_of(first_bac:ncol(hall_species)))
    new_names <- clean_species(colnames(abund)[-1])
    new_names <- ifelse(new_names %in% names(species_synonyms), species_synonyms[new_names], new_names)
    colnames(abund)[-1] <- new_names
    abund <- sum_duplicate_cols(abund)
  }

  meta <- meta_raw %>%
    transmute(sample_id = ID, subject_id = as.character(subject_id),
              diagnosis = diagnosis,
              time = days_from_first_collection / 7,
              age = NA_real_, antibiotic = NA_character_,
              # centre: SKST = Atlanta (different DNA kit), all others = Boston (Qiagen)
              site = ifelse(grepl("^SKST", subject_id), "Atlanta", "Boston"))

  # sequencing depth from the original metadata
  reads <- read.csv(file.path(raw_path, "cmd/HallAB_2017_metadata.csv")) %>%
    transmute(sample_id = X, reads = number_reads)
  meta <- meta %>% left_join(reads, by = "sample_id")
}

if (cohort == "hmp2") {
  meta_raw <- read.csv(file.path(path, "hmp2/hmp2_meta.csv"), check.names = FALSE)
  if (level == "genus") {
    hmp2_genus <- read_csv(file.path(path, "hmp2/hmp2_genus_merged.csv"))
    # genus columns are called "k__...|g__Bacteroides" -> keep only the genus name
    genus_cols <- grep("g__", colnames(hmp2_genus), value = TRUE)
    abund <- hmp2_genus %>%
      select(sample_id, all_of(genus_cols)) %>%
      rename_with(~ sub(".*g__", "", .x), all_of(genus_cols))
  } else {
    hmp2_species <- read_csv(file.path(path, "hmp2/hmp2_species_merged.csv"))
    # species columns are called "k__...|g__...|s__Faecalibacterium_prausnitzii" -> keep only the species name
    species_cols <- grep("s__", colnames(hmp2_species), value = TRUE)
    abund <- hmp2_species %>%
      select(sample_id, all_of(species_cols)) %>%
      rename_with(~ clean_species(sub(".*s__", "", .x)), all_of(species_cols))
  }

  meta <- meta_raw %>%
    transmute(sample_id = `External ID`,
              subject_id = as.character(`Participant ID`),
              diagnosis = diagnosis,
              time = week_num,
              age = consent_age, antibiotic = NA_character_,
              # centre = first letter of the Participant ID (C = Cedars-Sinai, M = MGH, ...)
              site = substr(`Participant ID`, 1, 1))

  # sequencing depth from the original metadata (metagenome rows only)
  reads <- read.csv(file.path(raw_path, "hmp2/hmp2_metadata_2018-08-20.csv"), sep = ";", check.names = FALSE) %>%
    filter(data_type == "metagenomics") %>%
    transmute(sample_id = `External ID`, reads = reads_filtered)
  meta <- meta %>% left_join(reads, by = "sample_id")
}

if (cohort == "franzosa") {
  meta_raw <- read.csv(file.path(path, "franzosa/franzosa_meta.csv"))
  franzosa_species <- read_csv(file.path(path, "franzosa/franzosa_features_merged.csv"))

  # Franzosa is only available at species level -> sum to genus
  # (genus = everything before the first "_", e.g. "Bacteroides_dorei" -> "Bacteroides";
  #  "Bacteroides_unclassified" thus also goes to Bacteroides)
  # with level = "species" the table is used directly (incl. "..._unclassified")
  first_bac <- which(colnames(franzosa_species) == "Methanobrevibacter_smithii")
  abund <- franzosa_species %>%
    select(sample_id, all_of(first_bac:ncol(franzosa_species)))
  if (level == "genus") {
    abund <- abund %>%
      pivot_longer(-sample_id, names_to = "species", values_to = "abundance") %>%
      mutate(genus = sub("_.*", "", species)) %>%
      group_by(sample_id, genus) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      pivot_wider(names_from = genus, values_from = abundance)
  } else {
    colnames(abund)[-1] <- clean_species(colnames(abund)[-1])
  }

  meta <- meta_raw %>%
    transmute(sample_id = sample_id, subject_id = as.character(sample_id),   # 1 sample per person
              diagnosis = Diagnosis, time = NA_real_,
              age = Age, antibiotic = antibiotic,
              # sub-cohort: PRISM or Validation (before the "|" in sample_id)
              site = sub("[|].*", "", sample_id),
              reads = NA_real_)   # no read counts available
}

# harmonise diagnosis labels (Franzosa: "Control")
meta <- meta %>% mutate(diagnosis = factor(tolower(diagnosis), levels = c("control", "cd")))

# remove excluded persons (see exclude_subjects above)
if (!is.na(exclude_pattern)) {
  exclude_subjects <- union(exclude_subjects, unique(grep(exclude_pattern, meta$subject_id, value = TRUE)))
}
if (length(exclude_subjects) > 0) {
  message("Excluding: ", paste(exclude_subjects, collapse = ", "), " (",
          sum(meta$subject_id %in% exclude_subjects), " Samples)")
  meta <- meta %>% filter(!subject_id %in% exclude_subjects)
}
if (any(meta$sample_id %in% exclude_samples)) {
  message("Excluding technically faulty samples: ",
          paste(intersect(meta$sample_id, exclude_samples), collapse = ", "))
  meta <- meta %>% filter(!sample_id %in% exclude_samples)
}
# too low sequencing depth (see min_reads in R/00_setup.R); NA = no value -> kept
low_reads <- meta %>% filter(!is.na(reads), reads < min_reads)
if (nrow(low_reads) > 0) {
  message("Excluding samples with < ", format(min_reads, big.mark = "'", scientific = FALSE), " reads: ",
          paste0(low_reads$sample_id, " (", low_reads$diagnosis, ", ", low_reads$reads, ")", collapse = ", "))
  meta <- meta %>% filter(!sample_id %in% low_reads$sample_id)
}
table(meta$diagnosis, useNA = "ifany")

# time per person since the first sample (0 = first sample)
# (depending on the person, Hall starts at 30-300 days, HMP2 at week 0-8)
meta <- meta %>%
  group_by(subject_id) %>%
  mutate(time = if (all(is.na(time))) time else time - min(time, na.rm = TRUE)) %>%
  ungroup()

# abundance matrix matching the metadata
length(intersect(meta$sample_id, abund$sample_id))
setdiff(meta$sample_id, abund$sample_id)     # only in metadata
setdiff(abund$sample_id, meta$sample_id)     # only in abundances

X_raw <- abund %>%
  filter(sample_id %in% meta$sample_id) %>%
  arrange(match(sample_id, meta$sample_id)) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# samples without a single count cannot be normalised (HMP2: CSM7KOOR)
rownames(X_raw)[rowSums(X_raw) == 0]
X_raw <- X_raw[rowSums(X_raw) > 0, , drop = FALSE]

meta <- meta %>% filter(sample_id %in% rownames(X_raw))
dim(X_raw); nrow(meta)

longitudinal <- max(table(meta$subject_id)) > 1
longitudinal

# baseline = earliest sample per person (for cross-sectional)
if (longitudinal && any(is.na(meta$time))) warning("time is missing for some samples -> baseline may be wrong!")
base_ids <- meta %>%
  arrange(subject_id, time, sample_id) %>%
  group_by(subject_id) %>% slice(1) %>% ungroup() %>%
  pull(sample_id)
meta$baseline <- meta$sample_id %in% base_ids

# ---- 1. Row sums + TSS ------------------------------------------------
dim(X_raw)
rs <- rowSums(X_raw)
summary(rs)            # 1 or 100?
hist(rs, breaks = 30)
mean(X_raw == 0)       # sparsity without filter

X_tss <- X_raw / rowSums(X_raw)

# ---- 2. Sparsity with thresholds 0 / 5 / 10 % -------------------------
prev <- colMeans(X_tss > 0)
hist(prev, breaks = 30, main = "Prevalence of genera", xlab = "share of samples")

thresholds <- c(0, 0.05, 0.10)

# thr = 0: only remove genera that occur in no sample
X_thr <- lapply(thresholds, function(t) {
  keep <- if (t == 0) prev > 0 else prev >= t
  Xf <- X_tss[, keep, drop = FALSE]
  Xf <- Xf[rowSums(Xf) > 0, , drop = FALSE]      # samples without any genus after the filter
  Xf / rowSums(Xf)                                # TSS again
})
names(X_thr) <- sprintf("thr_%02d", thresholds * 100)   # thr_00, thr_05, thr_10

# main analysis: 10 % filter (fewer zeros -> CLR depends less on the pseudocount);
# 0 % and 5 % run as sensitivity analysis. All single plots use main_thr.
main_thr <- "thr_10"

sparsity_tab <- tibble(
  threshold    = thresholds,
  n_samples    = sapply(X_thr, nrow),
  n_genera     = sapply(X_thr, ncol),
  sparsity     = sapply(X_thr, function(x) mean(x == 0)),
  samples_lost = nrow(X_tss) - sapply(X_thr, nrow)
)
sparsity_tab

# ---- 3. Samples with very few genera ----------------------------------
# number of genera per sample per threshold
few <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X), n_genera = rowSums(X > 0), threshold = nm)
}) %>%
  left_join(meta, by = "sample_id") %>%
  group_by(threshold) %>%
  mutate(few_genera = n_genera <= quantile(n_genera, 0.05)) %>%     # lowest 5 %, adjust cut-off if needed
  ungroup()

few %>% group_by(threshold, diagnosis) %>%
  summarise(median = median(n_genera), min = min(n_genera), .groups = "drop")

p_few <- ggplot(few, aes(n_genera, fill = diagnosis)) +
  geom_histogram(bins = 30, alpha = .7, position = "identity") +
  facet_wrap(~ threshold, scales = "free")
p_few
ggsave(file.path(outdir, "few_genera.png"), p_few, width = 9, height = 4)

# the samples themselves (incl. row sum before TSS, to spot technical outliers)
few %>% filter(few_genera) %>%
  mutate(total_raw = rs[sample_id]) %>%
  arrange(threshold, n_genera) %>%
  select(threshold, sample_id, subject_id, diagnosis, n_genera, total_raw)

# longitudinal: do the same persons appear several times?
few %>% filter(few_genera) %>%
  count(threshold, subject_id, diagnosis, name = "n_few_samples") %>%
  arrange(threshold, desc(n_few_samples))

# ---- 4. Alpha diversity -----------------------------------------------
alpha <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X),
         Richness = rowSums(X > 0),
         Shannon  = diversity(X, "shannon"),
         Simpson  = diversity(X, "simpson"),
         threshold = nm)
}) %>% left_join(meta, by = "sample_id")

# 4a. cross-sectional: baseline only, 1 sample per person
# (Wilcoxon assumes independent observations -> do not use all samples)
alpha %>% filter(baseline) %>%
  group_by(threshold, diagnosis) %>%
  summarise(n = n(), median_shannon = median(Shannon), median_richness = median(Richness),
            .groups = "drop")

alpha %>% filter(baseline) %>%
  group_by(threshold) %>%
  summarise(p_shannon  = wilcox.test(Shannon ~ diagnosis)$p.value,
            p_richness = wilcox.test(Richness ~ diagnosis)$p.value)

p_alpha <- alpha %>% filter(baseline) %>%
  ggplot(aes(diagnosis, Shannon, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .15, size = .8) +
  facet_wrap(~ threshold) + labs(title = "Shannon, baseline samples")
p_alpha
ggsave(file.path(outdir, "alpha_baseline.png"), p_alpha, width = 9, height = 4)

# 4a'. median per person: also 1 value per person, but from all samples
# (more robust than baseline, a single bad sample hardly counts)
if (longitudinal) {
  alpha_person <- alpha %>%
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_samples = n(), Shannon = median(Shannon), Richness = median(Richness),
              .groups = "drop")

  alpha_person %>%
    group_by(threshold, diagnosis) %>%
    summarise(n = n(), median_shannon = median(Shannon), median_richness = median(Richness),
              .groups = "drop")

  median_tab <- alpha_person %>%
    group_by(threshold) %>%
    summarise(p_shannon  = wilcox.test(Shannon ~ diagnosis)$p.value,
              p_richness = wilcox.test(Richness ~ diagnosis)$p.value)
  print(median_tab)

  p_alpha_med <- ggplot(alpha_person, aes(diagnosis, Shannon, fill = diagnosis)) +
    geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .15, size = .8) +
    facet_wrap(~ threshold) + labs(title = "Shannon, median per person")
  print(p_alpha_med)
  ggsave(file.path(outdir, "alpha_median_person.png"), p_alpha_med, width = 9, height = 4)
}

# 4b. longitudinal: all samples, person as random effect
if (longitudinal) {
  mixed_tab <- alpha %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      m <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d)
      vc <- as.data.frame(VarCorr(m))
      tibble(estimate_cd = fixef(m)[2],                  # difference CD - control
             p_mixed = coef(summary(m))[2, "Pr(>|t|)"],
             icc = vc$vcov[1] / sum(vc$vcov))            # share of variance between persons
    })
  print(mixed_tab)

  # spaghetti plot
  p_spag <- alpha %>% filter(threshold == main_thr) %>%
    arrange(subject_id, time) %>%
    ggplot(aes(time, Shannon, group = subject_id, colour = diagnosis)) +
    geom_line(alpha = .4) + geom_point(size = .8, alpha = .6) +
    labs(x = "weeks since first sample")
  print(p_spag)
  ggsave(file.path(outdir, "alpha_spaghetti.png"), p_spag, width = 8, height = 5)
}

# 4c. within a person: how stable is the microbiome over time?
if (longitudinal) {
  # (1) variation of alpha diversity: SD of Shannon per person (only persons with >= 3 samples)
  alpha_var <- alpha %>%
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_samples = n(), sd_shannon = sd(Shannon),
              range_shannon = max(Shannon) - min(Shannon), .groups = "drop") %>%
    filter(n_samples >= 3)

  alpha_var %>% group_by(threshold, diagnosis) %>%
    summarise(persons = n(), median_sd = median(sd_shannon), .groups = "drop")

  alpha_var_test <- alpha_var %>% group_by(threshold) %>%
    summarise(p_sd_shannon = wilcox.test(sd_shannon ~ diagnosis)$p.value)
  print(alpha_var_test)

  # (2) change in composition: Bray-Curtis distance
  #     0 = identical composition, 1 = no shared genera
  beta_pairs <- imap_dfr(X_thr, function(X, nm) {
    D <- as.matrix(vegdist(X, method = "bray"))
    m <- meta[match(rownames(X), meta$sample_id), ]
    idx <- which(upper.tri(D), arr.ind = TRUE)
    tibble(threshold = nm,
           s1 = rownames(X)[idx[, 1]], s2 = rownames(X)[idx[, 2]],
           p1 = m$subject_id[idx[, 1]], p2 = m$subject_id[idx[, 2]],
           d1 = m$diagnosis[idx[, 1]],  d2 = m$diagnosis[idx[, 2]],
           t1 = m$time[idx[, 1]],       t2 = m$time[idx[, 2]],
           bray = D[idx])
  })

  # within the same person vs. between different persons (same diagnosis)
  within_between <- beta_pairs %>%
    filter(d1 == d2) %>%
    mutate(diagnosis = d1,
           comparison = if_else(p1 == p2, "same person", "other person"))

  within_between %>% group_by(threshold, diagnosis, comparison) %>%
    summarise(pairs = n(), median_bray = median(bray), .groups = "drop")

  p_within <- within_between %>% filter(threshold == main_thr) %>%
    ggplot(aes(comparison, bray, fill = diagnosis)) +
    geom_boxplot(outlier.size = .3) +
    labs(x = NULL, y = "Bray-Curtis distance", title = "Samples of the same vs. different persons")
  print(p_within)
  ggsave(file.path(outdir, "beta_within_between.png"), p_within, width = 7, height = 5)

  # stability per person: distance between consecutive samples
  consecutive <- meta %>%
    filter(!is.na(time)) %>%
    arrange(subject_id, time) %>%
    group_by(subject_id) %>%
    mutate(next_id = lead(sample_id), gap_weeks = lead(time) - time) %>%
    ungroup() %>%
    filter(!is.na(next_id)) %>%
    select(subject_id, diagnosis, s_a = sample_id, s_b = next_id, gap_weeks)

  consecutive <- beta_pairs %>%
    select(threshold, s1, s2, bray) %>%
    bind_rows(rename(., s1 = s2, s2 = s1)) %>%           # pairs in both directions
    inner_join(consecutive, by = c("s1" = "s_a", "s2" = "s_b"))

  stability <- consecutive %>%
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_pairs = n(), median_bray_consec = median(bray), .groups = "drop")

  stability %>% group_by(threshold, diagnosis) %>%
    summarise(persons = n(), median = median(median_bray_consec), .groups = "drop")

  stability_test <- stability %>% group_by(threshold) %>%
    summarise(p_stability = wilcox.test(median_bray_consec ~ diagnosis)$p.value)
  print(stability_test)

  # does the change grow with the time between two samples?
  p_gap <- consecutive %>% filter(threshold == main_thr) %>%
    ggplot(aes(gap_weeks, bray, colour = diagnosis)) +
    geom_jitter(width = .2, alpha = .5) + geom_smooth(method = "lm", se = TRUE) +   # lm: Hall has only a few distinct gaps, loess fails there
    labs(x = "weeks between two consecutive samples", y = "Bray-Curtis distance")
  print(p_gap)
  ggsave(file.path(outdir, "beta_consecutive_gap.png"), p_gap, width = 7, height = 5)
}

# 4d. time course: does Shannon change over time, and differently in CD than in controls?
# time = weeks since the person's first sample
# m_null: diagnosis only; m_time: diagnosis * time (interaction = different slope)
# compared by likelihood-ratio test (REML = FALSE so that the models are comparable)
# m_slope: additionally a separate slope per person, (1 + time | subject_id); only converges with
# the optimiser "bobyqa" (the default optimiser gives a convergence warning)
if (longitudinal) {
  time_tab <- alpha %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      m_null <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d, REML = FALSE)
      m_time <- lmer(Shannon ~ diagnosis * time + (1 | subject_id), data = d, REML = FALSE)
      cf <- coef(summary(m_time))
      tibble(cd_at_start   = cf["diagnosiscd", "Estimate"],      # difference at the first sample
             p_cd_at_start = cf["diagnosiscd", "Pr(>|t|)"],
             slope_control = cf["time", "Estimate"],         # Shannon per week, control
             p_slope_control = cf["time", "Pr(>|t|)"],
             slope_cd = cf["time", "Estimate"] + cf["diagnosiscd:time", "Estimate"],
             p_interaction = cf["diagnosiscd:time", "Pr(>|t|)"],  # do the slopes differ?
             p_lrt = anova(m_null, m_time)$`Pr(>Chisq)`[2])     # does time improve the model?
    })
  print(time_tab)

  time_slope_tab <- alpha %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      m_slope <- tryCatch(
        lmer(Shannon ~ diagnosis * time + (1 + time | subject_id), data = d, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
        warning = function(w) NULL)
      if (is.null(m_slope) || isSingular(m_slope)) return(tibble(converged = FALSE))
      cf <- coef(summary(m_slope))
      tibble(converged = TRUE,
             slope_control = cf["time", "Estimate"], p_slope_control = cf["time", "Pr(>|t|)"],
             slope_cd = cf["time", "Estimate"] + cf["diagnosiscd:time", "Estimate"],
             p_interaction = cf["diagnosiscd:time", "Pr(>|t|)"])
    })
  print(time_slope_tab)

  # how many persons increase? (slope per person, only persons with >= 5 samples)
  person_slopes <- alpha %>% filter(threshold == main_thr) %>%
    group_by(subject_id, diagnosis) %>% filter(n() >= 5) %>%
    summarise(n = n(), slope = coef(lm(Shannon ~ time))[2], .groups = "drop")
  person_slopes %>% group_by(diagnosis) %>%
    summarise(persons = n(), increasing = sum(slope > 0), median_slope = median(slope))

  p_trend <- alpha %>% filter(threshold == main_thr) %>%
    ggplot(aes(time, Shannon, colour = diagnosis)) +
    geom_point(alpha = .25, size = .8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = "weeks since first sample", title = "Shannon over time (linear trends per diagnosis)")
  print(p_trend)
  ggsave(file.path(outdir, "alpha_time_trend.png"), p_trend, width = 8, height = 5)
}

# 4e. age as covariate: does the CD effect change when adjusting for age?
# (same persons in both models: only samples with age)
if (any(!is.na(alpha$age))) {
  age_tab <- alpha %>%
    filter(!is.na(age)) %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      if (longitudinal) {
        m0 <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d)
        m1 <- lmer(Shannon ~ diagnosis + age + (1 | subject_id), data = d)
      } else {
        m0 <- lm(Shannon ~ diagnosis, data = d)
        m1 <- lm(Shannon ~ diagnosis + age, data = d)
      }
      c0 <- coef(summary(m0)); c1 <- coef(summary(m1))
      pcol <- grep("^Pr", colnames(c1), value = TRUE)
      tibble(persons = n_distinct(d$subject_id),
             cd_without_age = c0["diagnosiscd", "Estimate"], p_without_age = c0["diagnosiscd", pcol],
             cd_with_age  = c1["diagnosiscd", "Estimate"], p_with_age  = c1["diagnosiscd", pcol],
             age_per_year = c1["age", "Estimate"],        p_age      = c1["age", pcol])
    })
  print(age_tab)
}

# 4f. sensitivity: without persons on antibiotics (only Franzosa has this information)
# antibiotics lower diversity independently of CD; in Franzosa only CD persons take them
if (any(alpha$antibiotic == "Yes", na.rm = TRUE)) {
  alpha %>% filter(threshold == main_thr, diagnosis == "cd") %>%
    group_by(antibiotic) %>%
    summarise(n = n(), median_shannon = median(Shannon), median_richness = median(Richness))

  abx_tab <- alpha %>%
    filter(baseline) %>%
    group_by(threshold) %>%
    summarise(p_shannon_all = wilcox.test(Shannon ~ diagnosis)$p.value,
              p_richness_all = wilcox.test(Richness ~ diagnosis)$p.value,
              n_without_abx = sum(antibiotic != "Yes" | is.na(antibiotic)),
              p_shannon_without_abx = wilcox.test(Shannon[antibiotic != "Yes" | is.na(antibiotic)] ~
                                                diagnosis[antibiotic != "Yes" | is.na(antibiotic)])$p.value,
              p_richness_without_abx = wilcox.test(Richness[antibiotic != "Yes" | is.na(antibiotic)] ~
                                                 diagnosis[antibiotic != "Yes" | is.na(antibiotic)])$p.value)
  print(abx_tab)
}

# 4g. centre or sub-cohort: does the CD effect depend on where the samples come from?
# HMP2: 13 of 14 controls from centre M, 19 of 31 CD from centre C
# Hall: all CD from Boston, 10 of 12 controls from Atlanta (completely confounded with centre)
# Franzosa: PRISM and Validation both contain CD and controls
# (1) CD effect with centre as covariate, (2) CD effect only in the centre with the most controls
site_tab <- alpha %>%
  group_by(threshold) %>%
  group_modify(function(d, key) {
    ok <- d %>% distinct(subject_id, diagnosis, site) %>% count(site, diagnosis) %>%
      group_by(site) %>% filter(n_distinct(diagnosis) == 2, min(n) >= 3) %>% pull(site) %>% unique()
    if (length(unique(d$site)) < 2 || length(ok) == 0) return(tibble(note = "centre cannot be separated from diagnosis"))
    main_site <- d %>% filter(site %in% ok, diagnosis == "control") %>% count(site) %>% slice_max(n, n = 1) %>% pull(site)
    fit <- function(f, dd) {
      m <- if (longitudinal) lmer(update(f, . ~ . + (1 | subject_id)), data = dd) else lm(f, data = dd)
      cf <- coef(summary(m)); pcol <- grep("^Pr", colnames(cf), value = TRUE)
      c(cf["diagnosiscd", "Estimate"], cf["diagnosiscd", pcol])
    }
    a0 <- fit(Shannon ~ diagnosis, d); a1 <- fit(Shannon ~ diagnosis + site, d)
    a2 <- fit(Shannon ~ diagnosis, filter(d, site == main_site))
    r0 <- fit(Richness ~ diagnosis, d); r1 <- fit(Richness ~ diagnosis + site, d)
    tibble(site_only = main_site,
           shannon_cd = a0[1], p_shannon = a0[2],
           shannon_cd_with_site = a1[1], p_shannon_with_site = a1[2],
           shannon_cd_site_only = a2[1], p_shannon_site_only = a2[2],
           richness_cd = r0[1], p_richness = r0[2],
           richness_cd_with_site = r1[1], p_richness_with_site = r1[2])
  })
print(site_tab)

# alpha diversity of the few-genera samples in comparison
alpha %>%
  left_join(select(few, sample_id, threshold, few_genera), by = c("sample_id", "threshold")) %>%
  group_by(threshold, few_genera) %>%
  summarise(n = n(), median_shannon = median(Shannon), .groups = "drop")

# ---- 5. PCA + UMAP at 0 / 5 / 10 % -----------------------------------
# CLR with pseudocount (half the smallest non-zero value)
clr_tf <- function(X) {
  pc <- min(X[X > 0]) / 2
  L <- log(X + pc)
  L - rowMeans(L)
}

ord <- imap_dfr(X_thr, function(X, nm) {
  Z <- clr_tf(X)
  pca <- prcomp(Z, center = TRUE, scale. = FALSE)
  set.seed(seed)
  um <- umap(Z, n_neighbors = min(15, nrow(Z) - 1), min_dist = 0.3)

  sc <- tibble(sample_id = rownames(Z), threshold = nm,
               PC1 = pca$x[, 1], PC2 = pca$x[, 2],
               UMAP1 = um[, 1], UMAP2 = um[, 2],
               var1 = 100 * summary(pca)$importance[2, 1],
               var2 = 100 * summary(pca)$importance[2, 2])
  # outliers in the PCA: distance to the median > median + 3 MAD
  d <- sqrt((sc$PC1 - median(sc$PC1))^2 + (sc$PC2 - median(sc$PC2))^2)
  sc$outlier <- d > median(d) + 3 * mad(d)
  sc
}) %>%
  left_join(meta, by = "sample_id") %>%
  left_join(select(few, sample_id, threshold, n_genera, few_genera),
            by = c("sample_id", "threshold"))

# explained variance per threshold (for the axis labels)
ord %>% distinct(threshold, var1, var2)

# PCA: colour = diagnosis, shape = outlier
p_pca <- ggplot(ord, aes(PC1, PC2, colour = diagnosis, shape = outlier)) +
  geom_point(alpha = .7) + facet_wrap(~ threshold, scales = "free")
p_pca
ggsave(file.path(outdir, "pca_thresholds.png"), p_pca, width = 11, height = 4)

# PCA: colour = number of genera (do few-genera samples lie at the edge?)
p_pca_n <- ggplot(ord, aes(PC1, PC2, colour = n_genera)) +
  geom_point(alpha = .8) + scale_colour_viridis_c() + facet_wrap(~ threshold, scales = "free")
p_pca_n
ggsave(file.path(outdir, "pca_n_genera.png"), p_pca_n, width = 11, height = 4)

# UMAP
p_umap <- ggplot(ord, aes(UMAP1, UMAP2, colour = diagnosis, shape = few_genera)) +
  geom_point(alpha = .8) + facet_wrap(~ threshold, scales = "free")
p_umap
ggsave(file.path(outdir, "umap_thresholds.png"), p_umap, width = 11, height = 4)

# longitudinal: connect samples of the same person in time order
if (longitudinal) {
  p_traj <- ord %>% filter(threshold == main_thr) %>%
    arrange(subject_id, time) %>%
    ggplot(aes(PC1, PC2, group = subject_id, colour = diagnosis)) +
    geom_path(alpha = .3) + geom_point(alpha = .7)
  print(p_traj)
  ggsave(file.path(outdir, "pca_subjects.png"), p_traj, width = 7, height = 5)
}

# PCoA with Bray-Curtis (sensitivity to the PCA on CLR)
# Bray-Curtis is not Euclidean -> no PCA but PCoA (classical MDS)
pcoa <- imap_dfr(X_thr, function(X, nm) {
  pc <- cmdscale(vegdist(X, method = "bray"), k = 2, eig = TRUE)
  expl <- 100 * pc$eig / sum(pc$eig[pc$eig > 0])
  tibble(sample_id = rownames(X), threshold = nm,
         PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2],
         var1 = expl[1], var2 = expl[2])
}) %>% left_join(meta, by = "sample_id")

pcoa %>% distinct(threshold, var1, var2)

p_pcoa <- ggplot(pcoa, aes(PCoA1, PCoA2, colour = diagnosis)) +
  geom_point(alpha = .7) + facet_wrap(~ threshold, scales = "free") +
  labs(title = "PCoA, Bray-Curtis")
p_pcoa
ggsave(file.path(outdir, "pcoa_bray_thresholds.png"), p_pcoa, width = 11, height = 4)

# beta diversity: do CD and control differ? (PERMANOVA)
# main analysis Aitchison (= Euclidean on CLR), sensitivity Bray-Curtis
# baseline samples only, so that persons do not count more than once
permanova <- imap_dfr(X_thr, function(X, nm) {
  ids <- intersect(rownames(X), meta$sample_id[meta$baseline])
  d <- meta[match(ids, meta$sample_id), ]
  dists <- list(aitchison = dist(clr_tf(X)[ids, ]),
                bray      = vegdist(X[ids, ], method = "bray"))
  imap_dfr(dists, function(D, dn) {
    set.seed(seed)
    a <- adonis2(D ~ diagnosis, data = d, permutations = 999)
    # dispersion: does the spread of the groups differ? (otherwise PERMANOVA can
    # also become significant only because of different spread)
    set.seed(seed)
    bd <- permutest(betadisper(D, d$diagnosis), permutations = 999)
    tibble(threshold = nm, distance = dn, n = length(ids),
           R2 = a$R2[1], p = a$`Pr(>F)`[1], p_dispersion = bd$tab$`Pr(>F)`[1])
  })
})
permanova %>% arrange(distance, threshold)

# PERMANOVA without persons on antibiotics (Franzosa only)
if (any(meta$antibiotic == "Yes", na.rm = TRUE)) {
  permanova_abx <- imap_dfr(X_thr, function(X, nm) {
    ids <- intersect(rownames(X), meta$sample_id[meta$baseline & (meta$antibiotic != "Yes" | is.na(meta$antibiotic))])
    d <- meta[match(ids, meta$sample_id), ]
    dists <- list(aitchison = dist(clr_tf(X)[ids, ]),
                  bray      = vegdist(X[ids, ], method = "bray"))
    imap_dfr(dists, function(D, dn) {
      set.seed(seed)
      a <- adonis2(D ~ diagnosis, data = d, permutations = 999)
      tibble(threshold = nm, distance = dn, n = length(ids), R2 = a$R2[1], p = a$`Pr(>F)`[1])
    })
  })
  print(permanova_abx %>% arrange(distance, threshold))
}

# PERMANOVA with centre/sub-cohort: diagnosis effect adjusted for centre (by = "margin")
base_meta <- meta %>% filter(baseline)
if (n_distinct(base_meta$site) > 1 &&
    all(table(base_meta$site, base_meta$diagnosis)[, "control"] < nrow(base_meta))) {
  permanova_site <- imap_dfr(X_thr, function(X, nm) {
    ids <- intersect(rownames(X), base_meta$sample_id)
    d <- meta[match(ids, meta$sample_id), ]
    imap_dfr(list(aitchison = dist(clr_tf(X)[ids, ]), bray = vegdist(X[ids, ], method = "bray")), function(D, dn) {
      set.seed(seed)
      a <- tryCatch(adonis2(D ~ site + diagnosis, data = d, by = "margin", permutations = 999), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      tibble(threshold = nm, distance = dn, R2_diagnosis = a["diagnosis", "R2"], p_diagnosis = a["diagnosis", "Pr(>F)"],
             R2_site = a["site", "R2"], p_site = a["site", "Pr(>F)"])
    })
  })
  print(permanova_site %>% arrange(distance, threshold))
}

# ---- 5b. Batch: check and correct (for visualisation only) ------------
# batch = site (Hall: Boston/Atlanta, HMP2: centre C/M, Franzosa: PRISM/Validation)
# decision path:
#   1. design: does every batch contain both diagnoses? If not, batch is confounded
#      with diagnosis and no correction can separate the two.
#   2. size: how much variation does batch explain compared with diagnosis? (permanova_site)
#   3. tests: batch as covariate in the model (4g, permanova_site), do NOT test on corrected data
#      (correction with diagnosis in the model makes p-values too optimistic).
#   4. visualisation: remove the batch effect while protecting the diagnosis effect, then redo PCA/UMAP.
#      Three Bioconductor methods for comparison:
#      - limma::removeBatchEffect (on CLR): estimates a batch mean per genus and subtracts it
#      - sva::ComBat (on CLR): additionally aligns the spread per batch (empirical Bayes)
#      - MMUPHin::adjust_batch (on relative abundances): ComBat variant for microbiome data
#        with many zeros (zeros stay zeros); then CLR as usual
#      All three ignore that samples of the same person belong together -> for plots only.
X_main <- X_thr[[main_thr]]
Z_main <- clr_tf(X_main)
m_main <- meta[match(rownames(Z_main), meta$sample_id), ]

batch_design <- m_main %>% distinct(subject_id, diagnosis, site) %>%
  count(site, diagnosis, name = "persons") %>%
  pivot_wider(names_from = diagnosis, values_from = persons, values_fill = 0)
batch_design
batch_confounded <- any(batch_design$control == 0 | batch_design$cd == 0)
if (n_distinct(m_main$site) > 1 && batch_confounded)
  warning("batch and diagnosis are partly confounded (one batch has only one diagnosis) -> correction unreliable")

if (n_distinct(m_main$site) > 1) {
  design <- model.matrix(~ diagnosis, data = m_main)      # diagnosis is protected

  # limma: genera in rows, samples in columns -> transpose
  Z_limma <- t(removeBatchEffect(t(Z_main), batch = m_main$site, design = design))

  Z_combat <- t(ComBat(dat = t(Z_main), batch = m_main$site, mod = design))

  mm_data <- m_main %>% mutate(site = factor(site)) %>% column_to_rownames("sample_id")
  mmuphin <- adjust_batch(feature_abd = t(X_main), batch = "site", covariates = "diagnosis",
                          data = mm_data, control = list(verbose = FALSE))
  X_mmuphin <- t(mmuphin$feature_abd_adj)
  X_mmuphin <- X_mmuphin / rowSums(X_mmuphin)
  Z_mmuphin <- clr_tf(X_mmuphin)

  Z_list <- list(before = Z_main, limma = Z_limma, ComBat = Z_combat, MMUPHin = Z_mmuphin)

  ord_batch <- map_dfr(Z_list, function(Z) {
    pca <- prcomp(Z, center = TRUE, scale. = FALSE)
    set.seed(seed)
    um <- umap(Z, n_neighbors = min(15, nrow(Z) - 1), min_dist = 0.3)
    tibble(sample_id = rownames(Z), PC1 = pca$x[, 1], PC2 = pca$x[, 2], UMAP1 = um[, 1], UMAP2 = um[, 2])
  }, .id = "correction") %>%
    mutate(correction = factor(correction, levels = names(Z_list))) %>%
    left_join(meta, by = "sample_id")

  p_ub <- ggplot(ord_batch, aes(UMAP1, UMAP2, colour = diagnosis, shape = site)) +
    geom_point(alpha = .7) + facet_wrap(~ correction, scales = "free", ncol = 2) +
    labs(title = paste0("UMAP before/after batch correction (", main_thr, ")"),
         subtitle = if (batch_confounded) "Note: batch and diagnosis partly confounded" else NULL)
  print(p_ub)
  ggsave(file.path(outdir, "umap_batch_before_after.png"), p_ub, width = 11, height = 9)

  p_pb <- ggplot(ord_batch, aes(PC1, PC2, colour = diagnosis, shape = site)) +
    geom_point(alpha = .7) + facet_wrap(~ correction, scales = "free", ncol = 2) +
    labs(title = paste0("PCA before/after batch correction (", main_thr, ")"))
  print(p_pb)
  ggsave(file.path(outdir, "pca_batch_before_after.png"), p_pb, width = 11, height = 9)

  # how much variation does batch explain before/after correction? (baseline, Aitchison)
  # descriptive only - NO test of the diagnosis effect on corrected data
  batch_r2 <- map_dfr(Z_list, function(Z) {
    ids <- intersect(rownames(Z), meta$sample_id[meta$baseline])
    d <- meta[match(ids, meta$sample_id), ]
    set.seed(seed)
    a <- adonis2(dist(Z[ids, ]) ~ site, data = d, permutations = 999)
    tibble(R2_batch = a$R2[1], p_batch = a$`Pr(>F)`[1])
  }, .id = "correction")
  print(batch_r2)
}

# main figures (UMAP + PCA): ComBat-corrected if every batch has both diagnoses;
# otherwise uncorrected (Hall: Atlanta has only controls -> correction would separate artificially)
use_combat <- n_distinct(m_main$site) > 1 && !batch_confounded
Z_fig <- if (use_combat) Z_combat else Z_main
label_fig <- if (use_combat) "ComBat-corrected, batch = site" else "without batch correction"
note_fig <- c(
  if (use_combat && min(as.matrix(batch_design[, c("control", "cd")])) < 3)
    "Note: one batch has fewer than 3 persons of one diagnosis",
  if (!use_combat && n_distinct(m_main$site) > 1)
    "no correction: batch and diagnosis confounded",
  if (longitudinal) "all time points; dense clusters are usually one person")

pca_h <- prcomp(Z_fig, center = TRUE, scale. = FALSE)
var_h <- round(100 * summary(pca_h)$importance[2, 1:2], 1)
set.seed(seed)
um_h <- umap(Z_fig, n_neighbors = min(15, nrow(Z_fig) - 1), min_dist = 0.3)
ord_fig <- tibble(sample_id = rownames(Z_fig), PC1 = pca_h$x[, 1], PC2 = pca_h$x[, 2],
                    UMAP1 = um_h[, 1], UMAP2 = um_h[, 2]) %>%
  left_join(meta, by = "sample_id")

p_uh <- ggplot(ord_fig, aes(UMAP1, UMAP2, colour = diagnosis, shape = site)) +
  geom_point(alpha = .75, size = 1.8) +
  labs(title = paste0("UMAP ", cohort, " (", main_thr, ", ", label_fig, ")"),
       subtitle = paste(note_fig, collapse = "\n"))
print(p_uh)
ggsave(file.path(outdir, "umap_main.png"), p_uh, width = 8, height = 6)

p_ph <- ggplot(ord_fig, aes(PC1, PC2, colour = diagnosis, shape = site)) +
  geom_point(alpha = .75, size = 1.8) +
  labs(x = paste0("PC1 (", var_h[1], " %)"), y = paste0("PC2 (", var_h[2], " %)"),
       title = paste0("PCA ", cohort, " (", main_thr, ", ", label_fig, ")"),
       subtitle = paste(note_fig, collapse = "\n"))
print(p_ph)
ggsave(file.path(outdir, "pca_main.png"), p_ph, width = 8, height = 6)

# UMAP coloured by person: are clusters driven by person rather than diagnosis?
if (longitudinal) {
  p_up <- ord %>% filter(threshold == main_thr) %>%
    ggplot(aes(UMAP1, UMAP2, colour = subject_id, shape = diagnosis)) +
    geom_point(alpha = .8) + guides(colour = "none") +
    labs(title = paste0("UMAP coloured by person (", main_thr, ")"))
  print(p_up)
  ggsave(file.path(outdir, "umap_person.png"), p_up, width = 7, height = 5)
}

# ---- 6. Who stands out? -----------------------------------------------
# at how many thresholds is a sample an outlier, and is it also a few-genera sample?
outliers <- ord %>%
  group_by(sample_id, subject_id, diagnosis) %>%
  summarise(n_thr_outlier = sum(outlier),
            ever_few_genera = any(few_genera),
            min_n_genera = min(n_genera),
            total_raw = rs[first(sample_id)], .groups = "drop") %>%
  filter(n_thr_outlier > 0) %>%
  arrange(desc(n_thr_outlier))
outliers

# sensitivity: alpha diversity test without the outliers (decide only after looking!)
# (without outliers some persons lose their baseline sample -> that person drops out)
out_ids <- outliers$sample_id
alpha %>% filter(baseline) %>%
  group_by(threshold) %>%
  summarise(n_with = n(), p_with = wilcox.test(Shannon ~ diagnosis)$p.value,
            n_without = sum(!sample_id %in% out_ids),
            p_without = wilcox.test(Shannon[!sample_id %in% out_ids] ~ diagnosis[!sample_id %in% out_ids])$p.value)

# save results
saveRDS(list(sparsity = sparsity_tab, few = few, alpha = alpha, ord = ord,
             outliers = outliers, permanova = permanova, pcoa = pcoa,
             mixed = if (longitudinal) mixed_tab else NULL,
             alpha_person = if (longitudinal) alpha_person else NULL,
             median_test = if (longitudinal) median_tab else NULL,
             alpha_var = if (longitudinal) alpha_var else NULL,
             stability = if (longitudinal) stability else NULL,
             time_model = if (longitudinal) time_tab else NULL,
             time_model_slopes = if (longitudinal) time_slope_tab else NULL,
             person_slopes = if (longitudinal) person_slopes else NULL,
             site = site_tab,
             permanova_site = if (exists("permanova_site")) permanova_site else NULL,
             batch_design = batch_design,
             batch_r2 = if (exists("batch_r2")) batch_r2 else NULL,
             ord_batch = if (exists("ord_batch")) ord_batch else NULL,
             ord_fig = ord_fig, fig_correction = label_fig,
             age_model = if (exists("age_tab")) age_tab else NULL,
             antibiotics = if (exists("abx_tab")) abx_tab else NULL,
             permanova_without_abx = if (exists("permanova_abx")) permanova_abx else NULL),
        file.path(outdir, "abundance_results.rds"))
