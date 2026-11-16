# --- Load abundance ---

# What this script does
# Loads the cleaned tables of one cohort and prepares the relative abundances for the
# diversity analyses. It is not run on its own: 05_alpha_diversity.R and 05_beta_diversity.R
# set "cohort" and "level" and then source this file, so both start from the same data.

# Decisions
# 1. Relative abundances (total sum scaling): samples differ in sequencing depth, so every
#    sample is divided by its own total.
# 2. Prevalence filter at genus level, see passes_filter in 00_setup.R. A genus stays if it
#    occurs in enough samples, and a species stays if its genus stays. Own decision.
# 3. Three filters: 10 % is the main analysis, 0 % and 5 % are sensitivity analyses.
# 4. After filtering, every sample is scaled back to a total of 1.
# 5. Samples among the lowest 5 % in the number of taxa are only flagged, not removed.

# --- Load the data ---
# Genus results go to results/<cohort>, species results to results/species/<cohort>.
outdir <- if (level == "genus") file.path(dir_results, cohort) else
  file.path(dir_results, "species", cohort)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

folder <- file.path(dir_processed, cohort)
meta   <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))

abund  <- read_csv(file.path(folder, paste0(cohort, "_", level, ".csv")), show_col_types = FALSE)

# abundance matrix in the same order as the metadata
X_raw <- abund %>%
  arrange(match(sample_id, meta$sample_id)) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# Hall and HMP2 have several samples per person, Franzosa has one.
# The 05 scripts run the analyses over time only if this is TRUE.
longitudinal <- max(table(meta$subject_id)) > 1

# --- TSS and prevalence filter ---
# Filtering taxa that are rare across samples reduces technical variation while keeping the
# differences between samples (Cao 2021). The exact cut-off is an own decision. Fewer zeros
# also make the CLR transformation in 05_beta_diversity.R depend less on the pseudocount.
# Results can depend on such pre-processing steps (Nearing 2022, shown for 16S data),
# so all analyses also run with 0 % and 5 %.
rs <- rowSums(X_raw)
X_tss <- X_raw / rs

thresholds <- c(0, 0.05, 0.10)

# thr = 0 only removes taxa that occur in no sample
X_thr <- lapply(thresholds, function(t) {
  Xf <- X_tss[, passes_filter(X_tss, cohort, level, t), drop = FALSE]
  # drops samples with no taxa left
  Xf <- Xf[rowSums(Xf) > 0, , drop = FALSE]
  # filtering lowers the sum, so scale back to 1
  Xf / rowSums(Xf)
})
names(X_thr) <- sprintf("thr_%02d", thresholds * 100)   # thr_00, thr_05, thr_10

main_thr <- "thr_10"

# --- Samples with very few taxa ---
# Flags the samples with the fewest taxa, so that the plots can show whether they stand out.
# The 5 % limit is an own decision for this check. n_genera counts species at species level.
few <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X), n_genera = rowSums(X > 0), threshold = nm)
}) %>%
  left_join(meta, by = "sample_id") %>%
  group_by(threshold) %>%
  mutate(few_genera = n_genera <= quantile(n_genera, 0.05)) %>%   
  ungroup()

# --- Helper ---
# shows a plot and save it
save_plot <- function(plot, name, width, height){
  print(plot)
  ggsave(file.path(outdir, name), plot, width = width, height = height)
}

# --- Literature ---
# Cao Q et al. (2021). Effects of rare microbiome taxa filtering on statistical analysis.
#   Front Microbiol 11:607325. doi:10.3389/fmicb.2020.607325
# Nearing JT et al. (2022). Microbiome differential abundance methods produce different results
#   across 38 datasets. Nat Commun 13:342. doi:10.1038/s41467-022-28034-z
