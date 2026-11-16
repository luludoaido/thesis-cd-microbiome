# --- 04 Load abundance ---
# Loads one cohort and prepares relative abundances for the diversity analyses. Not run on its
# own: the 05 scripts set "cohort" and "level" and source this file.

# Results: no files. Gives the 05 scripts meta, X_tss, X_thr (one matrix per filter: 0, 5, 10 %),
# main_thr, few (samples with very few taxa) and longitudinal (several samples per person)

# Methods: docs/methods.md, 1.3, 1.4, 7.1, 7.2

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
# prevalence filter (see methods 1.3) at 0, 5 and 10 % (see methods 7.1), then scaled back to 1 (see methods 1.4)
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
# lowest 5 % (see methods 7.2); n_genera counts species at species level.
few <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X), n_genera = rowSums(X > 0), threshold = nm)
}) %>%
  left_join(meta, by = "sample_id") %>%
  group_by(threshold) %>%
  mutate(few_genera = n_genera <= quantile(n_genera, 0.05)) %>%
  ungroup()

# --- Helper ---
# shows a plot and saves it
save_plot <- function(plot, name, width, height){
  print(plot)
  ggsave(file.path(outdir, name), plot, width = width, height = height)
}
