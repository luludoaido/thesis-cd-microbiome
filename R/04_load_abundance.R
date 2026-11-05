# --- Load abundance ---
# loads the cleaned tables of one cohort form 02_clean_*.R and applies the prevalece filter
# not run on its own - 05_alpha_diversity.R and 05_beta_diversity.R set "cohort" and "level" and then sources this file

# genus
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

longitudinal <- max(table(meta$subject_id)) > 1

# --- TSS and prevalence filter ---
rs <- rowSums(X_raw)
X_tss <- X_raw / rs

thresholds <- c(0, 0.05, 0.10)

# filter at genus level (see passes_filter in 00_setup.R); 
# thr = 0 only removes taxa that occur in no sample
X_thr <- lapply(thresholds, function(t) {
  Xf <- X_tss[, passes_filter(X_tss, cohort, level, t), drop = FALSE]
  # drops samples with no taxa left
  Xf <- Xf[rowSums(Xf) > 0, , drop = FALSE]      
  Xf / rowSums(Xf)
  # filtering lowers the sum, so scla back to 1
})
names(X_thr) <- sprintf("thr_%02d", thresholds * 100)   # thr_00, thr_05, thr_10

# main analysis: 10 % filter (fewer zeros -> CLR depends less on the pseudocount);
# 0 % and 5 % run as sensitivity analysis. All plots use main_thr.
main_thr <- "thr_10"

# --- Samples with very few genera ---
# number of genera per sample per threshold
few <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X), n_genera = rowSums(X > 0), threshold = nm)
}) %>%
  left_join(meta, by = "sample_id") %>%
  group_by(threshold) %>%
  mutate(few_genera = n_genera <= quantile(n_genera, 0.05)) %>%     # lowest 5 % per threshold
  ungroup()

# show a plot and save it
save_plot <- function(plot, name, width, height){
  print(plot)
  ggsave(file.path(outdir, name), plot, width = width, height = height)
}
