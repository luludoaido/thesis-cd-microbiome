# --- 11 Disease severity and taxa within CD ---
# Within CD only: which genera and how much dysbiosis go with higher calprotectin?
# Franzosa and HMP2 analysed separately and compared side by side; exploratory (see methods 13).

# Results (per cohort):
# - genera vs. log10 calprotectin (MaAsLin2, 95 % CI, BH q); HMP2 check with HBI
# - MD-index (Pascal and Gevers stool) vs. log10 calprotectin (lm / lmer)
# - agreement of the genus effects between the cohorts; forest plot of genera with q < 0.25

# Methods: docs/methods.md, section 13

library(tidyverse)
library(lme4)
library(lmerTest)
library(Maaslin2)
source("R/00_setup.R")

outdir <- file.path(dir_results, "hypotheses")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cohorts  <- c("franzosa", "hmp2")
main_thr <- 0.10
q_cut    <- 0.25   # MaAsLin2 default (see methods 13.1)

# --- Load: CD only, genus level, 10 % filter (see methods 1.3) ---
load_cd <- function(cohort) {
  folder <- file.path(dir_processed, cohort)
  meta <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE)
  if (cohort == "franzosa") meta <- meta %>% mutate(calpro = calprotectin, hbi = NA_real_)
  X <- read_csv(file.path(folder, paste0(cohort, "_genus.csv")), show_col_types = FALSE) %>%
    filter(sample_id %in% meta$sample_id) %>%
    arrange(match(sample_id, meta$sample_id)) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  X <- X / rowSums(X)
  X <- X[, passes_filter(X, cohort, "genus", main_thr), drop = FALSE]   # filter on all samples
  cd <- meta$diagnosis == "cd"
  # values below 1 set to 1 for the log, as in 9.2
  list(meta = meta[cd, ] %>% mutate(log_calpro = log10(pmax(calpro, 1))), X = X[cd, , drop = FALSE])
}
data <- set_names(map(cohorts, load_cd), cohorts)

# --- Genera vs. severity: MaAsLin2, same settings as H2b (see methods 11.3, 13.1) ---
run_severity <- function(cohort, predictor, label, keep = TRUE) {
  d <- data[[cohort]]
  m <- d$meta %>% filter(keep, !is.na(.data[[predictor]]))
  fixed <- c(predictor, if (n_distinct(m$site) > 1) "site")
  run_dir <- file.path(outdir, "maaslin2", paste0(cohort, "_severity_", label))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  Maaslin2(input_data = as.data.frame(d$X[m$sample_id, , drop = FALSE]),
           input_metadata = m %>% select(sample_id, all_of(fixed), subject_id) %>%
             column_to_rownames("sample_id") %>% as.data.frame(),
           output = run_dir, min_abundance = 0, min_prevalence = 0, min_variance = 0,
           normalization = "NONE", transform = "LOG", analysis_method = "LM",
           correction = "BH", standardize = FALSE, fixed_effects = fixed,
           random_effects = if (cohort == "hmp2") "subject_id" else NULL,
           plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  read_tsv(file.path(run_dir, "all_results.tsv"), show_col_types = FALSE) %>%
    filter(metadata == predictor) %>%
    transmute(cohort, analysis = label, predictor, genus = feature, samples = N,
              persons = n_distinct(m$subject_id), not_zero = N.not.0, coef, stderr,
              ci_low = coef - qnorm(0.975) * stderr, ci_high = coef + qnorm(0.975) * stderr,
              pval, qval = p.adjust(pval, "BH"))
}

fr_no_ab <- data$franzosa$meta$antibiotic %in% "No"
severity <- bind_rows(
  run_severity("franzosa", "log_calpro", "calprotectin"),
  run_severity("franzosa", "log_calpro", "calprotectin_no_antibiotics", keep = fr_no_ab),
  run_severity("hmp2", "log_calpro", "calprotectin"),
  run_severity("hmp2", "hbi", "hbi"))
write_csv(severity, file.path(outdir, "severity_taxa.csv"))
severity %>% filter(qval < q_cut) %>% arrange(cohort, analysis, qval) %>% print(n = Inf, width = Inf)

# --- Do the cohorts agree? genera tested in both, main analysis ---
agree <- severity %>% filter(analysis == "calprotectin") %>%
  select(cohort, genus, coef) %>% pivot_wider(names_from = cohort, values_from = coef) %>% drop_na()
agreement <- tibble(genera_in_both = nrow(agree),
                    same_direction = sum(sign(agree$franzosa) == sign(agree$hmp2)),
                    spearman_rho = cor(agree$franzosa, agree$hmp2, method = "spearman"),
                    spearman_p = cor.test(agree$franzosa, agree$hmp2, method = "spearman", exact = FALSE)$p.value)
print(agreement)
write_csv(agreement, file.path(outdir, "severity_agreement.csv"))

# --- MD-index vs. severity (Gevers 2014: MD-index rises with severity; see methods 13.2) ---
# both main versions of the MD-index (see methods 10.2)
md_scores <- read_csv(file.path(outdir, "dysbiosis_scores.csv"), show_col_types = FALSE) %>%
  select(sample_id, md_pascal, md_gevers_stool)
md_severity <- expand_grid(cohort = cohorts, version = c("md_pascal", "md_gevers_stool")) %>%
  pmap_dfr(function(cohort, version) {
    m <- data[[cohort]]$meta %>% inner_join(md_scores, by = "sample_id") %>%
      filter(!is.na(log_calpro)) %>% mutate(md = .data[[version]])
    fit <- if (cohort == "hmp2") lmer(md ~ log_calpro + site + (1 | subject_id), data = m)
           else lm(md ~ log_calpro + site, data = m)
    s <- coef(summary(fit))
    tibble(cohort, version, samples = nrow(m), persons = n_distinct(m$subject_id),
           estimate = s["log_calpro", "Estimate"], se = s["log_calpro", "Std. Error"],
           p = s["log_calpro", grep("^Pr", colnames(s))])
  })
print(md_severity)
write_csv(md_severity, file.path(outdir, "severity_md_index.csv"))

# --- Forest plot: genera with q < 0.25 in at least one cohort (main analysis) ---
hits <- severity %>% filter(analysis == "calprotectin", qval < q_cut) %>% pull(genus) %>% unique()
if (length(hits) > 0) {
  p_sev <- severity %>% filter(analysis == "calprotectin", genus %in% hits) %>%
    mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
           genus = fct_reorder(genus, coef),
           significant = if_else(qval < q_cut, "q < 0.25", "q ≥ 0.25")) %>%
    ggplot(aes(coef, genus, colour = cohort, shape = significant)) +
    geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y",
                  position = position_dodge(width = .6)) +
    geom_point(size = 2.4, position = position_dodge(width = .6)) +
    scale_colour_manual(values = unname(okabe_ito[c("orange", "sky_blue")])) +
    scale_shape_manual(values = c("q < 0.25" = 16, "q ≥ 0.25" = 1)) +
    labs(x = "change in log2 abundance per 10-fold higher calprotectin (95 % CI)", y = NULL,
         colour = NULL, shape = NULL, title = "Genera linked to calprotectin within CD") +
    theme_plot
  ggsave(file.path(outdir, "severity_taxa_forest.png"), p_sev,
         width = 8, height = 2 + 0.3 * length(hits), dpi = 200)
}
