# --- Fibre-associated taxa (SQ2a and SQ2b) ---
#
# What this script does
# SQ2a: is the Prevotella/Bacteroides ratio lower in CD than in controls?
# SQ2b: are F. prausnitzii and Roseburia lower in CD than in controls?
# Only Franzosa and HMP2. Run from the project folder.

# Decisions
# 1. Taxa were chosen in advance, so no prevalence filter for abundances and ratios.
# 2. Prevotella is absent from about two thirds of the samples, and a replacement value would
#    decide the result. SQ2a is therefore tested in parts: occurrence, ratio where both taxa
#    occur, and Prevotella-dominated samples.
# 3. SQ2b: MaAsLin2 is the main method. Abundances are TSS-scaled over all species beforehand,
#    so MaAsLin2 runs with normalization NONE, transform LOG, method LM, and no prevalence filter.
# 4. SQ2b check: linear models on CLR values, with the 10 % genus filter as in 05_beta_diversity.R.
# 5. Holm correction over the two SQ2b taxa, per cohort.

library(tidyverse)
library(lme4)
library(lmerTest)
library(Maaslin2)
source("R/00_setup.R")

outdir <- file.path(dir_results, "hypotheses")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("franzosa", "hmp2")
main_thr <- 0.10
hmp2_tax <- read_csv(file.path(dir_processed, "hmp2/hmp2_taxonomy.csv"), show_col_types = FALSE)

genus_of <- function(cohort, species) {
  if (cohort == "hmp2") hmp2_tax$genus[match(species, hmp2_tax$species)] else sub("_.*", "", species)
}

# --- Load ---
load_cohort <- function(cohort) {
  folder <- file.path(dir_processed, cohort)
  meta <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE) %>%
    mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))
  X <- read_csv(file.path(folder, paste0(cohort, "_species.csv")), show_col_types = FALSE) %>%
    filter(sample_id %in% meta$sample_id) %>%
    arrange(match(sample_id, meta$sample_id)) %>%
    column_to_rownames("sample_id") %>% as.matrix()
  list(meta = meta, X = X / rowSums(X))
}
data <- set_names(map(cohorts, load_cohort), cohorts)

fp_name <- "Faecalibacterium_prausnitzii"

taxa <- imap_dfr(data, function(d, cohort) {
  g <- genus_of(cohort, colnames(d$X))
  stopifnot(fp_name %in% colnames(d$X))
  tibble(sample_id = rownames(d$X),
         Prevotella  = rowSums(d$X[, g %in% "Prevotella", drop = FALSE]),
         Bacteroides = rowSums(d$X[, g %in% "Bacteroides", drop = FALSE]),
         F_prausnitzii = d$X[, fp_name],
         Roseburia = rowSums(d$X[, g %in% "Roseburia", drop = FALSE])) %>%
    left_join(d$meta %>% select(sample_id, subject_id, diagnosis, site, baseline, antibiotic), by = "sample_id") %>%
    mutate(cohort = cohort)
})

# --- Detection ---
detection <- taxa %>%
  pivot_longer(c(Prevotella, Bacteroides, F_prausnitzii, Roseburia), names_to = "taxon", values_to = "abund") %>%
  group_by(cohort, diagnosis, taxon) %>%
  summarise(samples = n(), persons = n_distinct(subject_id), prevalence_pct = mean(abund > 0) * 100,
            median_pct = median(abund) * 100, mean_pct = mean(abund) * 100, .groups = "drop") %>%
  arrange(cohort, taxon, diagnosis)
print(detection, n = Inf)
write_csv(detection, file.path(outdir, "sq2_taxa_detection.csv"))

fr <- taxa %>% filter(cohort == "franzosa")
hm <- taxa %>% filter(cohort == "hmp2")

# --- SQ2a, part 1: is Prevotella detected? ---
share <- function(x, dx, g) mean(x[dx == g]) * 100

sq2a_occurrence <- bind_rows(
  {
    d <- fr %>% mutate(present = Prevotella > 0)
    fit <- glm(present ~ diagnosis + site, data = d, family = binomial)
    tibble(cohort = "franzosa", unit = "samples", n = nrow(d),
           control_pct = share(d$present, d$diagnosis, "control"), cd_pct = share(d$present, d$diagnosis, "cd"),
           p_fisher = fisher.test(table(d$diagnosis, d$present))$p.value,
           model = "logistic: diagnosis + site",
           model_estimate = coef(summary(fit))["diagnosiscd", "Estimate"],
           model_p = coef(summary(fit))["diagnosiscd", "Pr(>|z|)"])
  },
  {
    d <- hm %>% mutate(present = Prevotella > 0)
    pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(ever = any(present), .groups = "drop")
    fit <- glmer(present ~ diagnosis + site + (1 | subject_id), data = d, family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    tibble(cohort = "hmp2", unit = "persons (ever detected)", n = nrow(pp),
           control_pct = share(pp$ever, pp$diagnosis, "control"), cd_pct = share(pp$ever, pp$diagnosis, "cd"),
           p_fisher = fisher.test(table(pp$diagnosis, pp$ever))$p.value,
           model = "logistic mixed: diagnosis + site + (1 | person)",
           model_estimate = fixef(fit)[["diagnosiscd"]],
           model_p = coef(summary(fit))["diagnosiscd", "Pr(>|z|)"])
  })
print(sq2a_occurrence, width = Inf)
write_csv(sq2a_occurrence, file.path(outdir, "sq2a_occurrence.csv"))

# --- SQ2a, part 2: ratio in samples with both taxa ---
both <- taxa %>% filter(Prevotella > 0, Bacteroides > 0) %>%
  mutate(log_ratio = log10(Prevotella / Bacteroides))

sq2a_ratio <- bind_rows(
  {
    d <- both %>% filter(cohort == "franzosa")
    s <- coef(summary(lm(log_ratio ~ diagnosis + site, data = d)))
    tibble(cohort = "franzosa", samples = nrow(d), persons = n_distinct(d$subject_id),
           n_control = sum(d$diagnosis == "control"), n_cd = sum(d$diagnosis == "cd"),
           median_control = median(d$log_ratio[d$diagnosis == "control"]),
           median_cd = median(d$log_ratio[d$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(log_ratio ~ diagnosis, data = d, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], model_p = s["diagnosiscd", "Pr(>|t|)"])
  },
  {
    d <- both %>% filter(cohort == "hmp2")
    pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(log_ratio), .groups = "drop")
    s <- coef(summary(lmer(log_ratio ~ diagnosis + site + (1 | subject_id), data = d)))
    tibble(cohort = "hmp2", samples = nrow(d), persons = n_distinct(d$subject_id),
           n_control = sum(pp$diagnosis == "control"), n_cd = sum(pp$diagnosis == "cd"),
           median_control = median(pp$v[pp$diagnosis == "control"]),
           median_cd = median(pp$v[pp$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(v ~ diagnosis, data = pp, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], model_p = s["diagnosiscd", "Pr(>|t|)"])
  })
print(sq2a_ratio, width = Inf)
write_csv(sq2a_ratio, file.path(outdir, "sq2a_ratio.csv"))

# --- SQ2a, part 3: Prevotella > Bacteroides ---
sq2a_dominated <- bind_rows(
  {
    d <- fr %>% mutate(dom = Prevotella > Bacteroides)
    tibble(cohort = "franzosa", unit = "samples", n = nrow(d),
           control_pct = share(d$dom, d$diagnosis, "control"), cd_pct = share(d$dom, d$diagnosis, "cd"),
           p_fisher = fisher.test(table(d$diagnosis, d$dom))$p.value)
  },
  {
    pp <- hm %>% group_by(subject_id, diagnosis) %>%
      summarise(dom = any(Prevotella > Bacteroides), .groups = "drop")
    tibble(cohort = "hmp2", unit = "persons (at least one sample)", n = nrow(pp),
           control_pct = share(pp$dom, pp$diagnosis, "control"), cd_pct = share(pp$dom, pp$diagnosis, "cd"),
           p_fisher = fisher.test(table(pp$diagnosis, pp$dom))$p.value)
  })
print(sq2a_dominated, width = Inf)
write_csv(sq2a_dominated, file.path(outdir, "sq2a_dominated.csv"))

# --- SQ2b: MaAsLin2 (main) ---
run_maaslin <- function(d, label, random) {
  run_dir <- file.path(outdir, "maaslin2", label)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  feat <- d %>% select(sample_id, F_prausnitzii, Roseburia) %>% column_to_rownames("sample_id") %>% as.data.frame()
  m <- d %>% select(sample_id, diagnosis, site, subject_id) %>%
    mutate(diagnosis = as.character(diagnosis)) %>% column_to_rownames("sample_id") %>% as.data.frame()
  Maaslin2(input_data = feat, input_metadata = m, output = run_dir,
           min_abundance = 0, min_prevalence = 0, min_variance = 0,
           normalization = "NONE", transform = "LOG", analysis_method = "LM",
           correction = "holm", standardize = FALSE,
           fixed_effects = c("diagnosis", "site"), random_effects = random,
           reference = "diagnosis,control",
           plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  # MaAsLin2 corrects over all terms including site, the hypothesis is only about diagnosis
  read_tsv(file.path(run_dir, "all_results.tsv"), show_col_types = FALSE) %>% filter(metadata == "diagnosis") %>%
    transmute(analysis = label, taxon = feature, value, samples = N, not_zero = N.not.0, coef, stderr, pval,
              p_holm = p.adjust(pval, "holm"))
}

sq2b_maaslin <- bind_rows(
  run_maaslin(fr, "franzosa", NULL),
  run_maaslin(fr %>% filter(antibiotic == "No"), "franzosa_without_antibiotics", NULL),
  run_maaslin(hm, "hmp2", "subject_id"))
print(sq2b_maaslin, n = Inf, width = Inf)
write_csv(sq2b_maaslin, file.path(outdir, "sq2b_maaslin2.csv"))

# --- SQ2b: CLR models (check) ---
# same CLR as in 05_beta_diversity.R
clr_tf <- function(X) {
  pc <- min(X[X > 0]) / 2
  L <- log(X + pc)
  L - rowMeans(L)
}

clr_values <- imap_dfr(data, function(d, cohort) {
  X <- d$X[, passes_filter(d$X, cohort, "species", main_thr), drop = FALSE]
  g <- genus_of(cohort, colnames(X))
  is_ros <- g %in% "Roseburia"
  # all Roseburia species become one part, because the hypothesis is about the genus
  parts <- cbind(X[, !is_ros, drop = FALSE], Roseburia = rowSums(X[, is_ros, drop = FALSE]))
  parts <- parts / rowSums(parts)
  Z <- clr_tf(parts)
  tibble(sample_id = rownames(Z), F_prausnitzii = Z[, fp_name], Roseburia = Z[, "Roseburia"], cohort = cohort)
}) %>% left_join(taxa %>% select(sample_id, subject_id, diagnosis, site, antibiotic), by = "sample_id")

sq2b_clr <- clr_values %>%
  pivot_longer(c(F_prausnitzii, Roseburia), names_to = "taxon", values_to = "clr") %>%
  group_split(cohort, taxon) %>%
  map_dfr(function(d) {
    if (d$cohort[1] == "franzosa") {
      s <- coef(summary(lm(clr ~ diagnosis + site, data = d)))
      v <- d %>% transmute(diagnosis, v = clr)
    } else {
      s <- coef(summary(lmer(clr ~ diagnosis + site + (1 | subject_id), data = d)))
      v <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(clr), .groups = "drop")
    }
    tibble(cohort = d$cohort[1], taxon = d$taxon[1], n = nrow(d),
           median_control = median(v$v[v$diagnosis == "control"]),
           median_cd = median(v$v[v$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(v ~ diagnosis, data = v, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], model_se = s["diagnosiscd", "Std. Error"],
           model_p = s["diagnosiscd", grep("^Pr", colnames(s), value = TRUE)])
  }) %>%
  group_by(cohort) %>% mutate(model_p_holm = p.adjust(model_p, "holm")) %>% ungroup()
print(sq2b_clr, n = Inf, width = Inf)
write_csv(sq2b_clr, file.path(outdir, "sq2b_clr.csv"))

# --- Plot ---
plot_data <- taxa %>%
  group_by(cohort, subject_id, diagnosis) %>%
  summarise(across(c(F_prausnitzii, Roseburia, Prevotella, Bacteroides), median), .groups = "drop") %>%
  pivot_longer(c(F_prausnitzii, Roseburia, Prevotella, Bacteroides), names_to = "taxon", values_to = "abund") %>%
  mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
         taxon = factor(taxon, c("F_prausnitzii", "Roseburia", "Prevotella", "Bacteroides"),
                        c("F. prausnitzii", "Roseburia", "Prevotella", "Bacteroides")),
         pct = pmax(abund * 100, 1e-4))

p <- ggplot(plot_data, aes(diagnosis, pct, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = .55) +
  geom_jitter(width = .15, alpha = .6, size = 1.2) +
  scale_y_log10() +
  scale_colour_manual(values = col_diag, guide = "none") +
  facet_grid(cohort ~ taxon) +
  labs(x = NULL, y = "relative abundance (%, log scale, zeros shown at the bottom)",
       title = "SQ2 taxa in CD and controls (median per person)") +
  theme_plot
ggsave(file.path(outdir, "sq2_taxa_cd_vs_control.png"), p, width = 10, height = 7, dpi = 200)

# H2b only: the two butyrate producers
plot_h2b <- plot_data %>%
  filter(taxon %in% c("F. prausnitzii", "Roseburia")) %>%
  mutate(taxon = droplevels(taxon),
         diagnosis = factor(diagnosis, c("control", "cd"), c("Control", "CD")))

p_h2b <- ggplot(plot_h2b, aes(diagnosis, pct, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = .55) +
  geom_jitter(width = .15, alpha = .6, size = 1.2) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = setNames(unname(col_diag), c("Control", "CD")), guide = "none") +
  facet_grid(cohort ~ taxon) +
  labs(x = NULL, y = "relative abundance (%, log scale)",
       title = "Butyrate producers in controls and CD (median per person)",
       caption = "Persons without the taxon are shown at the bottom") +
  theme_plot
ggsave(file.path(outdir, "sq2b_taxa_cd_vs_control.png"), p_h2b, width = 7, height = 6.5, dpi = 200)

# --- Literature ---
# Aitchison J (1986). The statistical analysis of compositional data. doi:10.1007/978-94-009-4109-0
# Mallick H et al. (2021). Multivariable association discovery in population-scale meta-omics
#   studies. PLoS Comput Biol 17(11):e1009442. doi:10.1371/journal.pcbi.1009442
