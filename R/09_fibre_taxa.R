# --- 09 H2: Fibre-associated taxa ---
# H2a: is the Prevotella/Bacteroides ratio lower in CD? H2b: are F. prausnitzii and Roseburia
# lower in CD? Franzosa and HMP2 only; the taxa are not prevalence-filtered (see methods 11.1).

# Results:
# - detection and abundance of Prevotella, Bacteroides, F. prausnitzii and Roseburia per group
# - H2a: Prevotella detected, log10 ratio where both occur, Prevotella-dominated samples (CD vs. control);
#   HMP2 always with a mixed model on all samples, one value per person as check
# - H2a, part 4: plant/animal taxa ratio after David 2014 (all samples), same tests as above
# - plots of the H2a ratios per person: both together and each alone (sq2a_*_cd_vs_control.png)
# - H2b: CD effect on F. prausnitzii and Roseburia, MaAsLin2 (main) and CLR models (check), Holm-adjusted
# - H2b extension (exploratory): E. rectale, A. hadrus, E. hallii, Coprococcus, same tests

# Methods: docs/methods.md, 1.10 and section 11

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

# --- SQ2a, part 1: is Prevotella detected? (see methods 11.2) ---
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

# --- SQ2a, part 2: ratio in samples with both taxa (see methods 11.2) ---
both <- taxa %>% filter(Prevotella > 0, Bacteroides > 0) %>%
  mutate(log_ratio = log10(Prevotella / Bacteroides))

sq2a_ratio <- bind_rows(
  {
    d <- both %>% filter(cohort == "franzosa")
    fit <- lm(log_ratio ~ diagnosis + site, data = d); s <- coef(summary(fit))
    tibble(cohort = "franzosa", samples = nrow(d), persons = n_distinct(d$subject_id),
           n_control = sum(d$diagnosis == "control"), n_cd = sum(d$diagnosis == "cd"),
           median_control = median(d$log_ratio[d$diagnosis == "control"]),
           median_cd = median(d$log_ratio[d$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(log_ratio ~ diagnosis, data = d, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
           model_p = s["diagnosiscd", "Pr(>|t|)"])
  },
  {
    d <- both %>% filter(cohort == "hmp2")
    pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(log_ratio), .groups = "drop")
    fit <- lmer(log_ratio ~ diagnosis + site + (1 | subject_id), data = d); s <- coef(summary(fit))
    tibble(cohort = "hmp2", samples = nrow(d), persons = n_distinct(d$subject_id),
           n_control = sum(pp$diagnosis == "control"), n_cd = sum(pp$diagnosis == "cd"),
           median_control = median(pp$v[pp$diagnosis == "control"]),
           median_cd = median(pp$v[pp$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(v ~ diagnosis, data = pp, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
           model_p = s["diagnosiscd", "Pr(>|t|)"])
  })
print(sq2a_ratio, width = Inf)
write_csv(sq2a_ratio, file.path(outdir, "sq2a_ratio.csv"))

# --- SQ2a, part 3: Prevotella > Bacteroides (see methods 11.2) ---
sq2a_dominated <- bind_rows(
  {
    d <- fr %>% mutate(dom = Prevotella > Bacteroides)
    tibble(cohort = "franzosa", unit = "samples", n = nrow(d),
           control_pct = share(d$dom, d$diagnosis, "control"), cd_pct = share(d$dom, d$diagnosis, "cd"),
           p_fisher = fisher.test(table(d$diagnosis, d$dom))$p.value)
  },
  {
    # all samples with person as random effect (see methods 1.5); Fisher per person as check
    d <- hm %>% mutate(dom = Prevotella > Bacteroides)
    pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(dom = any(dom), .groups = "drop")
    fit <- glmer(dom ~ diagnosis + site + (1 | subject_id), data = d, family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    tibble(cohort = "hmp2", unit = "persons (at least one sample)", n = nrow(pp),
           control_pct = share(pp$dom, pp$diagnosis, "control"), cd_pct = share(pp$dom, pp$diagnosis, "cd"),
           p_fisher = fisher.test(table(pp$diagnosis, pp$dom))$p.value,
           model = "logistic mixed: diagnosis + site + (1 | person)",
           model_estimate = fixef(fit)[["diagnosiscd"]],
           model_p = coef(summary(fit))["diagnosiscd", "Pr(>|z|)"])
  })
print(sq2a_dominated, width = Inf)
write_csv(sq2a_dominated, file.path(outdir, "sq2a_dominated.csv"))

# --- SQ2a, part 4: plant/animal taxa ratio after David 2014 (see methods 11.7) ---
# log10 [(Roseburia + E. rectale + R. bromii) / (Alistipes + Bilophila + Bacteroides)];
# zeros replaced by half of the smallest non-zero value of the cohort, as for the MD-index (methods 10.3)
log_ratio_rep <- function(num, den, factor = 0.5) {
  log10(ifelse(num > 0, num, min(num[num > 0]) * factor) /
        ifelse(den > 0, den, min(den[den > 0]) * factor))
}

david_species <- c("Eubacterium_rectale", "Ruminococcus_bromii")
david <- imap_dfr(data, function(d, cohort) {
  g <- genus_of(cohort, colnames(d$X))
  stopifnot(david_species %in% colnames(d$X))
  tibble(sample_id = rownames(d$X),
         plant  = rowSums(d$X[, g %in% "Roseburia" | colnames(d$X) %in% david_species, drop = FALSE]),
         animal = rowSums(d$X[, g %in% c("Alistipes", "Bilophila", "Bacteroides"), drop = FALSE])) %>%
    mutate(david = log_ratio_rep(plant, animal),
           david_small = log_ratio_rep(plant, animal, factor = 0.05),
           cohort = cohort)
}) %>% left_join(taxa %>% select(sample_id, cohort, subject_id, diagnosis, site, antibiotic),
                 by = c("sample_id", "cohort"))

# Franzosa: Wilcoxon (main), lm with centre (check); HMP2: lmer (main), Wilcoxon on the median
# per person (check) (see methods 1.11)
test_ratio <- function(d, y, label) {
  d$y <- d[[y]]
  if (d$cohort[1] == "franzosa") {
    fit <- lm(y ~ diagnosis + site, data = d)
    v <- d$y; dx <- d$diagnosis
  } else {
    fit <- lmer(y ~ diagnosis + site + (1 | subject_id), data = d)
    pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(y), .groups = "drop")
    v <- pp$v; dx <- pp$diagnosis
  }
  s <- coef(summary(fit))
  tibble(cohort = d$cohort[1], measure = y, analysis = label, samples = nrow(d),
         persons = n_distinct(d$subject_id),
         zeros_plant = sum(d$plant == 0), zeros_animal = sum(d$animal == 0),
         median_control = median(v[dx == "control"]), median_cd = median(v[dx == "cd"]),
         p_wilcoxon = wilcox.test(v[dx == "cd"], v[dx == "control"], exact = FALSE)$p.value,
         model_estimate = s["diagnosiscd", "Estimate"], ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
         model_p = s["diagnosiscd", grep("^Pr", colnames(s))])
}

dv_fr <- david %>% filter(cohort == "franzosa")
dv_hm <- david %>% filter(cohort == "hmp2")
sq2a_david <- bind_rows(
  test_ratio(dv_fr, "david", "main"),
  test_ratio(dv_fr, "david_small", "smaller replacement of zeros"),
  test_ratio(dv_fr %>% filter(antibiotic == "No"), "david", "without antibiotics"),
  test_ratio(dv_hm, "david", "main"),
  test_ratio(dv_hm, "david_small", "smaller replacement of zeros"))
print(sq2a_david, width = Inf)
write_csv(sq2a_david, file.path(outdir, "sq2a_david.csv"))
write_csv(david, file.path(outdir, "sq2a_david_scores.csv"))

# --- SQ2b: MaAsLin2, main (see methods 11.3, 11.5) ---
run_maaslin <- function(d, label, random, features = c("F_prausnitzii", "Roseburia")) {
  run_dir <- file.path(outdir, "maaslin2", label)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  feat <- d %>% select(sample_id, all_of(features)) %>% column_to_rownames("sample_id") %>% as.data.frame()
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
    transmute(analysis = label, taxon = feature, value, samples = N, not_zero = N.not.0, coef, stderr,
              ci_low = coef - qnorm(0.975) * stderr, ci_high = coef + qnorm(0.975) * stderr,  # normal approximation
              pval, p_holm = p.adjust(pval, "holm"))
}

sq2b_maaslin <- bind_rows(
  run_maaslin(fr, "franzosa", NULL),
  run_maaslin(fr %>% filter(antibiotic == "No"), "franzosa_without_antibiotics", NULL),
  run_maaslin(hm, "hmp2", "subject_id"))
print(sq2b_maaslin, n = Inf, width = Inf)
write_csv(sq2b_maaslin, file.path(outdir, "sq2b_maaslin2.csv"))

# --- SQ2b: CLR models, check (see methods 11.4, 11.5) ---
# same CLR as in 05_beta_diversity.R (see methods 7.7)
clr_tf <- function(X) {
  pc <- min(X[X > 0]) / 2
  L <- log(X + pc)
  L - rowMeans(L)
}

# CLR parts after the 10 % genus filter; all species of `genus` become one part,
# because the hypothesis is about the genus
clr_parts <- function(d, cohort, genus) {
  X <- d$X[, passes_filter(d$X, cohort, "species", main_thr), drop = FALSE]
  in_genus <- genus_of(cohort, colnames(X)) %in% genus
  parts <- cbind(X[, !in_genus, drop = FALSE], rowSums(X[, in_genus, drop = FALSE]))
  colnames(parts)[ncol(parts)] <- genus
  clr_tf(parts / rowSums(parts))
}

clr_values <- imap_dfr(data, function(d, cohort) {
  Z <- clr_parts(d, cohort, "Roseburia")
  tibble(sample_id = rownames(Z), F_prausnitzii = Z[, fp_name], Roseburia = Z[, "Roseburia"], cohort = cohort)
}) %>% left_join(taxa %>% select(sample_id, subject_id, diagnosis, site, antibiotic), by = "sample_id")

# CD effect on CLR values: lm (Franzosa) or lmer (HMP2), Wilcoxon as check, Holm per cohort
clr_tests <- function(values, taxa_names) values %>%
  pivot_longer(all_of(taxa_names), names_to = "taxon", values_to = "clr") %>%
  filter(!is.na(clr)) %>%
  group_split(cohort, taxon) %>%
  map_dfr(function(d) {
    if (d$cohort[1] == "franzosa") {
      fit <- lm(clr ~ diagnosis + site, data = d); s <- coef(summary(fit))
      v <- d %>% transmute(diagnosis, v = clr)
    } else {
      fit <- lmer(clr ~ diagnosis + site + (1 | subject_id), data = d); s <- coef(summary(fit))
      v <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(clr), .groups = "drop")
    }
    tibble(cohort = d$cohort[1], taxon = d$taxon[1], n = nrow(d),
           median_control = median(v$v[v$diagnosis == "control"]),
           median_cd = median(v$v[v$diagnosis == "cd"]),
           p_wilcoxon = wilcox.test(v ~ diagnosis, data = v, exact = FALSE)$p.value,
           model_estimate = s["diagnosiscd", "Estimate"], model_se = s["diagnosiscd", "Std. Error"],
           ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
           model_p = s["diagnosiscd", grep("^Pr", colnames(s), value = TRUE)])
  }) %>%
  group_by(cohort) %>% mutate(model_p_holm = p.adjust(model_p, "holm")) %>% ungroup()

sq2b_clr <- clr_tests(clr_values, c("F_prausnitzii", "Roseburia"))
print(sq2b_clr, n = Inf, width = Inf)
write_csv(sq2b_clr, file.path(outdir, "sq2b_clr.csv"))

# --- SQ2b extension: other butyrate producers, exploratory (see methods 11.6) ---
# Present in both cohorts; Holm over these four taxa, separately from the main hypothesis
extra <- c("E_rectale", "A_hadrus", "E_hallii", "Coprococcus")
extra_species <- c(E_rectale = "Eubacterium_rectale", A_hadrus = "Anaerostipes_hadrus", E_hallii = "Eubacterium_hallii")

butyrate_extra <- imap_dfr(data, function(d, cohort) {
  g <- genus_of(cohort, colnames(d$X))
  stopifnot(extra_species %in% colnames(d$X))
  tibble(sample_id = rownames(d$X),
         E_rectale = d$X[, extra_species[["E_rectale"]]], A_hadrus = d$X[, extra_species[["A_hadrus"]]],
         E_hallii = d$X[, extra_species[["E_hallii"]]],
         Coprococcus = rowSums(d$X[, g %in% "Coprococcus", drop = FALSE]))
}) %>% left_join(taxa %>% select(sample_id, subject_id, diagnosis, site, antibiotic, cohort), by = "sample_id")

sq2b_extra_maaslin <- bind_rows(
  run_maaslin(filter(butyrate_extra, cohort == "franzosa"), "franzosa_extra", NULL, extra),
  run_maaslin(filter(butyrate_extra, cohort == "hmp2"), "hmp2_extra", "subject_id", extra))
print(sq2b_extra_maaslin, n = Inf, width = Inf)
write_csv(sq2b_extra_maaslin, file.path(outdir, "sq2b_extra_maaslin2.csv"))

# CLR check; a species removed by the genus filter stays NA and is not tested
clr_extra <- imap_dfr(data, function(d, cohort) {
  Z <- clr_parts(d, cohort, "Coprococcus")
  col <- function(s) if (s %in% colnames(Z)) Z[, s] else NA_real_
  tibble(sample_id = rownames(Z), E_rectale = col(extra_species[["E_rectale"]]),
         A_hadrus = col(extra_species[["A_hadrus"]]), E_hallii = col(extra_species[["E_hallii"]]),
         Coprococcus = Z[, "Coprococcus"], cohort = cohort)
}) %>% left_join(taxa %>% select(sample_id, subject_id, diagnosis, site, antibiotic), by = "sample_id")

sq2b_extra_clr <- clr_tests(clr_extra, extra)
print(sq2b_extra_clr, n = Inf, width = Inf)
write_csv(sq2b_extra_clr, file.path(outdir, "sq2b_extra_clr.csv"))

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

# H2a: both ratios, median per person; Prevotella/Bacteroides only in samples with both taxa
# (see methods 11.2, 11.7); p = main test (Franzosa Wilcoxon, HMP2 lmer, see methods 1.11)
ratio_labels <- c(pb = "Prevotella / Bacteroides\n(samples with both taxa)",
                  david = "Plant / animal taxa\n(David 2014, all samples)")
plot_h2a <- bind_rows(
  both %>% transmute(cohort, subject_id, diagnosis, measure = "pb", value = log_ratio),
  david %>% transmute(cohort, subject_id, diagnosis, measure = "david", value = david)) %>%
  group_by(cohort, measure, subject_id, diagnosis) %>%
  summarise(value = median(value), .groups = "drop")

main_p <- bind_rows(
  sq2a_ratio %>% transmute(cohort, measure = "pb", p = if_else(cohort == "franzosa", p_wilcoxon, model_p)),
  sq2a_david %>% filter(analysis == "main") %>%
    transmute(cohort, measure = "david", p = if_else(cohort == "franzosa", p_wilcoxon, model_p))) %>%
  mutate(label = paste0(if_else(cohort == "franzosa", "Wilcoxon ", "lmer "),
                        if_else(p < 0.001, "p < 0.001", paste("p =", signif(p, 2)))))

tidy_h2a <- function(d) d %>%
  mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
         measure = factor(measure, names(ratio_labels), ratio_labels))

# one plot for the given ratios: both ratios as cohort x ratio grid, one ratio with cohorts side by side
h2a_plot <- function(measures, title) {
  d <- tidy_h2a(filter(plot_h2a, measure %in% measures)) %>%
    mutate(diagnosis = factor(diagnosis, c("control", "cd"), c("Control", "CD")))
  ggplot(d, aes(diagnosis, value, colour = diagnosis)) +
    geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, width = .55) +
    geom_jitter(width = .15, alpha = .6, size = 1.2) +
    geom_text(data = tidy_h2a(filter(main_p, measure %in% measures)), aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.4, size = 3.6) +
    scale_y_continuous(expand = expansion(mult = c(.05, .15))) +   # room for the p label
    scale_colour_manual(values = setNames(unname(col_diag), c("Control", "CD")), guide = "none") +
    (if (length(measures) > 1) facet_grid(cohort ~ measure) else facet_wrap(~ cohort)) +
    labs(x = NULL, y = "log10 ratio (median per person)", title = title,
         caption = "HMP2 lmer adjusted for centre and person; above 0 = more of the first group") +
    theme_plot
}
ggsave(file.path(outdir, "sq2a_ratios_cd_vs_control.png"),
       h2a_plot(c("pb", "david"), "H2a: diet-related taxa ratios in controls and CD"), width = 8, height = 7, dpi = 200)
ggsave(file.path(outdir, "sq2a_pb_cd_vs_control.png"),
       h2a_plot("pb", "H2a: Prevotella / Bacteroides (samples with both taxa)"), width = 7, height = 4.5, dpi = 200)
ggsave(file.path(outdir, "sq2a_david_cd_vs_control.png"),
       h2a_plot("david", "H2a: plant / animal taxa ratio (David 2014)"), width = 7, height = 4.5, dpi = 200)

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

# extension: the four other butyrate producers
plot_extra <- butyrate_extra %>%
  group_by(cohort, subject_id, diagnosis) %>%
  summarise(across(all_of(extra), median), .groups = "drop") %>%
  pivot_longer(all_of(extra), names_to = "taxon", values_to = "abund") %>%
  mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
         taxon = factor(taxon, extra, c("E. rectale", "A. hadrus", "E. hallii", "Coprococcus")),
         diagnosis = factor(diagnosis, c("control", "cd"), c("Control", "CD")),
         pct = pmax(abund * 100, 1e-4))

p_extra <- ggplot(plot_extra, aes(diagnosis, pct, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = .55) +
  geom_jitter(width = .15, alpha = .6, size = 1.2) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = setNames(unname(col_diag), c("Control", "CD")), guide = "none") +
  facet_grid(cohort ~ taxon) +
  labs(x = NULL, y = "relative abundance (%, log scale)",
       title = "Other butyrate producers in controls and CD (median per person)",
       caption = "Persons without the taxon are shown at the bottom") +
  theme_plot
ggsave(file.path(outdir, "sq2b_extra_cd_vs_control.png"), p_extra, width = 10, height = 6.5, dpi = 200)
