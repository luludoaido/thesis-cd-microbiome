# --- 08 H1: Dysbiosis ---
# Compares CD and controls in two dysbiosis measures: the MD-index (Gevers) and the median
# Bray-Curtis distance to the controls (Lloyd-Price). Franzosa and HMP2 only.

# Results:
# - MD-index and Bray-Curtis score per sample
# - MD-index in four versions: Pascal 2017 stool signature and Gevers 2014 restricted to taxa
#   found in stool (both main), all 13 taxa of Gevers Table S2I and the 7 of the summary (checks)
# - CD effect on both scores: Franzosa Wilcoxon (main), lm with centre as check; HMP2 lmer (main),
#   Wilcoxon on the median per person as check; estimate, SE, 95 % CI, p; with sensitivity analyses
# - share of dysbiotic samples (above the 90th percentile of the controls) in CD and controls
#   (Franzosa: Fisher; HMP2: logistic mixed model, Wilcoxon per person as check)
# - detection of the Gevers taxa; agreement of the two measures

# Methods: docs/methods.md, 1.10 and section 10

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
source("R/00_setup.R")

outdir <- file.path(dir_results, "hypotheses")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("franzosa", "hmp2")
main_thr <- 0.10

# MD-index taxa (see methods 10.2): Gevers 2014 Table S2I (Greengenes names mapped to the
# taxonomy used here; short = the 7 taxa of the summary) and the stool signature of Pascal 2017
md_taxa <- tribble(
  ~set,     ~taxon,                  ~direction, ~rank,    ~match,                            ~short,
  "gevers", "Enterobacteriaceae",    "up",       "family", "Enterobacteriaceae",              TRUE,
  "gevers", "Pasteurellaceae",       "up",       "family", "Pasteurellaceae",                 TRUE,
  "gevers", "Veillonellaceae",       "up",       "family", "Veillonellaceae",                 TRUE,
  "gevers", "Fusobacteriaceae",      "up",       "family", "Fusobacteriaceae",                TRUE,
  "gevers", "Neisseriaceae",         "up",       "family", "Neisseriaceae",                   FALSE,
  "gevers", "Micrococcaceae",        "up",       "family", "Micrococcaceae",                  FALSE,
  "gevers", "Gemellaceae",           "up",       "genus",  "Gemella",                         FALSE,
  "gevers", "Bacteroidales",         "down",     "order",  "Bacteroidales",                   TRUE,
  "gevers", "Clostridiales",         "down",     "order",  "Clostridiales",                   TRUE,
  "gevers", "Erysipelotrichaceae",   "down",     "order",  "Erysipelotrichales",              TRUE,
  "gevers", "Bifidobacteriaceae",    "down",     "family", "Bifidobacteriaceae",              FALSE,
  "gevers", "Verrucomicrobiaceae",   "down",     "family", "Akkermansiaceae",                 FALSE,
  "gevers", "Coriobacteriaceae",     "down",     "order",  "Coriobacteriales|Eggerthellales", FALSE,
  "pascal", "Escherichia",           "up",       "genus",  "Escherichia",                     FALSE,
  "pascal", "Fusobacterium",         "up",       "genus",  "Fusobacterium",                   FALSE,
  "pascal", "Faecalibacterium",      "down",     "genus",  "Faecalibacterium",                FALSE,
  "pascal", "Anaerostipes",          "down",     "genus",  "Anaerostipes",                    FALSE,
  "pascal", "Collinsella",           "down",     "genus",  "Collinsella",                     FALSE,
  "pascal", "Methanobrevibacter",    "down",     "genus",  "Methanobrevibacter",              FALSE,
  "pascal", "Christensenellaceae",   "down",     "family", "Christensenellaceae",             FALSE,
  "pascal", "Peptostreptococcaceae", "down",     "family", "Peptostreptococcaceae",           FALSE) %>%
  mutate(id = paste(set, taxon, sep = ":"))

min_prev <- 0.10   # Gevers taxa must be found in >= 10 % of all samples of every cohort (see methods 10.2)

# --- Taxonomy (see methods 10.7) ---
hmp2_tax <- read_csv(file.path(dir_processed, "hmp2/hmp2_taxonomy.csv"), show_col_types = FALSE)

# Franzosa "genera" that are really higher levels, or that HMP2 does not have
manual_lookup <- tribble(
  ~genus,                  ~family,                 ~order,
  "Bacteroidales",         NA,                      "Bacteroidales",
  "Clostridiales",         NA,                      "Clostridiales",
  "Clostridiaceae",        "Clostridiaceae",        "Clostridiales",
  "Lachnospiraceae",       "Lachnospiraceae",       "Clostridiales",
  "Peptostreptococcaceae", "Peptostreptococcaceae", "Clostridiales",
  "Ruminococcaceae",       "Ruminococcaceae",       "Clostridiales",
  "Subdoligranulum",       "Ruminococcaceae",       "Clostridiales",
  "Erysipelotrichaceae",   "Erysipelotrichaceae",   "Erysipelotrichales",
  "Burkholderiales",       NA,                      "Burkholderiales",
  "Methanosphaera",        "Methanobacteriaceae",   "Methanobacteriales",
  "Propionibacterium",     "Propionibacteriaceae",  "Propionibacteriales")

genus_lookup <- bind_rows(manual_lookup,
                          hmp2_tax %>% filter(!genus %in% manual_lookup$genus) %>%
                            count(genus, family, order) %>% group_by(genus) %>%
                            slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>% select(-n))

# check: genera with more than one family in HMP2, which could be misassigned in Franzosa
ambiguous <- hmp2_tax %>% distinct(genus, family) %>% count(genus) %>% filter(n > 1)
cat("Genera with more than one family in the HMP2 taxonomy:", nrow(ambiguous), "\n")
print(ambiguous)

taxonomy_of <- function(cohort, species) {
  if (cohort == "hmp2") {
    hmp2_tax[match(species, hmp2_tax$species), c("family", "order")]
  } else {
    genus_lookup[match(sub("_.*", "", species), genus_lookup$genus), c("family", "order")]
  }
}

# --- Load ---
load_cohort <- function(cohort) {
  folder <- file.path(dir_processed, cohort)
  meta <- read_csv(file.path(folder, paste0(cohort, "_meta.csv")), show_col_types = FALSE) %>%
    mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))
  read_X <- function(level) {
    d <- read_csv(file.path(folder, paste0(cohort, "_", level, ".csv")), show_col_types = FALSE) %>%
      filter(sample_id %in% meta$sample_id) %>%
      arrange(match(sample_id, meta$sample_id)) %>%
      column_to_rownames("sample_id") %>% as.matrix()
    d / rowSums(d)
  }
  list(meta = meta, species = read_X("species"), genus = read_X("genus"))
}
data <- set_names(map(cohorts, load_cohort), cohorts)

# --- MD-index (see methods 10.1-10.4) ---
# abundance of each MD-index taxon per sample
md_parts <- function(cohort, X) {
  tax <- taxonomy_of(cohort, colnames(X)) %>% mutate(genus = sub("_.*", "", colnames(X)))
  cat(cohort, ": mean share of abundance without an order:",
      round(mean(rowSums(X[, is.na(tax$order), drop = FALSE])) * 100, 2), "%\n")
  parts <- map(set_names(md_taxa$id), function(i) {
    g <- md_taxa[md_taxa$id == i, ]
    rowSums(X[, tax[[g$rank]] %in% strsplit(g$match, "|", fixed = TRUE)[[1]], drop = FALSE])
  })
  tibble(sample_id = rownames(X), !!!parts)
}
parts <- imap(data, ~ md_parts(.y, .x$species))

# found in stool: share of samples > 0 in every cohort, all diagnoses together (see methods 10.2)
stool_ok <- map_dfr(parts, ~ summarise(.x, across(all_of(md_taxa$id), ~ mean(.x > 0))), .id = "cohort") %>%
  pivot_longer(-cohort, names_to = "id", values_to = "prev") %>%
  group_by(id) %>% summarise(min_prev_pct = min(prev) * 100, stool = min(prev) >= min_prev)
md_taxa <- md_taxa %>% left_join(stool_ok, by = "id")
print(md_taxa %>% select(set, taxon, direction, min_prev_pct, stool), n = Inf)

# Pascal is a stool signature and is used as published; Gevers only with taxa found in stool
md_versions <- list(
  md_pascal       = md_taxa$set == "pascal",
  md_gevers_stool = md_taxa$set == "gevers" & md_taxa$stool,
  md_gevers13     = md_taxa$set == "gevers",
  md_gevers7      = md_taxa$set == "gevers" & md_taxa$short)
md_main <- c("md_pascal", "md_gevers_stool")

md_of <- function(p, sel, factor = 0.5) {
  sum_of <- function(dir) rowSums(as.matrix(p[, md_taxa$id[sel & md_taxa$direction == dir], drop = FALSE]))
  md_index(sum_of("up"), sum_of("down"), factor)
}

md_index <- function(up, down, factor = 0.5) {
  rep_up <- min(up[up > 0]) * factor
  rep_down <- min(down[down > 0]) * factor
  log10(ifelse(up > 0, up, rep_up) / ifelse(down > 0, down, rep_down))
}

# --- Distance-based score (see methods 10.5) ---
# same CLR as in 05_beta_diversity.R (see methods 7.7)
clr_tf <- function(X) {
  pc <- min(X[X > 0]) / 2
  L <- log(X + pc)
  L - rowMeans(L)
}

# median distance of each sample to the reference controls, without samples of the same person
distance_score <- function(cohort, d, level, distance = "bray", thr = 0, after_week = NULL) {
  X <- d[[level]]
  meta <- d$meta
  keep <- if (thr == 0) colSums(X > 0) > 0 else passes_filter(X, cohort, level, thr)
  Xf <- X[, keep, drop = FALSE]
  Xf <- Xf / rowSums(Xf)
  D <- if (distance == "bray") as.matrix(vegdist(Xf, method = "bray")) else as.matrix(dist(clr_tf(Xf)))
  in_ref <- meta$diagnosis == "control"
  if (!is.null(after_week)) in_ref <- in_ref & meta$time > after_week
  sapply(seq_len(nrow(D)), function(i) {
    median(D[i, in_ref & meta$subject_id != meta$subject_id[i]])
  })
}

# "after_week" only applies to HMP2, Franzosa has one time point per person
score_variants <- function(cohort, d) {
  wk <- if (cohort == "hmp2") 20 else NULL
  tibble(bray            = distance_score(cohort, d, "species", thr = main_thr, after_week = wk),
         bray_all_ctrl   = distance_score(cohort, d, "species", thr = main_thr),
         bray_no_filter  = distance_score(cohort, d, "species", after_week = wk),
         bray_genus10    = distance_score(cohort, d, "genus", thr = main_thr, after_week = wk),
         aitch_species10 = distance_score(cohort, d, "species", "aitchison", main_thr, wk),
         aitch_genus10   = distance_score(cohort, d, "genus", "aitchison", main_thr, wk))
}
dist_measures <- c("bray", "bray_all_ctrl", "bray_no_filter", "bray_genus10", "aitch_species10", "aitch_genus10")

scores <- imap_dfr(data, function(d, cohort) {
  p <- parts[[cohort]]
  p %>%
    mutate(!!!map(md_versions, ~ md_of(p, .x)),
           md_pascal_small       = md_of(p, md_versions$md_pascal, 0.05),
           md_gevers_stool_small = md_of(p, md_versions$md_gevers_stool, 0.05)) %>%
    left_join(d$meta %>% select(sample_id, subject_id, diagnosis, site, time, baseline, antibiotic),
              by = "sample_id") %>%
    bind_cols(score_variants(cohort, d)) %>%
    mutate(cohort = cohort)
})
write_csv(scores, file.path(outdir, "dysbiosis_scores.csv"))

# --- Detection of the MD-index taxa ---
detection <- scores %>%
  pivot_longer(all_of(md_taxa$id), names_to = "id", values_to = "abund") %>%
  group_by(cohort, diagnosis, id) %>%
  summarise(samples = n(), prevalence_pct = mean(abund > 0) * 100,
            median_pct = median(abund) * 100, .groups = "drop") %>%
  left_join(md_taxa %>% select(id, set, taxon, direction, stool), by = "id") %>%
  mutate(id = factor(id, md_taxa$id)) %>%
  arrange(cohort, id, diagnosis) %>%
  select(cohort, set, taxon, direction, stool, diagnosis, samples, prevalence_pct, median_pct)
print(detection, n = Inf)
write_csv(detection, file.path(outdir, "dysbiosis_taxa_detection.csv"))

# --- Tests (see methods 10.8) ---
median_diff <- function(y, dx) median(y[dx == "cd"]) - median(y[dx == "control"])

test_franzosa <- function(d, y, label) {
  d$y <- d[[y]]
  fit <- lm(y ~ diagnosis + site, data = d)
  tibble(cohort = "franzosa", measure = y, analysis = label, samples = nrow(d),
         persons = n_distinct(d$subject_id),
         median_control = median(d$y[d$diagnosis == "control"]),
         median_cd = median(d$y[d$diagnosis == "cd"]),
         p_wilcoxon = wilcox.test(y ~ diagnosis, data = d, exact = FALSE)$p.value,
         model_estimate = coef(summary(fit))["diagnosiscd", "Estimate"],
         model_se = coef(summary(fit))["diagnosiscd", "Std. Error"],
         ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
         model_p = coef(summary(fit))["diagnosiscd", "Pr(>|t|)"])
}

# HMP2: main test lmer on all samples; check: Wilcoxon on the median per person (see methods 1.5)
test_hmp2 <- function(d, y, label) {
  d$y <- d[[y]]
  pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = median(y), .groups = "drop")
  fit <- lmer(y ~ diagnosis + site + (1 | subject_id), data = d)
  tibble(cohort = "hmp2", measure = y, analysis = label, samples = nrow(d),
         persons = n_distinct(d$subject_id),
         median_control = median(pp$v[pp$diagnosis == "control"]),
         median_cd = median(pp$v[pp$diagnosis == "cd"]),
         p_wilcoxon = wilcox.test(v ~ diagnosis, data = pp, exact = FALSE)$p.value,
         model_estimate = coef(summary(fit))["diagnosiscd", "Estimate"],
         model_se = coef(summary(fit))["diagnosiscd", "Std. Error"],
         ci_low = ci_cd(fit)[1], ci_high = ci_cd(fit)[2],
         model_p = coef(summary(fit))["diagnosiscd", "Pr(>|t|)"])
}

fr <- scores %>% filter(cohort == "franzosa")
hm <- scores %>% filter(cohort == "hmp2")

tests <- bind_rows(
  map_dfr(md_main, ~ test_franzosa(fr, .x, "main")),
  map_dfr(md_main, ~ test_franzosa(fr, paste0(.x, "_small"), "smaller replacement of zeros")),
  map_dfr(md_main, ~ test_franzosa(fr %>% filter(antibiotic == "No"), .x, "without antibiotics")),
  test_franzosa(fr, "md_gevers13", "check: all 13 Gevers taxa"),
  test_franzosa(fr, "md_gevers7", "check: 7 taxa of the Gevers summary"),
  map_dfr(dist_measures, ~ test_franzosa(fr, .x, "main")),
  test_franzosa(fr %>% filter(antibiotic == "No"), "bray", "without antibiotics"),
  map_dfr(md_main, ~ test_hmp2(hm, .x, "main")),
  map_dfr(md_main, ~ test_hmp2(hm, paste0(.x, "_small"), "smaller replacement of zeros")),
  test_hmp2(hm, "md_gevers13", "check: all 13 Gevers taxa"),
  test_hmp2(hm, "md_gevers7", "check: 7 taxa of the Gevers summary"),
  map_dfr(dist_measures, ~ test_hmp2(hm, .x, "main")))
print(tests, n = Inf, width = Inf)
write_csv(tests, file.path(outdir, "dysbiosis_tests.csv"))

# --- Dysbiotic samples (see methods 10.6) ---
dysbiotic <- bind_rows(
  tidyr::expand_grid(q = c(0.80, 0.90, 0.95), m = "bray"),
  tidyr::expand_grid(q = 0.90, m = setdiff(dist_measures, "bray"))
) %>% pmap_dfr(function(q, m) {
  map_dfr(c("franzosa", "hmp2"), function(ch) {
    d <- scores %>% filter(cohort == ch)
    limit <- quantile(d[[m]][d$diagnosis == "control"], q)
    d <- d %>% mutate(dysbiotic = .data[[m]] > limit)
    if (ch == "franzosa") {
      p <- fisher.test(table(d$diagnosis, d$dysbiotic))$p.value
      p_check <- NA_real_
      unit <- "samples (Fisher)"
      s <- d %>% group_by(diagnosis) %>% summarise(v = mean(dysbiotic) * 100)
    } else {
      # several samples per person: logistic mixed model on all samples (see methods 1.5);
      # check: Wilcoxon on the share of dysbiotic time points per person
      fit <- glmer(dysbiotic ~ diagnosis + site + (1 | subject_id), data = d, family = binomial,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      p <- coef(summary(fit))["diagnosiscd", "Pr(>|z|)"]
      pp <- d %>% group_by(subject_id, diagnosis) %>% summarise(v = mean(dysbiotic) * 100, .groups = "drop")
      p_check <- wilcox.test(v ~ diagnosis, data = pp, exact = FALSE)$p.value
      unit <- "logistic mixed model; check: Wilcoxon on share per person"
      s <- pp %>% group_by(diagnosis) %>% summarise(v = mean(v))
    }
    tibble(cohort = ch, measure = m, percentile = q * 100, limit = unname(limit), unit = unit,
           control_pct = s$v[s$diagnosis == "control"], cd_pct = s$v[s$diagnosis == "cd"],
           p = p, p_check = p_check)
  })
})
print(dysbiotic, n = Inf, width = Inf)
write_csv(dysbiotic, file.path(outdir, "dysbiosis_dysbiotic_share.csv"))

# how similar are the two measures?
print(scores %>% group_by(cohort) %>%
        summarise(rho_pascal_bray = cor(md_pascal, bray, method = "spearman"),
                  rho_gevers_stool_bray = cor(md_gevers_stool, bray, method = "spearman"),
                  rho_pascal_gevers_stool = cor(md_pascal, md_gevers_stool, method = "spearman")))

# --- Plot ---
plot_measures <- c("MD-index Pascal (log10)", "MD-index Gevers stool (log10)", "Bray-Curtis score")
p <- scores %>%
  group_by(cohort, subject_id, diagnosis) %>%
  summarise(`MD-index Pascal (log10)` = median(md_pascal),
            `MD-index Gevers stool (log10)` = median(md_gevers_stool),
            `Bray-Curtis score` = median(bray), .groups = "drop") %>%
  pivot_longer(all_of(plot_measures), names_to = "measure", values_to = "value") %>%
  mutate(cohort = recode(cohort, franzosa = "Franzosa", hmp2 = "HMP2"),
         measure = factor(measure, plot_measures)) %>%
  ggplot(aes(diagnosis, value, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = .55) +
  geom_jitter(width = .15, alpha = .6, size = 1.3) +
  scale_colour_manual(values = col_diag, guide = "none") +
  facet_grid(measure ~ cohort, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "Dysbiosis in CD and controls (median per person)") +
  theme_plot
ggsave(file.path(outdir, "dysbiosis_cd_vs_control.png"), p, width = 7, height = 9, dpi = 200)
