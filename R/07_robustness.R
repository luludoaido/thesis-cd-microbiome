library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
source("R/00_setup.R")

# =====================================================================
# Robustness checks for the main analysis (06)
# 1. Clinical scores: calprotectin + HBI for all samples (plot)
# 2. Franzosa: which samples are the Shannon outliers? (mean -0.37 vs. median 0.30)
# 3. Does sample size affect Shannon or its SD? (reads, samples per person)
# 4. Time trend HMP2: why does Shannon rise only in controls? Model stratified by time
# 5. Hall: is "Boston only" usable?
# reads the finished results of R/06_diversity.R (genus, abundance_results.rds)
# main analysis: 10 % filter (thr_10)
# =====================================================================

res_dir <- dir_results
path    <- dir_processed
outdir  <- file.path(res_dir, "robustness")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# colours col_diag and theme_plot come from R/00_setup.R

res <- map(set_names(c("franzosa", "hmp2", "hall")), ~ readRDS(file.path(res_dir, .x, "abundance_results.rds")))
main_thr <- "thr_10"
alpha <- map(res, ~ .x$alpha %>% filter(threshold == main_thr) %>%
               mutate(evenness = Shannon / log(Richness)))   # Pielou J: evenness only

# small helper: coefficient + p from a (g)lmer/lm as a table row
coef_row <- function(m, term) {
  s <- coef(summary(m))
  p_col <- grep("^Pr", colnames(s))
  tibble(term = term, estimate = s[term, "Estimate"], se = s[term, "Std. Error"], p = s[term, p_col])
}

# ---- 1. Clinical scores: calprotectin + HBI for all samples ------------
# categories:
#   calprotectin < 50 normal | 50-250 grey zone | > 250 active inflammation (µg/g)
#   HBI < 5 remission | 5-7 mild | 8-16 moderate | > 16 severe (only meaningful for CD)
calpro_cat <- function(x) cut(x, c(-Inf, 50, 250, Inf), right = FALSE,
                              labels = c("< 50 normal", "50-250 grey zone", "> 250 active"))
hbi_cat <- function(x) cut(x, c(-Inf, 5, 8, 17, Inf), right = FALSE,
                           labels = c("< 5 remission", "5-7 mild", "8-16 moderate", "> 16 severe"))

# only samples that are in the analysis (same exclusions)
fr_meta <- read.csv(file.path(path, "franzosa/franzosa_meta.csv"))
ha_meta <- read.csv(file.path(path, "hall/hall_meta.csv"))
hm_clin <- read_csv(file.path(res_dir, "hmp2_clinical/hmp2_clinical_per_sample.csv"), show_col_types = FALSE)

clinical <- bind_rows(
  alpha$franzosa %>% select(sample_id, subject_id, diagnosis) %>%
    left_join(fr_meta %>% transmute(sample_id, calpro = Fecal.Calprotectin, hbi = NA_real_), by = "sample_id") %>%
    mutate(cohort = "Franzosa"),
  alpha$hmp2 %>% select(sample_id, subject_id, diagnosis) %>%
    left_join(hm_clin %>% select(sample_id, calpro, hbi), by = "sample_id") %>%
    mutate(cohort = "HMP2"),
  alpha$hall %>% select(sample_id, subject_id, diagnosis, site) %>%
    left_join(ha_meta %>% transmute(sample_id = ID, calpro = calprotectin, hbi = HBI), by = "sample_id") %>%
    mutate(cohort = "Hall")
) %>%
  mutate(cohort = factor(cohort, c("Franzosa", "HMP2", "Hall")),
         calpro_cat = calpro_cat(calpro), hbi_cat = hbi_cat(hbi))

clinical_tab <- clinical %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(), persons = n_distinct(subject_id),
            calpro_n = sum(!is.na(calpro)), calpro_median = median(calpro, na.rm = TRUE),
            calpro_lt50 = sum(calpro_cat == levels(calpro_cat)[1], na.rm = TRUE),
            calpro_50_250 = sum(calpro_cat == levels(calpro_cat)[2], na.rm = TRUE),
            calpro_gt250 = sum(calpro_cat == levels(calpro_cat)[3], na.rm = TRUE),
            hbi_n = sum(!is.na(hbi)), hbi_median = median(hbi, na.rm = TRUE),
            hbi_remission = sum(hbi_cat == levels(hbi_cat)[1], na.rm = TRUE),
            hbi_mild = sum(hbi_cat == levels(hbi_cat)[2], na.rm = TRUE),
            hbi_moderate = sum(hbi_cat == levels(hbi_cat)[3], na.rm = TRUE),
            hbi_severe = sum(hbi_cat == levels(hbi_cat)[4], na.rm = TRUE), .groups = "drop")
clinical_tab
write_csv(clinical_tab, file.path(outdir, "clinical_categories.csv"))

# HBI values in controls would be a data error
clinical %>% filter(diagnosis == "control", !is.na(hbi))

# plot A: each point = one sample, background = category
# (log scale; values < 5 are set to 5 so that 0 can be shown)
bands_cal <- tibble(ymin = c(5, 50, 250), ymax = c(50, 250, 5000),
                    cat = factor(levels(clinical$calpro_cat), levels(clinical$calpro_cat)))
# light versions of the severity colours as background bands
col_cal <- adjustcolor(col_severity[c(1, 3, 4)], alpha.f = .25) %>% set_names(levels(clinical$calpro_cat))
d <- clinical %>% filter(!is.na(calpro)) %>% mutate(calpro_plot = pmin(pmax(calpro, 5), 5000))
p_cal <- ggplot(d, aes(diagnosis, calpro_plot)) +
  geom_rect(data = bands_cal, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = cat),
            inherit.aes = FALSE) +
  geom_hline(yintercept = c(50, 250), linetype = "dashed", colour = okabe_ito[["grey"]]) +
  geom_jitter(aes(colour = diagnosis), width = .25, height = 0, size = 1.4, alpha = .6) +
  geom_boxplot(width = .35, outlier.shape = NA, fill = NA, colour = "black") +
  scale_y_log10(breaks = c(5, 10, 50, 100, 250, 1000, 5000)) +
  scale_fill_manual(values = col_cal) + scale_colour_manual(values = col_diag, guide = "none") +
  facet_wrap(~ cohort) +
  labs(x = NULL, y = "faecal calprotectin (µg/g, log10)", fill = NULL,
       title = "Calprotectin - every sample") + theme_plot
p_cal
ggsave(file.path(outdir, "clinical_calprotectin_all.png"), p_cal, width = 8.5, height = 4.8, dpi = 200)

bands_hbi <- tibble(ymin = c(-.5, 4.5, 7.5, 16.5), ymax = c(4.5, 7.5, 16.5, 30),
                    cat = factor(levels(clinical$hbi_cat), levels(clinical$hbi_cat)))
col_hbi <- adjustcolor(col_severity, alpha.f = .25) %>% set_names(levels(clinical$hbi_cat))
d <- clinical %>% filter(!is.na(hbi), diagnosis == "cd")
p_hbi <- ggplot(d, aes(cohort, hbi)) +
  geom_rect(data = bands_hbi, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = cat),
            inherit.aes = FALSE) +
  geom_jitter(colour = col_diag["cd"], width = .25, height = .15, size = 1.4, alpha = .6) +
  geom_boxplot(width = .35, outlier.shape = NA, fill = NA, colour = "black") +
  scale_fill_manual(values = col_hbi) +
  coord_cartesian(ylim = c(0, max(d$hbi) + 1)) +
  labs(x = NULL, y = "HBI (CD samples only)", fill = NULL,
       title = "HBI - every CD sample") + theme_plot
p_hbi
ggsave(file.path(outdir, "clinical_hbi_all.png"), p_hbi, width = 6, height = 4.8, dpi = 200)

# plot B: shares per category (samples), n = number of samples with value
shares <- bind_rows(
  clinical %>% filter(!is.na(calpro)) %>% count(cohort, diagnosis, cat = calpro_cat) %>% mutate(score = "Calprotectin"),
  clinical %>% filter(!is.na(hbi), diagnosis == "cd") %>% count(cohort, diagnosis, cat = hbi_cat) %>% mutate(score = "HBI (CD)")
) %>%
  group_by(score, cohort, diagnosis) %>% mutate(share = n / sum(n), n_tot = sum(n)) %>% ungroup() %>%
  mutate(group = paste0(diagnosis, "\nn=", n_tot),
         cat = factor(cat, c(levels(clinical$calpro_cat), levels(clinical$hbi_cat))))
p_share <- ggplot(shares, aes(group, share, fill = cat)) +
  geom_col(width = .7, colour = "white") +
  geom_text(aes(label = ifelse(share >= .06, n, "")), position = position_stack(vjust = .5), size = 3.3) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(col_severity[c(1, 3, 4)], col_severity) %>%
                      set_names(levels(shares$cat))) +
  facet_grid(~ score + cohort, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "share of samples", fill = NULL) + theme_plot +
  guides(fill = guide_legend(ncol = 4))
p_share
ggsave(file.path(outdir, "clinical_categories_shares.png"), p_share, width = 10, height = 5.2, dpi = 200)

# plot C: HMP2 + Hall per person over time (how does activity change?)
d <- clinical %>% filter(cohort != "Franzosa", diagnosis == "cd", !is.na(hbi)) %>%
  left_join(bind_rows(alpha$hmp2, alpha$hall) %>% select(sample_id, time), by = "sample_id")
p_hbi_time <- ggplot(d, aes(time, hbi, group = subject_id)) +
  geom_hline(yintercept = c(4.5, 7.5, 16.5), linetype = "dashed", colour = okabe_ito[["grey"]]) +
  geom_line(colour = col_diag["cd"], alpha = .35) + geom_point(colour = col_diag["cd"], size = 1, alpha = .6) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(x = "weeks since first sample", y = "HBI", title = "HBI per CD person over time") + theme_plot
ggsave(file.path(outdir, "clinical_hbi_time.png"), p_hbi_time, width = 8.5, height = 4.2, dpi = 200)

# ---- 2. Franzosa: which samples are the outliers? --------------------
fr <- alpha$franzosa %>%
  left_join(fr_meta %>% select(sample_id, calpro = Fecal.Calprotectin, immunosuppressant, mesalamine, steroids),
            by = "sample_id") %>%
  group_by(diagnosis) %>%
  mutate(q1 = quantile(Shannon, .25), q3 = quantile(Shannon, .75),
         outlier_low = Shannon < q1 - 1.5 * (q3 - q1),      # boxplot rule (whiskers)
         outlier_high = Shannon > q3 + 1.5 * (q3 - q1),
         rank_low = rank(Shannon)) %>%
  ungroup()

fr_out <- fr %>%
  filter(outlier_low | outlier_high | (diagnosis == "cd" & rank_low <= 10)) %>%
  arrange(Shannon) %>%
  select(sample_id, diagnosis, site, Shannon, Richness, evenness, outlier_low, outlier_high,
         antibiotic, steroids, immunosuppressant, calpro, age)
print(fr_out, n = Inf)
write_csv(fr_out, file.path(outdir, "franzosa_shannon_outliers.csv"))

# are the outliers also conspicuous in the PCA/ordinations?
res$franzosa$outliers %>% filter(sample_id %in% fr_out$sample_id)

# who are they compared with the other CD samples? (antibiotics, calprotectin, richness)
fr %>% filter(diagnosis == "cd") %>%
  group_by(outlier_low) %>%
  summarise(n = n(), antibiotics = sum(antibiotic == "Yes", na.rm = TRUE),
            steroids_yes = sum(steroids == "Yes", na.rm = TRUE),
            calpro_median = median(calpro, na.rm = TRUE), richness_median = median(Richness),
            evenness_median = median(evenness))

# mean vs. median: how much do the outliers pull the effect?
fr_effects <- bind_rows(
  tibble(analysis = "difference of medians",
         estimate = median(fr$Shannon[fr$diagnosis == "cd"]) - median(fr$Shannon[fr$diagnosis == "control"]), p = NA),
  tibble(analysis = "Hodges-Lehmann (Wilcoxon)",
         estimate = wilcox.test(Shannon ~ diagnosis, data = fr, conf.int = TRUE)$estimate * -1,
         p = wilcox.test(Shannon ~ diagnosis, data = fr)$p.value),
  coef_row(lm(Shannon ~ diagnosis + site, data = fr), "diagnosiscd") %>%
    transmute(analysis = "lm + sub-cohort", estimate, p),
  coef_row(lm(Shannon ~ diagnosis + site, data = filter(fr, !outlier_low)), "diagnosiscd") %>%
    transmute(analysis = "lm + sub-cohort, without low outliers", estimate, p),
  coef_row(lm(Shannon ~ diagnosis + site, data = filter(fr, antibiotic != "Yes")), "diagnosiscd") %>%
    transmute(analysis = "lm + sub-cohort, without antibiotics", estimate, p)
)
fr_effects
write_csv(fr_effects, file.path(outdir, "franzosa_effect_with_without_outliers.csv"))

p_fr <- ggplot(fr, aes(diagnosis, Shannon, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = .5, colour = "black") +
  geom_jitter(aes(shape = antibiotic == "Yes"), width = .2, size = 1.8, alpha = .6) +
  ggrepel::geom_text_repel(data = filter(fr, outlier_low | outlier_high),
                           aes(label = sub(".*[|]", "", sample_id)), size = 3, colour = "black") +
  scale_colour_manual(values = col_diag, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c("no antibiotics", "antibiotics")) +
  facet_wrap(~ site) +
  labs(x = NULL, y = "Shannon (genus)", shape = NULL, title = "Franzosa: who are the outliers?") + theme_plot
p_fr
ggsave(file.path(outdir, "franzosa_outliers.png"), p_fr, width = 7.5, height = 4.8, dpi = 200)

# ---- 3. Does sample size affect Shannon / SD? ---------------------------
# 3a. sequencing depth (reads) -> Shannon? More reads = more rare genera detected.
#     Franzosa has no read counts -> HMP2 + Hall only
long <- bind_rows(alpha$hmp2 %>% mutate(cohort = "HMP2"), alpha$hall %>% mutate(cohort = "Hall")) %>%
  mutate(log_reads = log10(reads))

reads_tab <- long %>%
  group_by(cohort) %>%
  group_modify(function(d, key) {
    m0 <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d)
    m1 <- lmer(Shannon ~ diagnosis + log_reads + (1 | subject_id), data = d)
    mr <- lmer(Richness ~ log_reads + (1 | subject_id), data = d)
    mreads <- lmer(log_reads ~ diagnosis + (1 | subject_id), data = d)
    tibble(rho_shannon_reads = cor(d$Shannon, d$log_reads, method = "spearman"),
           rho_richness_reads = cor(d$Richness, d$log_reads, method = "spearman"),
           shannon_per_log10reads = coef_row(m1, "log_reads")$estimate,
           p_shannon_reads = coef_row(m1, "log_reads")$p,
           richness_per_log10reads = coef_row(mr, "log_reads")$estimate,
           p_richness_reads = coef_row(mr, "log_reads")$p,
           cd_effect_without_reads = coef_row(m0, "diagnosiscd")$estimate,
           cd_effect_with_reads = coef_row(m1, "diagnosiscd")$estimate,
           p_cd_with_reads = coef_row(m1, "diagnosiscd")$p,
           reads_median_cd = median(d$reads[d$diagnosis == "cd"]),
           reads_median_control = median(d$reads[d$diagnosis == "control"]),
           p_reads_cd_vs_control = coef_row(mreads, "diagnosiscd")$p)
  }) %>% ungroup()
glimpse(reads_tab)
write_csv(reads_tab, file.path(outdir, "reads_vs_diversity.csv"))

p_reads <- long %>%
  pivot_longer(c(Shannon, Richness), names_to = "metric") %>%
  ggplot(aes(reads / 1e6, value, colour = diagnosis)) +
  geom_point(size = 1.2, alpha = .45) +
  geom_smooth(method = "loess", se = FALSE, linewidth = .9) +
  scale_x_log10() + scale_colour_manual(values = col_diag) +
  facet_grid(metric ~ cohort, scales = "free") +
  labs(x = "reads per sample (millions, log10)", y = NULL, colour = NULL,
       title = "Does sequencing depth drive diversity?") + theme_plot
ggsave(file.path(outdir, "reads_vs_diversity.png"), p_reads, width = 8.5, height = 6, dpi = 200)

# 3b. samples per person -> SD of Shannon per person?
#     The sample SD is systematically too small for small n: E[s] = c4(n) * sigma
#     (c4(2) = 0.80, c4(5) = 0.94, c4(10) = 0.97) and very uncertain.
#     If CD and controls have different numbers of samples per person, this can bias the SD comparison.
c4 <- function(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))

sd_person <- long %>%
  group_by(cohort, subject_id, diagnosis) %>%
  summarise(n_samples = n(), sd_shannon = sd(Shannon), .groups = "drop") %>%
  filter(n_samples >= 2) %>%
  mutate(sd_corrected = sd_shannon / c4(n_samples))

sd_tab <- sd_person %>%
  group_by(cohort) %>%
  summarise(n_persons = n(),
            samples_median_cd = median(n_samples[diagnosis == "cd"]),
            samples_median_control = median(n_samples[diagnosis == "control"]),
            p_samples_cd_vs_control = wilcox.test(n_samples ~ diagnosis)$p.value,
            rho_n_vs_sd = cor(n_samples, sd_shannon, method = "spearman"),
            p_rho = cor.test(n_samples, sd_shannon, method = "spearman", exact = FALSE)$p.value,
            sd_median_cd = median(sd_shannon[diagnosis == "cd"]),
            sd_median_control = median(sd_shannon[diagnosis == "control"]),
            p_sd = wilcox.test(sd_shannon ~ diagnosis)$p.value,
            p_sd_c4_corrected = wilcox.test(sd_corrected ~ diagnosis)$p.value,
            # linear model: SD ~ diagnosis + number of samples
            p_sd_adj_n = coef_row(lm(sd_shannon ~ diagnosis + log(n_samples)), "diagnosiscd")$p)
glimpse(sd_tab)

# correction by subsampling: exactly k samples per person (random), 1000 repetitions
k_sub <- 4
sub_tab <- map_dfr(c("HMP2", "Hall"), function(co) {
  d <- long %>% filter(cohort == co) %>% group_by(subject_id) %>% filter(n() >= k_sub) %>% ungroup()
  set.seed(seed)   # per cohort, so that Hall does not depend on HMP2
  map_dfr(1:1000, function(i) {
    s <- d %>% group_by(subject_id, diagnosis) %>% slice_sample(n = k_sub) %>%
      summarise(sd_shannon = sd(Shannon), .groups = "drop")
    tibble(diff = median(s$sd_shannon[s$diagnosis == "cd"]) - median(s$sd_shannon[s$diagnosis == "control"]),
           p = wilcox.test(sd_shannon ~ diagnosis, data = s)$p.value)
  }) %>%
    summarise(cohort = co, k = k_sub, persons = n_distinct(d$subject_id),
              diff_median = median(diff), diff_q025 = quantile(diff, .025), diff_q975 = quantile(diff, .975),
              share_p_lt_005 = mean(p < .05))
})
sub_tab
sd_tab <- sd_tab %>% left_join(sub_tab, by = "cohort")
write_csv(sd_tab, file.path(outdir, "sd_vs_samples_per_person.csv"))

p_sd <- ggplot(sd_person, aes(n_samples, sd_shannon, colour = diagnosis)) +
  geom_point(size = 2.2, alpha = .75) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, colour = "black", linewidth = .7) +
  scale_colour_manual(values = col_diag) +
  facet_wrap(~ cohort, scales = "free_x") +
  labs(x = "samples per person", y = "SD of Shannon per person", colour = NULL,
       title = "Does the number of samples drive the SD?") + theme_plot
ggsave(file.path(outdir, "sd_vs_samples_per_person.png"), p_sd, width = 8.5, height = 4.5, dpi = 200)

# ---- 4. Time trend: why does Shannon rise only in controls? -------------
# candidates: (a) richness/evenness, (b) sequencing depth over time,
# (c) centre (13 of 14 HMP2 controls from M), (d) stratified by time
time_fit <- function(f, d) lmer(f, data = d, REML = FALSE,
                                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

trend_tab <- long %>%
  group_by(cohort) %>%
  group_modify(function(d, key) {
    out <- list(
      Shannon   = time_fit(Shannon ~ diagnosis * time + (1 | subject_id), d),
      Richness  = time_fit(Richness ~ diagnosis * time + (1 | subject_id), d),
      Evenness  = time_fit(evenness ~ diagnosis * time + (1 | subject_id), d),
      log_reads = time_fit(log_reads ~ diagnosis * time + (1 | subject_id), d),
      `Shannon + reads` = time_fit(Shannon ~ diagnosis * time + log_reads + (1 | subject_id), d)
    )
    if (key$cohort == "HMP2") {
      dm <- filter(d, site == "M")
      out$`Shannon, centre M only` <- time_fit(Shannon ~ diagnosis * time + (1 | subject_id), dm)
    }
    imap_dfr(out, function(m, nm) {
      # slope control = time; slope CD = time + interaction
      ctl <- coef_row(m, "time")
      L <- setNames(rep(0, length(fixef(m))), names(fixef(m))); L[c("time", "diagnosiscd:time")] <- 1
      cd <- contest(m, L, joint = FALSE)
      tibble(outcome = nm,
             slope_control = ctl$estimate, p_control = ctl$p,
             slope_cd = cd$Estimate, p_cd = cd$`Pr(>|t|)`,
             p_interaction = coef_row(m, "diagnosiscd:time")$p)
    })
  }) %>% ungroup()
print(trend_tab, width = Inf)
write_csv(trend_tab, file.path(outdir, "time_trend_explanations.csv"))

# stratified by time: CD effect in time windows (is the difference there at the start or does it develop?)
win_breaks <- c(0, 10, 20, 30, 40, Inf)
strat_tab <- long %>%
  mutate(window = cut(time, win_breaks, right = FALSE,
                      labels = c("0-9", "10-19", "20-29", "30-39", "40+"))) %>%
  group_by(cohort, window) %>%
  group_modify(function(d, key) {
    if (n_distinct(d$diagnosis) < 2) return(tibble())
    m <- if (max(table(d$subject_id)) > 1) lmer(Shannon ~ diagnosis + (1 | subject_id), data = d)
         else lm(Shannon ~ diagnosis, data = d)
    r <- coef_row(m, "diagnosiscd")
    tibble(samples = nrow(d), cd_persons = n_distinct(d$subject_id[d$diagnosis == "cd"]),
           control_persons = n_distinct(d$subject_id[d$diagnosis == "control"]),
           cd_effect = r$estimate, se = r$se, p = r$p)
  }) %>% ungroup()
print(strat_tab, n = Inf)
write_csv(strat_tab, file.path(outdir, "time_stratified.csv"))

p_strat <- ggplot(strat_tab, aes(window, cd_effect)) +
  geom_hline(yintercept = 0, colour = okabe_ito[["grey"]]) +
  geom_pointrange(aes(ymin = cd_effect - 1.96 * se, ymax = cd_effect + 1.96 * se), colour = col_diag["cd"]) +
  geom_text(aes(label = paste0(cd_persons, "/", control_persons)),
            y = -Inf, vjust = -.6, size = 3, colour = "black") +
  facet_wrap(~ cohort) +
  labs(x = "weeks since first sample", y = "Shannon CD - control (95 % CI)",
       title = "CD effect per time window", caption = "numbers = CD / control persons in the window") + theme_plot
ggsave(file.path(outdir, "time_stratified.png"), p_strat, width = 8.5, height = 4.5, dpi = 200)

# Shannon, richness, evenness, reads over time (HMP2): what changes?
p_trend <- long %>% filter(cohort == "HMP2") %>%
  pivot_longer(c(Shannon, Richness, evenness, log_reads), names_to = "metric") %>%
  mutate(metric = factor(metric, c("Shannon", "Richness", "evenness", "log_reads"),
                         c("Shannon", "Richness", "Evenness (Pielou)", "log10 reads"))) %>%
  ggplot(aes(time, value, colour = diagnosis)) +
  geom_point(size = .9, alpha = .3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_colour_manual(values = col_diag) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(x = "weeks since first sample", y = NULL, colour = NULL, title = "HMP2: what changes over time?") +
  theme_plot
ggsave(file.path(outdir, "hmp2_time_metrics.png"), p_trend, width = 11, height = 4.2, dpi = 200)

# ---- 5. Hall: is "Boston only" usable? --------------------------------
bo <- alpha$hall %>% filter(site == "Boston")
bo %>% group_by(diagnosis) %>%
  summarise(persons = n_distinct(subject_id), samples = n(),
            samples_per_person = paste(range(table(subject_id)), collapse = "-"))

m_bo  <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = bo)
m_bor <- lmer(Richness ~ diagnosis + (1 | subject_id), data = bo)
bo_person <- bo %>% group_by(subject_id, diagnosis) %>%
  summarise(Shannon = median(Shannon), Richness = median(Richness), .groups = "drop")
boston_tab <- bind_rows(
  coef_row(m_bo, "diagnosiscd") %>% mutate(analysis = "Shannon, mixed model"),
  coef_row(m_bor, "diagnosiscd") %>% mutate(analysis = "Richness, mixed model"),
  tibble(analysis = "Shannon, median per person (Wilcoxon)",
         estimate = median(bo_person$Shannon[bo_person$diagnosis == "cd"]) -
           median(bo_person$Shannon[bo_person$diagnosis == "control"]),
         p = wilcox.test(Shannon ~ diagnosis, data = bo_person, exact = TRUE)$p.value)
)
boston_tab

# power: which CD effect could be detected at all with 3 controls?
# simulation with the variance components from Boston (between / within person) and the
# actual number of samples per person; test = mixed model ~ t-test on person means
vc <- as.data.frame(VarCorr(m_bo))
sd_between <- vc$sdcor[vc$grp == "subject_id"]; sd_within <- vc$sdcor[vc$grp == "Residual"]
n_i <- bo %>% count(subject_id, diagnosis)
power_sim <- function(delta, n_sim = 2000) {
  set.seed(seed)   # same random numbers for every effect -> smooth power curve
  mean(replicate(n_sim, {
    mu <- rnorm(nrow(n_i), 0, sd_between) + rnorm(nrow(n_i), 0, sd_within / sqrt(n_i$n)) +
      delta * (n_i$diagnosis == "cd")
    t.test(mu[n_i$diagnosis == "cd"], mu[n_i$diagnosis == "control"], var.equal = TRUE)$p.value < .05
  }))
}
power_tab <- tibble(delta = seq(0, 0.8, by = 0.05)) %>% mutate(power = map_dbl(delta, power_sim))
power_tab
mde_80 <- power_tab$delta[which(power_tab$power >= .8)[1]]
cat("Boston: SD between persons", round(sd_between, 2), "| within", round(sd_within, 2),
    "| 80 % power only from |delta| =", mde_80, "\n")
# for comparison: CD effects in the other cohorts (Franzosa -0.37, HMP2 approx. -0.2)

p_bo <- ggplot(bo, aes(reorder(subject_id, Shannon, median), Shannon, colour = diagnosis)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 1.2, alpha = .6) +
  scale_colour_manual(values = col_diag) +
  labs(x = "person (Boston only)", y = "Shannon (genus)", colour = NULL,
       title = sprintf("Hall Boston: 9 CD vs 3 controls: 80 %% power only for |diff| >= %.2f", mde_80)) +
  theme_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "hall_boston.png"), p_bo, width = 8.5, height = 4.8, dpi = 200)

p_pow <- ggplot(power_tab, aes(delta, power)) +
  geom_hline(yintercept = .8, linetype = "dashed", colour = okabe_ito[["grey"]]) +
  geom_vline(xintercept = c(.2, .37), linetype = "dotted", colour = col_diag["cd"]) +
  annotate("text", x = c(.2, .37), y = .05, label = c("HMP2", "Franzosa"), hjust = -.1, size = 3.3,
           colour = col_diag["cd"]) +
  geom_line(colour = "black") + geom_point(colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "true Shannon difference CD - control", y = "power (α = 0.05)",
       title = "Hall Boston: power with 3 controls") + theme_plot
ggsave(file.path(outdir, "hall_boston_power.png"), p_pow, width = 6.5, height = 4.5, dpi = 200)

saveRDS(list(clinical = clinical, clinical_tab = clinical_tab, franzosa_outliers = fr_out, franzosa_effects = fr_effects,
             reads = reads_tab, sd = sd_tab, trend = trend_tab, strat = strat_tab,
             boston = boston_tab, boston_power = power_tab),
        file.path(outdir, "robustness_results.rds"))
