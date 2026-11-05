# --- Alpha diversity ---
# Does the diversity within a sample differ between CD and controls? 
# looks at both levels (species and genus)
# Started by run_levels.R, which sets cohort and level.

library(tidyverse)
library(vegan)
library(lme4)      
library(lmerTest)  
source("R/00_setup.R")

source("R/04_load_abundance.R")

# --- Samples with very few genera ---
p_few <- ggplot(few, aes(n_genera, fill = diagnosis)) +
  geom_histogram(bins = 30, alpha = .7, position = "identity") +
  facet_wrap(~ threshold, scales = "free") +
  theme_plot

save_plot(p_few, "few_genera.png", 9, 4)

# --- Alpha diversity per sample ---
alpha <- imap_dfr(X_thr, function(X, nm) {
  tibble(sample_id = rownames(X),
         Richness = rowSums(X > 0),
         Shannon  = diversity(X, "shannon"),
         Simpson  = diversity(X, "simpson"),
         threshold = nm)
}) %>% left_join(meta, by = "sample_id")

# --- CD vs. control: first sample of each person ---
baseline_tab <- alpha %>% 
  filter(baseline) %>%
  group_by(threshold) %>%
  summarise(p_shannon_baseline = wilcox.test(Shannon ~diagnosis)$p.value,
            p_richness_baseline = wilcox.test(Richness ~diagnosis)$p.value)
print(baseline_tab)

alpha_tests <- baseline_tab

p_alpha <- alpha %>% 
  filter(baseline) %>%
  ggplot(aes(diagnosis, Shannon, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .15, size = .8) +
  facet_wrap(~ threshold) +
  labs(title = "Shannon, baseline samples") +
  theme_plot

save_plot(p_alpha, "alpha_baseline.png", 9, 4)

# --- CD vs. control: all samples (Hall and HMP2) ---
if (longitudinal){
  #median per person
  alpha_person <- alpha %>% 
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_samples = n(), Shannon = median(Shannon), Richness = median(Richness),
              .groups = "drop")
  
  median_tab <- alpha_person %>%
    group_by(threshold) %>%
    summarise(p_shannon = wilcox.test(Shannon ~ diagnosis)$p.value,
              p_richness = wilcox.test(Richness ~ diagnosis)$p.value)
  print(median_tab)
  
  p_alpha_med <- ggplot(alpha_person, aes(diagnosis, Shannon, fill = diagnosis))+
    geom_boxplot(outlier.shape = NA, alpha = .6) + geom_jitter(width = .15, size = .8) +
    facet_wrap(~ threshold) +
    labs(title = "Shannon, median per person") +
    theme_plot
  
  save_plot(p_alpha_med, "alpha_median_person.png", 9, 4)
  
  #mixed model: all samples with person as random effect
  mixed_tab <- alpha %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      m <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d)
      vc <- as.data.frame(VarCorr(m))
      tibble(estimate_cd = fixef(m)[2], p_mixed = coef(summary(m))[2, "Pr(>|t|)"],
             icc = vc$vcov[1]/sum(vc$vcov)) #calculates the share of the variance between people
    })
  print(mixed_tab)
  
  alpha_tests <- alpha_tests %>%
    left_join(median_tab, by = "threshold") %>%
    left_join(mixed_tab, by = "threshold")
  
  p_spag <- alpha %>% 
    filter(threshold == main_thr) %>%
    arrange(subject_id, time) %>%
    ggplot(aes(time, Shannon, group = subject_id, colour = diagnosis)) +
    geom_line(alpha = .4) +
    geom_point(size = .8, alpha = .6) +
    labs(x = "weeks since first sample") +
    theme_plot
  save_plot(p_spag, "alpha_spaghetti.png", 8, 5)
}

# --- Stability within a person (Hall and HMP2) ---
if (longitudinal){
  #variation of Shannon, for people with at least 3 samples
  alpha_var <- alpha %>%
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_samples = n(), sd_shannon = sd(Shannon), range_shannon = max(Shannon) - min(Shannon),
              .groups = "drop") %>%
    filter(n_samples >= 3)
  
  alpha_var_test <- alpha_var %>%
    group_by(threshold) %>%
    summarise(p_sd_shannon = wilcox.test(sd_shannon ~diagnosis)$p.value)
  
  print(alpha_var_test)
  # Bray-Curtis distance between all pairs of samples (0 = identical, 1 = nothing shared)
  beta_pairs <- imap_dfr(X_thr, function(X,nm){
    D <- as.matrix(vegdist(X, method = "bray"))
    m <- meta[match(rownames(X), meta$sample_id), ]
    idx <- which(upper.tri(D), arr.ind = TRUE)
    tibble(threshold = nm, s1 = rownames(X)[idx[, 1]], 
           s2 = rownames(X)[idx[, 2]],
           p1 = m$subject_id[idx[, 1]],
           p2 = m$subject_id[idx[, 2]],
           d1 = m$diagnosis[idx[, 1]],
           d2 = m$diagnosis[idx[, 2]],
           bray = D[idx])
  })
  # same person vs. other person with the same diagnosis
  within_between <- beta_pairs %>%
    filter(d1 == d2) %>%
    mutate(diagnosis = d1, 
           comparison = if_else(p1 == p2, "same person", "other person"))
  
  p_within <- within_between %>%
    filter(threshold == main_thr) %>%
    ggplot(aes(comparison, bray, fill = diagnosis)) +
    geom_boxplot(outlier.size = .3) +
    labs(x = NULL, y = "Bray-Curtis distance", title = "Samples of the same vs. different persons") +
    theme_plot
  
  save_plot(p_within, "beta_within_between.png", 7, 5)
  
  #distance between two samples in a row of the same person
  next_sample <- meta %>%
    filter(!is.na(time)) %>%
    arrange(subject_id, time) %>%
    group_by(subject_id) %>%
    mutate(next_id = lead(sample_id), gap_weeks = lead(time) - time) %>%
    ungroup() %>%
    filter(!is.na(next_id)) %>%
    select(subject_id, diagnosis, s_a = sample_id, s_b = next_id, gap_weeks)
  
  consecutive <- beta_pairs %>%
    select(threshold, s1, s2, bray) %>%
    bind_rows(rename(., s1 = s2, s2 = s1)) %>%
    inner_join(next_sample, by = c("s1" = "s_a", "s2" = "s_b"))
  
  stability <- consecutive %>%
    group_by(threshold, subject_id, diagnosis) %>%
    summarise(n_pairs = n(), median_bray_consec = median(bray), .groups = "drop")
  
  stability_test <- stability %>%
    group_by(threshold) %>%
    summarise(p_stability = wilcox.test(median_bray_consec ~diagnosis)$p.value)
  print(stability_test)
  write_csv(left_join(alpha_var_test, stability_test, by = "threshold"), 
            file.path(outdir, "stability_tests.csv"))
  
  # lm
  p_gap <- consecutive %>% filter(threshold == main_thr) %>%
    ggplot(aes(gap_weeks, bray, colour = diagnosis)) +
    geom_jitter(width = .2, alpha = .5) + 
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = "weeks between two consecutive samples", y = "Bray-Curtis distance") + 
    theme_plot
  save_plot(p_gap, "beta_consecutive_gap.png", 7, 5)
}
  
# --- Change over time (Hall and HMP2) ---
if (longitudinal) {
  # time = weeks since the first sample; the interaction tests if 
  # CD and controls change differently
  # REML = FALSE so that the two models can be compared
  time_tab <- alpha %>%
    group_by(threshold) %>%
    group_modify(function(d, key) {
      m_null <- lmer(Shannon ~ diagnosis + (1 | subject_id), data = d, REML = FALSE)
      m_time <- lmer(Shannon ~ diagnosis * time + (1 | subject_id), data = d, REML = FALSE)
      cf <- coef(summary(m_time))
      tibble(cd_at_start   = cf["diagnosiscd", "Estimate"],
             p_cd_at_start = cf["diagnosiscd", "Pr(>|t|)"],
             slope_control = cf["time", "Estimate"],
             p_slope_control = cf["time", "Pr(>|t|)"],
             slope_cd = cf["time", "Estimate"] + cf["diagnosiscd:time", "Estimate"],
             p_interaction = cf["diagnosiscd:time", "Pr(>|t|)"],
             p_lrt = anova(m_null, m_time)$`Pr(>Chisq)`[2])
    })
  print(time_tab)
  write_csv(time_tab, file.path(outdir, "time_model.csv"))
  
  # same model with a separate slope per person
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
  write_csv(time_slope_tab, file.path(outdir, "time_model_slopes.csv"))
  
  # how many people with at least 5 samples increase over time?
  person_slopes <- alpha %>% filter(threshold == main_thr) %>%
    group_by(subject_id, diagnosis) %>% filter(n() >= 5) %>%
    summarise(n = n(), slope = coef(lm(Shannon ~ time))[2], .groups = "drop")
  print(person_slopes %>% group_by(diagnosis) %>%
          summarise(persons = n(), increasing = sum(slope > 0), median_slope = median(slope)))
  
  p_trend <- alpha %>% filter(threshold == main_thr) %>%
    ggplot(aes(time, Shannon, colour = diagnosis)) +
    geom_point(alpha = .25, size = .8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(x = "weeks since first sample", title = "Shannon over time (linear trends per diagnosis)") + theme_plot
  save_plot(p_trend, "alpha_time_trend.png", 8, 5)
}

# --- Age: does the cd effect change when adjusting for age? ---
# Hall has no age
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
  write_csv(age_tab, file.path(outdir, "age_model.csv"))
}

# --- Antibiotics: same test without people on antibiotics ---
# only Franzosa
if (any(alpha$antibiotic == "Yes", na.rm = TRUE)) {
  abx_tab <- alpha %>%
    filter(baseline) %>%
    mutate(no_abx = antibiotic != "Yes" | is.na(antibiotic)) %>%
    group_by(threshold) %>%
    summarise(p_shannon_all = wilcox.test(Shannon ~ diagnosis)$p.value,
              p_richness_all = wilcox.test(Richness ~ diagnosis)$p.value,
              n_without_abx = sum(no_abx),
              p_shannon_without_abx = wilcox.test(Shannon[no_abx] ~ diagnosis[no_abx])$p.value,
              p_richness_without_abx = wilcox.test(Richness[no_abx] ~ diagnosis[no_abx])$p.value)
  print(abx_tab)
  write_csv(abx_tab, file.path(outdir, "without_antibiotics.csv"))
}

# --- Centre: does the cd effect depend on where the samples come from? ---
site_tab <- alpha %>%
  group_by(threshold) %>%
  group_modify(function(d, key) {
    # centres with at least 3 persons of both diagnoses
    ok <- d %>% distinct(subject_id, diagnosis, site) %>% count(site, diagnosis) %>%
      group_by(site) %>% filter(n_distinct(diagnosis) == 2, min(n) >= 3) %>% pull(site) %>% unique()
    if (length(unique(d$site)) < 2 || length(ok) == 0) return(tibble(note = "centre cannot be separated from diagnosis"))
    main_site <- d %>% filter(site %in% ok, diagnosis == "control") %>% count(site) %>% slice_max(n, n = 1) %>% pull(site)
    
    # CD effect and its p-value
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

# --- Save ---
write_csv(alpha_tests, file.path(outdir, "alpha_tests.csv"))
write_csv(site_tab, file.path(outdir, "alpha_by_centre.csv"))
saveRDS(list(alpha = alpha, few = few), file.path(outdir, "alpha_results.rds")) 