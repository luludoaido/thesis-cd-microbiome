# --- EDA of the clinical data (Franzosa, Hall, HMP2) ---

# What this script does
# Describes calprotectin, HBI and medication in the three cohorts, and how much data each cohort
# has. It only describes the clinical data; the link to diversity is in 07_robustness.R.

# Decisions
# 1. Comparisons use one value per person (the median in Hall and HMP2), because samples of the
#    same person are not independent.
# 2. HMP2: the time window for matching clinical values is checked with 2, 1 and 0 weeks.
# 3. Calprotectin is shown on a log scale, because it is strongly right-skewed.
# 4. HBI is only analysed in CD.

library(tidyverse)
source("R/00_setup.R")

outdir <- file.path(dir_results, "eda_clinical")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

franzosa <- read_csv(file.path(dir_processed, "franzosa/franzosa_meta.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "Franzosa", sample_id, person = subject_id, diagnosis,
            calprotectin, antibiotic, immunosuppressant, mesalamine, steroids)

hall <- read_csv(file.path(dir_processed, "hall/hall_meta.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "Hall", sample_id, person = subject_id, diagnosis, calprotectin, HBI)

hmp2 <- read_csv(file.path(dir_processed, "hmp2/hmp2_meta.csv"), show_col_types = FALSE) %>%
  transmute(cohort = "HMP2", sample_id, person = subject_id, diagnosis, calprotectin = calpro, HBI = hbi,
            antibiotic_since_last_visit = abx, calpro_how, hbi_how, calpro_gap, hbi_gap)

meta <- bind_rows(franzosa, hall, hmp2) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))

# --- HMP2: how were the clinical values matched and does the time window matter? ---
hmp2_clinical <- meta %>% filter(cohort == "HMP2")

# number of samples per source
print(table(hbi = hmp2_clinical$hbi_how, useNA = "ifany"))
print(table(calprotectin = hmp2_clinical$calpro_how, useNA = "ifany"))

# Sensitivity analysis: repeat the comparison with shorter windows
# 2 = main analysis, 1 = at most 1 week apart, 0 = same visit only
gap_sensitivity <- map_dfr(c(2,1,0), function(g) {
  calpro_g <- hmp2_clinical %>% filter(!is.na(calprotectin),
                                       calpro_gap <= g)
  per_person <- calpro_g %>%
    group_by(person, diagnosis) %>%
    summarise(calprotectin = median(calprotectin), .groups = "drop")
  hbi_g <- hmp2_clinical %>% filter(diagnosis == "cd", !is.na(HBI), hbi_gap <= g)
  tibble(max_gap_weeks = g,
         calpro_samples = nrow(calpro_g),
         calpro_person_cd = sum(per_person$diagnosis == "cd"),
         calpro_person_control = sum(per_person$diagnosis == "control"),
         calpro_median_cd = median(calpro_g$calprotectin[calpro_g$diagnosis == "cd"]),
         calpro_median_control = median(calpro_g$calprotectin[calpro_g$diagnosis == "control"]),
         calpro_p_wilcoxon = wilcox.test(calprotectin ~ diagnosis, data = per_person)$p.value, 
         hbi_sample_cd = nrow(hbi_g),
         hbi_share_remission = round(mean(hbi_g$HBI < 5), 2))
})
print(gap_sensitivity, width = Inf)
write_csv(gap_sensitivity, file.path(outdir, "hmp2_gap_sensitivity.csv"))

# --- Calprotectin ---
# Note: in Hall only 3 of the 12 controls have a value, so the comparison is uncertain.

calpro <- meta %>% filter(!is.na(calprotectin))

calpro_table <- calpro %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(), person = n_distinct(person), median = median(calprotectin),
            q1 = quantile(calprotectin, .25), q3 = quantile(calprotectin, .75), .groups = "drop")

print(calpro_table)
write_csv(calpro_table, file.path(outdir, "calprotectin.csv"))

# CD vs. control with the Wilcoxon rank-sum test (Wilcoxon 1945, Mann & Whitney 1947), which needs
# no normal distribution. One value per person (median for Hall and HMP2).
calpro_test <- calpro %>%
  group_by(cohort, person, diagnosis) %>%
  summarise(calprotectin = median(calprotectin), .groups = "drop") %>%
  group_by(cohort) %>%
  summarise(p_wilcoxon = wilcox.test(calprotectin ~diagnosis)$p.value)

print(calpro_test)
write_csv(calpro_test, file.path(outdir, "calprotectin_test.csv"))

# Log scale because strongly right-skewed; values of 0 are dropped.
# Lines = commonly used cut-offs: there is no single accepted value. Below 50 is the most sensitive
# range (Rokkas 2018), 250 indicates large ulcers (D'Haens 2012). 
p_calpro <- ggplot(filter(calpro, calprotectin > 0),
                   aes(calprotectin, fill = diagnosis)) +
  geom_histogram(bins = 30, alpha = .7, position = "identity")+
  geom_vline(xintercept = c(50, 250), linetype = "dashed")+
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE))+
  scale_fill_manual(values = col_diag) +
  facet_wrap( ~cohort, scales = "free_y") +
  labs(x = "faecal calprotectin (log10)", y = "samples", fill = NULL) +
  theme_plot

print(p_calpro)  
ggsave(file.path(outdir, "calprotectin.png"), p_calpro, width = 10, height = 4, dpi = 200)

# --- HBI (Hall and HMP2, CD only) ---
# HBI is a symptom score for Crohn's disease (Harvey & Bradshaw 1980)
hbi <- meta %>% filter(diagnosis == "cd", !is.na(HBI))

# Common categories: <5 remission, 5-7 mild, 8-16 moderate, >16 severe.
# A score of 4 or less matches a CDAI of 150 or less, the usual definition of remission (Vermeire 2010).
hbi_table <- hbi %>%
  mutate(activity = cut(HBI, c( -Inf, 4, 7, 16, Inf),
                        labels = c("remission", "mild", "moderate", "severe"))) %>% 
  count(cohort, activity, name = "samples") %>%
  group_by(cohort) %>%
  mutate(share = round(samples/ sum(samples), 2)) %>%
  ungroup()

print(hbi_table)
write_csv(hbi_table, file.path(outdir, "hbi.csv"))  

p_hbi <- ggplot(hbi, aes(HBI)) +
  geom_histogram(binwidth = 1, fill = col_diag[["cd"]]) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  facet_wrap(~ cohort, scales = "free_y") +
  labs(x = "HBI (CD only)", y = "samples")+
  theme_plot

print(p_hbi)  
ggsave(file.path(outdir, "hbi.png"), p_hbi, width = 8, height = 4, dpi = 200)  
  
# Do HBI and calprotectin agree? (CD samples with both values)
# They need not: HBI measures symptoms, and calprotectin was not related to symptom scores in CD (D'Haens 2012).
hbi %>%
  filter(!is.na(calprotectin)) %>% 
  group_by(cohort) %>%
  summarise(samples = n(),
            rho = cor(HBI, calprotectin, method = "spearman"),
            p = cor.test(HBI, calprotectin, method = "spearman", exact = FALSE)$p.value) %>%
  print()

# --- Medication ---
# Franzosa: medication at the time of sampling
medication <- meta %>%
  filter(cohort == "Franzosa") %>%
  select(diagnosis, antibiotic, immunosuppressant, mesalamine, steroids) %>%
  pivot_longer(-diagnosis, names_to = "medication", values_to = "value")%>%
  group_by(medication, diagnosis) %>% 
  summarise(n_yes = sum(value == "Yes", na.rm = TRUE), n_value = sum(!is.na(value)),
            share_yes = round(n_yes /n_value, 2), .groups = "drop")

print(medication)
write_csv(medication, file.path(outdir, "medication_franzosa.csv"))
  
# HMP2: antibiotic since last visit
antibiotics_hmp2 <- meta %>%
  filter(cohort == "HMP2") %>%
  group_by(diagnosis) %>%
  summarise(samples = n(),
            n_yes = sum(antibiotic_since_last_visit == "Yes", na.rm = TRUE),
            n_values = sum(!is.na(antibiotic_since_last_visit)),
            person_with_yes = n_distinct(person[antibiotic_since_last_visit %in% "Yes"]))

print(antibiotics_hmp2)  
write_csv(antibiotics_hmp2, file.path(outdir, "antibiotics_hmp2.csv"))


# --- Missing Values ---
# 1 = cohort does not have this variable, output = the share of NA per cohort
missing_values <- meta %>%
  group_by(cohort) %>%
  summarise(across(c(calprotectin, HBI, antibiotic, antibiotic_since_last_visit, immunosuppressant, mesalamine, steroids),
                   ~ round(mean(is.na(.)), 2)))

print(missing_values, width = Inf)
write_csv(missing_values, file.path(outdir, "missing_values.csv"))

# --- Literature ---
# D'Haens G et al. (2012). Fecal calprotectin is a surrogate marker for endoscopic lesions in
#   inflammatory bowel disease. Inflamm Bowel Dis 18(12):2218-2224. doi:10.1002/ibd.22917
# Harvey RF, Bradshaw JM (1980). A simple index of Crohn's-disease activity. Lancet 315(8167):514.
#   doi:10.1016/s0140-6736(80)92767-1
# Mann HB, Whitney DR (1947). On a test of whether one of two random variables is stochastically
#   larger than the other. Ann Math Stat 18(1):50-60. doi:10.1214/aoms/1177730491
# Rokkas T et al. (2018). Fecal calprotectin in assessing inflammatory bowel disease endoscopic activity:
#   a diagnostic accuracy meta-analysis. J Gastrointestin Liver Dis 27(3):299-306. doi:10.15403/jgld.2014.1121.273.pti
# Vermeire S et al. (2010). Correlation between the Crohn's disease activity and Harvey-Bradshaw indices
#   in assessing Crohn's disease severity. Clin Gastroenterol Hepatol 8(4):357-363. doi:10.1016/j.cgh.2010.01.001
# Wilcoxon F (1945). Individual comparisons by ranking methods. Biometrics Bull 1(6):80-83. doi:10.2307/3001968
