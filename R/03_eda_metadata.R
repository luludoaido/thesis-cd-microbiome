library(tidyverse)
library(readxl)
source("R/00_setup.R")

# =====================================================================
# EDA of the metadata (Hall, HMP2, Franzosa)
# 1-4: samples, diagnosis, age, NAs
# 5-9: longitudinal, calprotectin, HBI, medication, NA overview
# =====================================================================

path <- dir_processed

meta_franzosa <- read.csv(file.path(path, "franzosa/franzosa_meta.csv"))
meta_hall     <- read.csv(file.path(path, "hall/hall_meta.csv"))
meta_hmp2     <- read.csv(file.path(path, "hmp2/hmp2_meta.csv"), check.names = FALSE)

# remove technically faulty Hall sample (exclude_samples, see R/00_setup.R)
meta_hall <- meta_hall %>% filter(!ID %in% exclude_samples)

# remove samples with too low sequencing depth (min_reads, see R/00_setup.R); Franzosa has no read counts
low_hall <- low_read_samples("hall")
low_hmp2 <- low_read_samples("hmp2")
message("Hall < 1 million reads: ", paste(intersect(meta_hall$ID, low_hall), collapse = ", "))
message("HMP2 < 1 million reads: ", paste(intersect(meta_hmp2$`External ID`, low_hmp2), collapse = ", "))
meta_hall <- meta_hall %>% filter(!ID %in% low_hall)
meta_hmp2 <- meta_hmp2 %>% filter(!`External ID` %in% low_hmp2)

# Franzosa writes "Control", Hall/HMP2 "control" -> harmonise
meta_franzosa <- meta_franzosa %>% mutate(Diagnosis = recode(Diagnosis, Control = "control"))


# ---- 1. samples vs. people ----
# Franzosa: cross-sectional, 1 sample per person -> sample_id = person
n_samples_franzosa <- nrow(meta_franzosa)
n_people_franzosa  <- n_distinct(meta_franzosa$sample_id)

n_people_hall  <- n_distinct(meta_hall$subject_id)
n_samples_hall <- nrow(meta_hall)

sum_samples_hall <- meta_hall %>% count(subject_id) %>% summarise(
  min = min(n), median = median(n), max = max(n))   # samples per person

n_people_hmp2  <- n_distinct(meta_hmp2$`Participant ID`)
n_samples_hmp2 <- nrow(meta_hmp2)

sum_samples_hmp2 <- meta_hmp2 %>% count(`Participant ID`) %>% summarise(
  min = min(n), median = median(n), max = max(n))   # samples per person

# ---- 2. disease vs. control ----

tab_franzosa <- meta_franzosa %>%
  group_by(Diagnosis) %>%
  summarise(
    People  = n_distinct(sample_id),
    Samples = n(),
    .groups = "drop"
  ) %>%
  mutate(`Samples per person` = round(Samples / People, 1))

tab_hall <- meta_hall %>%
  group_by(diagnosis) %>%
  summarise(
    People  = n_distinct(subject_id),
    Samples = n(),
    .groups = "drop"
  ) %>%
  mutate(`Samples per person` = round(Samples / People, 1))

tab_hmp2 <- meta_hmp2 %>%
  group_by(diagnosis) %>%
  summarise(
    People  = n_distinct(`Participant ID`),
    Samples = n(),
    .groups = "drop"
  ) %>%
  mutate(`Samples per person` = round(Samples / People, 1))

tab_franzosa
tab_hall
tab_hmp2

# ---- 3. Age ----
# Hall only has age_category (no numeric age)
meta_franzosa %>% group_by(Diagnosis) %>%
  summarise(n = n(), median = median(Age, na.rm = TRUE), mean = mean(Age, na.rm = TRUE),
            q1 = quantile(Age, .25, na.rm = TRUE),
            q3 = quantile(Age, .75, na.rm = TRUE),
            min = min(Age, na.rm = TRUE), max = max(Age, na.rm = TRUE))

# HMP2: age per person (not per sample, otherwise persons with many samples count more)
meta_hmp2 %>% distinct(`Participant ID`, diagnosis, consent_age) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), median = median(consent_age, na.rm = TRUE), mean = mean(consent_age, na.rm = TRUE),
            q1 = quantile(consent_age, .25, na.rm = TRUE),
            q3 = quantile(consent_age, .75, na.rm = TRUE),
            min = min(consent_age, na.rm = TRUE), max = max(consent_age, na.rm = TRUE))

meta_hall %>% distinct(subject_id, diagnosis, age_category) %>% count(diagnosis, age_category)

# ---- 4. NAs ----
meta_franzosa %>% summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "spalte", values_to = "n_na") %>%
  arrange(desc(n_na))

meta_hmp2 %>% summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "spalte", values_to = "n_na") %>%
  arrange(desc(n_na))

meta_hall %>% summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "spalte", values_to = "n_na") %>%
  arrange(desc(n_na))


# =====================================================================
# All cohorts together (5-9)
# =====================================================================

# ---- Preparation: same column names per cohort -----------------------
# .person / .diag / .time, so that the code below need not be written three times
meta_list <- list(
  hall = meta_hall %>%
    mutate(.person = as.character(subject_id), .diag = diagnosis,
           .time = days_from_first_collection / 7),          # days -> weeks
  hmp2 = meta_hmp2 %>%
    mutate(.person = as.character(`Participant ID`), .diag = diagnosis,
           .time = week_num),
  franzosa = meta_franzosa %>%
    mutate(.person = as.character(sample_id), .diag = Diagnosis,
           .time = NA_real_)                                  # 1 sample per person
)

# time per person since the first sample (as in R/06_diversity.R)
# (depending on the person, Hall starts at 30-300 days, HMP2 at week 0-8)
meta_list[c("hall", "hmp2")] <- lapply(meta_list[c("hall", "hmp2")], function(m) {
  m %>% group_by(.person) %>% mutate(.time = .time - min(.time)) %>% ungroup()
})

# ---- 5. Longitudinal: what do the time points look like? -------------
# time columns: Hall = days_from_first_collection (+ visit_number), HMP2 = week_num (+ visit_num)
# Franzosa has none

# samples per person + follow-up duration per person
samples_pp <- imap_dfr(meta_list, function(m, nm) {
  m %>% group_by(.person, .diag) %>%
    summarise(n_samples = n(),
              followup_weeks = if (all(is.na(.time))) NA_real_ else max(.time, na.rm = TRUE) - min(.time, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(cohort = nm)
})

samples_pp %>%
  group_by(cohort, .diag) %>%
  summarise(people = n(),
            min = min(n_samples), median = median(n_samples), max = max(n_samples),
            median_followup_weeks = median(followup_weeks),
            max_followup_weeks = max(followup_weeks), .groups = "drop")

ggplot(filter(samples_pp, cohort != "franzosa"), aes(n_samples, fill = .diag)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_wrap(~ cohort, scales = "free") +
  labs(x = "samples per person", y = "persons", fill = "diagnosis")

# gaps between consecutive samples (weeks)
gaps <- imap_dfr(meta_list[c("hall", "hmp2")], function(m, nm) {
  m %>% arrange(.person, .time) %>% group_by(.person) %>%
    mutate(gap_weeks = .time - lag(.time)) %>% ungroup() %>%
    filter(!is.na(gap_weeks)) %>% transmute(cohort = nm, .diag, gap_weeks)
})
gaps %>% group_by(cohort) %>%
  summarise(median = median(gap_weeks), q1 = quantile(gap_weeks, .25),
            q3 = quantile(gap_weeks, .75), max = max(gap_weeks))

# time points per person (each row is one person)
imap_dfr(meta_list[c("hall", "hmp2")], function(m, nm) mutate(m, cohort = nm) %>%
           select(cohort, .person, .diag, .time)) %>%
  ggplot(aes(.time, fct_reorder(.person, .time, .fun = max), colour = .diag)) +
  geom_line(alpha = .4) + geom_point(size = .8) +
  facet_wrap(~ cohort, scales = "free") +
  labs(x = "weeks since first sample", y = "person", colour = "diagnosis") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ---- 6. Fecal Calprotectin --------------------------------------------
# NA = cohort has no value (HMP2 metadata contain no calprotectin)
calpro_col <- c(hall = "calprotectin", hmp2 = NA, franzosa = "Fecal.Calprotectin")

calpro <- imap_dfr(meta_list, function(m, nm) {
  col <- calpro_col[[nm]]
  if (is.na(col)) return(NULL)
  tibble(cohort = nm, person = m$.person, diagnosis = m$.diag,
         calpro = suppressWarnings(as.numeric(m[[col]])))
})

# how many samples / persons have a value at all?
calpro %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(), with_value = sum(!is.na(calpro)),
            persons_with_value = n_distinct(person[!is.na(calpro)]),
            median = median(calpro, na.rm = TRUE),
            q1 = quantile(calpro, .25, na.rm = TRUE),
            q3 = quantile(calpro, .75, na.rm = TRUE), .groups = "drop")

# CD vs. control (Franzosa: 1 sample per person -> Wilcoxon ok;
# Hall: several samples per person -> descriptive only, or median per person)
calpro %>% filter(cohort == "franzosa") %>%
  wilcox.test(calpro ~ diagnosis, data = .)

calpro %>% filter(cohort == "hall", !is.na(calpro)) %>%
  group_by(person, diagnosis) %>% summarise(calpro = median(calpro), .groups = "drop") %>%
  wilcox.test(calpro ~ diagnosis, data = .)

# distribution (log scale because strongly right-skewed; values = 0 are dropped here)
calpro %>% filter(!is.na(calpro), calpro > 0) %>%
  ggplot(aes(calpro, fill = diagnosis)) +
  geom_histogram(bins = 30, alpha = .7, position = "identity") +
  geom_vline(xintercept = c(50, 250), linetype = "dashed") +   # common cut-offs (µg/g)
  scale_x_log10() + facet_wrap(~ cohort, scales = "free") +
  labs(x = "faecal calprotectin (log10)")

# ---- 7. Harvey-Bradshaw-Index (HBI) -----------------------------------
# only Hall has HBI
hbi_col <- c(hall = "HBI", hmp2 = NA, franzosa = NA)

hbi <- imap_dfr(meta_list, function(m, nm) {
  col <- hbi_col[[nm]]
  if (is.na(col)) return(NULL)
  tibble(cohort = nm, person = m$.person, diagnosis = m$.diag,
         hbi = suppressWarnings(as.numeric(m[[col]])))
})

hbi %>%
  group_by(cohort, diagnosis) %>%
  summarise(samples = n(), with_value = sum(!is.na(hbi)),
            persons_with_value = n_distinct(person[!is.na(hbi)]),
            median = median(hbi, na.rm = TRUE), .groups = "drop")

# HBI only makes sense for CD -> controls with an HBI value would be a data error
hbi %>% filter(diagnosis == "control", !is.na(hbi))

# common categories: <5 remission, 5-7 mild, 8-16 moderate, >16 severe
hbi %>% filter(diagnosis == "CD", !is.na(hbi)) %>%
  mutate(aktivitaet = cut(hbi, c(-Inf, 4, 7, 16, Inf),
                          labels = c("remission", "mild", "moderate", "severe"))) %>%
  count(cohort, aktivitaet)

hbi %>% filter(diagnosis == "CD", !is.na(hbi)) %>%
  ggplot(aes(hbi)) + geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 4.5, linetype = "dashed") +           # remission cut-off
  labs(x = "HBI (CD only)", y = "samples")

# do HBI and calprotectin agree? (Hall, same samples)
meta_hall %>% filter(diagnosis == "CD", !is.na(HBI), !is.na(calprotectin)) %>%
  { print(cor.test(.$HBI, .$calprotectin, method = "spearman")); . } %>%
  ggplot(aes(HBI, calprotectin)) + geom_jitter(width = .15) + scale_y_log10()

# ---- 8. Medication ----------------------------------------------------
# only Franzosa has medication data (Hall and HMP2 metadata: none)
med_cols <- list(hall = character(0), hmp2 = character(0),
                 franzosa = c("antibiotic", "immunosuppressant", "mesalamine", "steroids"))

# per cohort and column: how many samples / persons have a value?
med_avail <- imap_dfr(meta_list, function(m, nm) {
  map_dfr(med_cols[[nm]], function(col) {
    has <- !is.na(m[[col]]) & m[[col]] != ""
    tibble(cohort = nm, variable = col,
           samples_total = nrow(m),
           samples_with_value = sum(has),
           persons_total = n_distinct(m$.person),
           persons_with_value = n_distinct(m$.person[has]))
  })
})
med_avail

# share of "Yes" per medication and diagnosis
imap_dfr(meta_list, function(m, nm) {
  if (length(med_cols[[nm]]) == 0) return(NULL)
  m %>% select(.diag, all_of(med_cols[[nm]])) %>%
    pivot_longer(-.diag, names_to = "medication", values_to = "value") %>%
    mutate(cohort = nm)
}) %>%
  group_by(cohort, medication, .diag) %>%
  summarise(n_yes = sum(value == "Yes", na.rm = TRUE), n_value = sum(!is.na(value)),
            share_yes = round(n_yes / n_value, 2), .groups = "drop")

# ---- 9. Overview: variable x cohort, how much is missing? -------------
# share of NA per column and cohort (empty cell = cohort does not have the column)
na_overview <- imap_dfr(meta_list, function(m, nm) {
  m %>% select(-starts_with(".")) %>%
    summarise(across(everything(), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "share_na") %>%
    mutate(cohort = nm)
}) %>%
  pivot_wider(names_from = cohort, values_from = share_na)
na_overview %>% arrange(variable)
