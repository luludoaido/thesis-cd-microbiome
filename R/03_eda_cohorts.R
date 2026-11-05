# --- EDA of the cohorts (Franzosa, HMP2, Hall) ---
# Who is in the data? CD vs. Control, samples per person, age, time points.

library(tidyverse)
source("R/00_setup.R")

outdir <- file.path(dir_results, "eda_cohorts")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Meta data ---
# the cleaned metadata from the 02_clean scripts; exclusions and time are already applied
read_meta <- function(cohort_id, name) {
  read_csv(file.path(dir_processed, cohort_id, paste0(cohort_id, "_meta.csv")), show_col_types = FALSE) %>%
    mutate(cohort = name, person = subject_id)
}

# Make one table, easier for analysis
meta <- bind_rows(read_meta("franzosa", "Franzosa"), read_meta("hall", "Hall"), read_meta("hmp2", "HMP2")) %>%
  select(cohort, sample_id, person, diagnosis, age, age_category, time) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("control", "cd")))

per_person <- meta %>%
  group_by(cohort, person, diagnosis) %>%
  summarise(n_samples = n(), followup_weeks = max(time) - min(time), .groups = "drop")

samples_person <- per_person %>%
  group_by(cohort, diagnosis) %>%
  summarise(persons = n(), samples = sum(n_samples), min_samples = min(n_samples), 
            median_samples = median(n_samples), max_samplex = max(n_samples), 
            median_followup_weeks = median(followup_weeks), max_followup_weeks = max(followup_weeks), 
            .groups = "drop")

print(samples_person, width = Inf)
write_csv(samples_person, file.path(outdir, "samples_persons.csv"))

# --- How many samples per person, in CD and controls ---
person_by_samples <- per_person %>%
  count(cohort, n_samples, diagnosis, name = "person") %>%
  pivot_wider(names_from = diagnosis, values_from = person, values_fill = 0) %>%
  select(cohort, samples_per_person = n_samples, control, cd)

print(person_by_samples, n = Inf)
write_csv(person_by_samples, file.path(outdir, "person_by_samples.csv"))

# panel titles with the number of people
panel_titles <- samples_person %>%
  select(cohort, diagnosis, persons) %>%
  pivot_wider(names_from = diagnosis, values_from = persons) %>%
  mutate(title = paste0(cohort, " (", control, " control, ", cd, " CD)"))
panel_titles <- set_names(panel_titles$title, panel_titles$cohort)

p_samples <- per_person%>%
  count(cohort, n_samples, diagnosis, name = "persons") %>%
  ggplot(aes(n_samples, persons, fill = diagnosis)) +
  geom_col(position = position_dodge(preserve = "single", width = .9), width = .9) +
  scale_x_continuous(breaks = function(x)
    unique(floor(pretty(x)))) +
  scale_fill_manual(values = col_diag) +
  facet_wrap(~ cohort, scales = "free", labeller = as_labeller(panel_titles)) +
  labs(x = "samples per person", y = "people", fill = NULL) +
  theme_plot

print(p_samples)
ggsave(file.path(outdir, "person_by_samples.png"), p_samples, 
       width = 10, height = 4, dpi = 200)

# --- Age ---
age_table <- meta %>%
  filter(!is.na(age)) %>%
  distinct(cohort, person, diagnosis, age) %>%
  group_by(cohort, diagnosis) %>%
  summarise(persons = n(), median = median(age), mean = mean(age),
            q1 = quantile(age, .25), q3 = quantile(age, .75),
            min = min(age), max = max(age), .groups = "drop")
print(age_table)
write_csv(age_table, file.path(outdir, "age.csv"))

meta %>%
  filter(cohort == "Hall") %>%
  distinct(person, diagnosis, age_category) %>%
  count(diagnosis, age_category) %>%
  print()

# --- Time points (Hall and HMP2) ---
longitudinal <- meta %>% filter(cohort != "Franzosa")

longitudinal %>%
  arrange(person, time) %>%
  group_by(person) %>%
  mutate(gap_weeks = time - lag(time)) %>%
  ungroup() %>%
  filter(!is.na(gap_weeks)) %>%
  group_by(cohort) %>%
  summarise(median = median(gap_weeks), q1 = quantile(gap_weeks, .25),
            q3 = quantile(gap_weeks, .75), max = max(gap_weeks)) %>%
  print()

# each row as one person
p_time <- ggplot(longitudinal, aes(time, fct_reorder(person, time, .fun = max), colour = diagnosis)) +
  geom_line(alpha = .4) +
  geom_point(size = .8) +
  scale_color_manual(values = col_diag) +
  facet_wrap(~ cohort, scales = "free") +
  labs(x = "weeks since first sample", y = "person", colour = NULL) +
  theme_plot +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
print(p_time)
ggsave(file.path(outdir, "time_points.png"), p_time, width = 9, height = 5, dpi = 200)
