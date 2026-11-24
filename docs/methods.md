# Methods and decisions

How the data are selected and analysed, why, and where it is in the code.
Every entry has a number. The scripts refer to it as `(see methods 4.3)`.

- New entries go at the end of the matching section with the next free number. Numbers are never reused or renumbered.
- Changes go under **Changes** with the date, so the history stays visible.
- Rules that apply to several scripts are variables or functions in `R/00_setup.R`.
- Short source names (e.g. Gevers 2014) are listed in full under [Literature](#literature).

## Overview: from raw data to analysis

As of 2026-10-28, genus level, main analysis (10 % prevalence filter).

| Cohort | Raw data | After selection (metadata file) | In the analysis | Persons (CD / control) |
|---|---|---|---|---|
| Franzosa 2019 | 220 samples | 144 (without UC) | 144 | 88 / 56 (1 sample per person) |
| Hall 2017 | 259 samples, 32 persons | 162 (without UC) | 158 | 9 / 12 |
| HMP2 | 1638 metagenome samples, 130 persons | 578 (CD + controls, adults) | 571 | 30 / 14 |

| Section | Topic | Scripts |
|---|---|---|
| 0 | Project decisions | – |
| 1 | Rules for all cohorts | `00_setup.R`, all |
| 2 | Franzosa 2019 | `02_clean_franzosa.R` |
| 3 | Hall 2017 | `01_download_cmd.R`, `02_clean_hall.R` |
| 4 | HMP2 | `02_clean_hmp2.R` |
| 5 | Clinical data | `03_eda_clinical.R`, `07_robustness.R` |
| 6 | Cohort description | `03_eda_cohorts.R` |
| 7 | Diversity | `04_load_abundance.R`, `05_alpha_diversity.R`, `05_beta_diversity.R` |
| 8 | Most and least abundant taxa | `06_extremes.R` |
| 9 | Robustness | `07_robustness.R` |
| 10 | H1: dysbiosis | `08_dysbiosis.R` |
| 11 | H2: fibre-associated taxa | `09_fibre_taxa.R` |
| 12 | Diet- and CD-associated taxa | `10_species_cd_diet.R` |
| 13 | Disease severity within CD | `11_severity_taxa.R` |

---

## 0 Project decisions

### 0.1 Case-control comparison
- **Decision:** CD is compared with healthy controls; the thesis does not only describe CD.
- **Why:** fits the assignment; SQ1 (dysbiosis) is a group comparison by definition.
- **Date:** 2026-08-20 (meeting)

### 0.2 Timeline
- **Decision:** submission on 2027-01-07, four days before the official deadline (2027-01-11).
- **Date:** 2026-08-26

### 0.3 Search criteria for datasets
- **Decision:** paired same-stool samples are not needed, because groups are compared, not samples. Any metabolomics platform is fine if SCFAs are well annotated. Diet metadata are not used.
- **Why:** untargeted HMP2 data annotate butyrate and propionate with high confidence; diet is approached through microbiome proxies.
- **Date:** 2026-09-18

### 0.4 Excluded dataset candidates
- **Decision:** about 50 candidates → 13 with both omics → 8 checked in detail → 2 meet every criterion (HMP2, Franzosa). Reasons per candidate in `Datensatz_Tracker_CD_Thesis_v3.xlsx`.
- **Date:** 2026-09-21

### 0.5 No regional balance
- **Decision:** the goal of 2 North American, 2 European and 2 Asian cohorts is dropped.
- **Why:** only cohorts from the USA and the Netherlands remained.
- **Date:** 2026-09-22

### 0.6 Measured butyrate from HMP2
- **Decision:** HMP2 is the cohort for measured butyrate and propionate (relative peak intensities, not mM), compared within the cohort only.
- **Why:** the only cohort with reliably annotated butyrate (HMDB0000039) and propionate (HMDB0000237); acetate is missing.
- **Status:** part of H3, work in progress (0.7). The HMP2 metabolomics are not yet in `data/processed`.
- **Date:** 2026-09-18

### 0.7 SQ3 / H3
- **Plan so far:** H3b (enzymes, EC level) as main analysis; H3a (measured butyrate and propionate) as exploratory check whether enzymes and metabolites agree.
- **Why:** measured butyrate rests on 14 adult controls. Franzosa has enzymes but no pathways, so the enzyme level keeps H3b comparable across both cohorts.
- **Status:** work in progress; how exactly H3 is done is not decided yet.
- **Date:** 2026-09-22; status updated 2026-11-16

### 0.8 No own processing of raw reads
- **Decision:** IBDomics and Sheba are not processed from raw reads (MetaPhlAn + HUMAnN) for now.
- **Why:** it would only extend SQ1/SQ2/H3b, needs a cluster and about 1.5–2 TB.
- **Status:** deferred; only if time is left.
- **Date:** 2026-09-29

### 0.9 Few controls for measured butyrate
- **Decision:** butyrate stays exploratory and is reported with effect size and CI. 1000IBD (EGA) is not requested.
- **Why:** EGA access takes weeks to months. 1000IBD (238 / 255) could detect d = 0.25.
- **Status:** deferred.
- **Date:** 2026-09-29

### 0.10 Cohorts are analysed separately
- **Decision:** each cohort is analysed on its own (HMP2 mixed model, Franzosa standard model with sub-cohort) and the results are compared side by side. No meta-analysis and no pooling of raw data. Measured butyrate (HMP2 only) is reported alone.
- **Why:** different processing (2.1), different patient groups, and HMP2 would dominate (≈ 571 vs. 144 samples).
- **Date:** 2026-09-22
- **Changes:** 2026-11-16: the planned fixed-effect meta-analysis and the pooled variant are dropped.

### 0.11 Which data for which question
- **Decision:** SQ1, SQ2 and H3b need only shotgun metagenomes; H3a needs measured metabolites. SQ1 and SQ2 use Franzosa and HMP2 (1.10); SQ3 uses HMP2 and Franzosa (work in progress, 0.7).
- **Date:** 2026-09-22
- **Changes:** 2026-11-09: Hall no longer used for the hypotheses.

### 0.12 Other curatedMetagenomicData cohorts
- **Decision:** Nielsen 2014, LifeLinesDeep 2016, VilaAV 2018 and one more cohort are not used.
- **Why:** they only distinguish "IBD", so CD cannot be separated; VilaAV also has no healthy controls.
- **Date:** 2026-09-29

### 0.13 Data handling
- **Decision:** data are downloaded and kept locally; raw and cleaned tables as csv; per cohort separate tables for metadata, taxa and function; code in a GitHub repository.
- **Why:** formats differ between projects (long vs. wide, several IDs per person).
- **Date:** 2026-09-29

### 0.14 Order of the analyses and machine learning
- **Decision:** variability and diversity first, then taxa, then metabolites. PCA and UMAP are compared. Machine learning only if more power is needed at the end.
- **Why:** diversity is the basis for the other analyses.
- **Date:** 2026-09-29

---

## 1 Rules for all cohorts

### 1.1 CD and controls only, no UC
- **Rule:** samples with ulcerative colitis (UC) are removed.
- **Why:** the thesis compares Crohn's disease with healthy controls.
- **Code:** the three `02_clean_*.R` scripts.
- **Effect:** Franzosa −76 samples, Hall −97 samples, HMP2 −459 samples.
- **Changes:** done by hand in Excel until 2026-10-28, in code since then.

### 1.2 Minimum sequencing depth: 1 million reads
- **Rule:** samples with fewer than `min_reads` = 1,000,000 reads are excluded.
- **Why:** with very few reads only the 2–3 most abundant genera are found, which distorts richness and Shannon. Example: Hall p9216_mo2 has 48,368 reads and 2 genera, the other samples of the same person have 41–49 genera. The limit is an own decision; shallow shotgun profiles were tested down to 0.5 million reads.
- **Source:** Hillmann 2018.
- **Read count:** Hall `number_reads` (curatedMetagenomicData), HMP2 `reads_filtered` (after QC and removal of human reads). Franzosa has no read counts, so the rule does not apply there.
- **Code:** `min_reads`, `finish_meta()` in `R/00_setup.R`.
- **Effect:** Hall −3 samples (p9216_mo2, p9281_mo6, SKST037_3_G103026), HMP2 −6 samples. No person is lost completely.

### 1.3 Prevalence filter at genus level
- **Rule:** a genus is kept if it occurs in at least 10 % of the samples of a cohort. At species level, a species stays if its genus stays, with no further filter on the species. Species without a genus cannot be filtered and stay. 0 % and 5 % run as sensitivity analyses.
- **Why:** very rare taxa are mostly noise and make the CLR transformation unstable. Filtering at genus level keeps species that are rare themselves but belong to a common genus. Results can depend on such steps, hence the sensitivity analyses. The exact cut-off is an own decision.
- **Source:** Cao 2021; Nearing 2022 (shown for 16S data).
- **Code:** `passes_filter()` in `R/00_setup.R`, used in `04`, `06`, `08`, `09`.
- **Effect:** species level: Franzosa 191 of 201 species, HMP2 286 of 436.
- **Changes:** 2026-11-05: species are filtered through their genus. Before, species were filtered on their own.

### 1.4 Relative abundances
- **Rule:** every sample is divided by its own total (total sum scaling). After filtering, it is scaled back to a total of 1.
- **Why:** samples differ in sequencing depth.
- **Code:** `R/04_load_abundance.R`.

### 1.5 Repeated samples per person
- **Rule:** in Hall and HMP2 the main test uses all samples: a mixed model with person as random effect, or for PERMANOVA a permutation of whole persons (7.9). One value per person (earliest sample `baseline`, or median) is used as check, and for measures that are defined per person anyway (SD, stability, 7.4, 9.5). Descriptive comparisons use one value per person (median). Franzosa has one sample per person.
- **Why:** several samples of the same person are not independent. In HMP2 the ICC was 0.58, so 141 samples correspond to about 53 independent observations.
- **Covariates:** sex is not used (decided 2026-11-16).
- **Source:** Bates 2015 (lme4); Kuznetsova 2017 (lmerTest, p-values by Satterthwaite).
- **Code:** `finish_meta()` in `R/00_setup.R` sets `baseline` and `time` (weeks since the first sample).
- **Changes:** 2026-11-16: mixed model / person permutation became the main test in `03`, `05`, `08` and `09`; before, several tests used only the earliest sample or one value per person, and the HBI–calprotectin correlation in `03` treated all samples as independent.

### 1.6 Centre and batch correction
- **Rule:** in tests, the centre (`site`) is a covariate. Corrected data are used for figures only (PCA/UMAP), and only if every centre has both diagnoses. ComBat is used for the main figures; limma and MMUPHin are compared.
- **Why:** correcting with diagnosis in the model and then testing makes p-values too optimistic. In Hall, Atlanta has only controls, so a correction would separate artificially; Hall stays uncorrected.
- **Source:** Nygaard 2016; Johnson 2007 (ComBat); Ritchie 2015 (limma); Ma 2022 (MMUPHin).
- **Centre:** Franzosa 2.3, Hall 3.3, HMP2 4.7.
- **Code:** batch section in `R/05_beta_diversity.R`.

### 1.7 Fixed seed
- **Rule:** `seed <- 42`, set again directly before permutation tests, UMAP and simulations. The value has no meaning.
- **Why:** same results on every run. Two runs of all scripts gave identical files (checked 2026-10-28).
- **Code:** `R/00_setup.R`, `05_beta_diversity.R`, `07_robustness.R`.

### 1.8 Two-group test without normal distribution
- **Rule:** comparisons of two groups with one value per person use the Wilcoxon rank-sum test.
- **Why:** the values (diversity, abundances, clinical scores) are often not normally distributed.
- **Source:** Wilcoxon 1945; Mann & Whitney 1947.

### 1.9 Colours
- **Rule:** Okabe-Ito palette.
- **Why:** readable with colour blindness.
- **Source:** Wong 2011.
- **Code:** `okabe_ito`, `col_diag` in `R/00_setup.R`.

### 1.10 Hypothesis tests only with Franzosa and HMP2
- **Rule:** the hypothesis scripts (`08`, `09`) use Franzosa and HMP2 only. Hall stays in the descriptive and diversity analyses.
- **Why:** part of the Hall cohort was processed with a different lab kit, so the samples within the cohort are no longer comparable and the data cannot be used for the hypotheses.
- **Changes:** Hall was used for SQ1/SQ2 from 2026-09-29 (9 CD / 12 control persons after exclusions). Not used for the hypotheses since 2026-11-09; reason recorded 2026-11-16.

### 1.11 Main test and checks
- **Rule:** continuous outcomes in Franzosa (one sample per person): Wilcoxon rank-sum test is the main test; lm with centre is the check and gives the effect size. Hall and HMP2: mixed model on all samples (1.5), Wilcoxon on the median per person as check. Every model effect is reported with a 95 % CI (t quantile, Satterthwaite df for lmer; `ci_cd()` in `R/00_setup.R`). Exceptions: binary outcomes (Fisher / logistic models, 10.6, 11.2) and H2b (MaAsLin2, 11.3).
- **Why:** the values are often skewed; Wilcoxon needs no normal distribution. The centre was checked in Franzosa: no effect on Shannon (p = 1.0) and MD-index (p = 0.31, 13 taxa; 0.47 with the 7 taxa of the summary), but on richness (p < 0.001), Bray-Curtis score (p = 0.03) and Aitchison score (p < 0.001). The CD effect stays the same after adjustment, so a Wilcoxon stratified by centre (van Elteren) is not needed. Residuals of lm/lmer are not normal but nearly symmetric; both are robust to this.
- **Source:** Wilcoxon 1945; Mann & Whitney 1947; Lumley 2002; Schielzeth 2020; Satterthwaite 1946.
- **Changes:** 2026-11-16: HMP2 Wilcoxon check in `08` changed from the earliest sample to the median per person; 95 % CI added in `08` and `09`.

---

## 2 Franzosa 2019

### 2.1 Data source and processing
- **Data:** Supplementary Datasets 1, 2, 4 and 6 of Franzosa 2019.
- **Processing:** species and enzymes come from MetaPhlAn 2 / HUMAnN 2 (bioBakery 2), HMP2 and Hall from bioBakery 3. The reference databases differ, so absolute values (especially richness) are not directly comparable with the other cohorts. Only direction and effect size of CD vs. control within each cohort are compared. Confirmed in the "Code availability" section of Franzosa 2019.
- **Source:** Franzosa 2019; Truong 2015.

### 2.2 One sample per person
- **Rule:** `sample_id` is also used as `subject_id`. Time and reads stay empty, so the read limit (1.2) does not apply.

### 2.3 Centre
- **Rule:** the sub-cohort (PRISM or Validation) is the centre.

### 2.4 Genus = sum of species
- **Rule:** the genus is everything before the first `_` of the species name; species of one genus are summed.

### 2.5 Butyrate and propionate
- **Rule:** metabolite clusters HILIC-neg_Cluster_0013 (butyrate / isobutyrate) and HILIC-neg_Cluster_0003 (propionate) from Dataset 1.
- **Why:** butyrate and isobutyrate are not separated in this dataset, so they stay together as "butyrate". The values are peak intensities, not concentrations.
- **Use:** for SCFAs, Franzosa is used for propionate only. The merged butyrate/isobutyrate peak is at most an exploratory check, because isobutyrate comes from protein fermentation and can hide a drop in butyrate.

### 2.6 No further exclusions
- **Rule:** the 3 samples with very low Shannon (boxplot rule) are flagged but kept. `07_robustness.R` runs the comparison without them (9.3).

---

## 3 Hall 2017

### 3.1 Data source
- **Data:** curatedMetagenomicData, version 2021-03-31, `HallAB_2017.relative_abundance`. `01_download_cmd.R` downloads and saves the data unchanged.
- **Why:** a fixed version gives the same data on a later download. Relative abundances instead of estimated counts, as for the other cohorts.
- **Source:** Hall 2017; Pasolli 2017.
- **Changes:** until 2026-10-28 the raw file had been edited in Excel (UC removed, `control` entered). Since then it is unchanged and the selection is in the code. The taxa values are identical.

### 3.2 Controls
- **Rule:** controls have no `disease_subtype` (NA) in curatedMetagenomicData; `02_clean_hall.R` sets `control`.

### 3.3 Centre
- **Rule:** subject IDs starting with SKST are from Atlanta, all others from Boston.

### 3.4 Time in weeks
- **Rule:** time is converted from days to weeks, the unit used for HMP2.

### 3.5 Species renamed
- **Rule:** Hall uses a newer taxonomy, so species are renamed to the names used in HMP2 and Franzosa. Only cases relevant for the hypotheses or frequent ones; extend the list if needed. Columns with the same name after renaming are added up.

### 3.6 Species without a genus
- **Rule:** 23 species have no genus in the taxonomy, among them *[Eubacterium] rectale*, several *[Clostridium]* and *[Eubacterium]* species and "Firmicutes bacterium CAG:…". They are dropped when summing to genus level but kept at species level.
- **Open:** should they be assigned to a genus (e.g. *Agathobacter* for *E. rectale*, an important butyrate producer)?

### 3.7 Genus-only rows
- **Rule:** the 3 rows classified only to genus level are added to that genus.

### 3.8 SKST006_10_G102994 excluded
- **Rule:** in `exclude_samples` in `R/00_setup.R`.
- **Why:** 100 % Streptococcus, 1 genus and Shannon = 0 despite 7.3 million reads; biologically implausible. It is the only sample of this person, so person SKST006 (control) is lost.

---

## 4 HMP2 (IBDMDB)

### 4.1 Data source
- **Data:** ibdmdb.org, bioBakery 3 products: `hmp2_metadata_2018-08-20.csv`, `taxonomic_profiles_3.tsv`, `pathabundances_3.tsv`, `ecs_3.tsv`.
- **Source:** Lloyd-Price 2019; Beghini 2021.

### 4.2 Metagenomes of CD and controls only
- **Rule:** `data_type == "metagenomics"`, CD and nonIBD. nonIBD is called `control`.

### 4.3 Adults with known age only
- **Rule:** `hmp2_min_age <- 18`, based on `consent_age`. Persons without age are dropped.
- **Why:** comparability with the other cohorts; own decision of 21.09.2026, without a sensitivity analysis including all ages. Without an age it cannot be checked whether a person is an adult.
- **Effect:** 43 children dropped (Cincinnati 26, MGH Pediatrics 11, Emory 6), plus 5 persons without age (C3007, C3031, E5006, E5023, M2082; 18 samples).
- **Changes:** the selection used to be done in Excel; 4 persons without age were listed in the code without reason and C3031 (CD, 11 samples) was included. Since 2026-10-28 the rule applies to all 5. Shannon CD vs. control −0.255 (p = 0.036) before, −0.248 (p = 0.044) after.

### 4.4 Clinical values from the nearest visit
- **Rule:** HBI, calprotectin and antibiotics come from the same sample, otherwise from the nearest visit of the same person at most `max_gap` = 2 weeks away. Columns `_how` and `_gap` show where the value comes from.
- **Why:** these values are often recorded at a different visit than the stool sample. The 2 weeks are an own decision; `03_eda_clinical.R` repeats the comparison with 1 and 0 weeks.

### 4.5 MSM6J2MF removed from the feature tables
- **Why:** no enzyme data. Removed from all feature tables so they contain the same samples; it stays in `hmp2_meta.csv`.

### 4.6 CSM7KOOR excluded
- **Rule:** `hmp2_empty_samples` in `R/00_setup.R`.
- **Why:** not a single count (row sum 0), cannot be normalised.

### 4.7 Centre
- **Rule:** first letter of the Participant ID.

### 4.8 Butyrate pathways and enzymes
- **Rule:** total rows only, not the rows stratified by species. Pathways CENTFERM-PWY, PWY-5676 and PWY-5677; enzymes EC 2.3.1.19 (phosphate butyryltransferase), 2.7.2.7 (butyrate kinase) and 2.8.3.8 (acetate CoA-transferase).
- **Why:** single genes such as *but* and *buk* can be misannotated, so whole pathways are more robust. The selection is an own decision; succinate fermentation is not covered by Vital 2014.
- **Source:** Vital 2014.
- **Propionate:** EC 4.1.1.41, 4.2.1.28 and 4.2.1.54 were planned; not in the code, H3 is work in progress (0.7).

---

## 5 Clinical data

### 5.1 Calprotectin cut-offs
- **Rule:** 150 µg/g is the main cut-off for low vs. high activity; 50 and 250 µg/g are sensitivity analyses. Shown on a log scale; values of 0 are dropped in the plot.
- **Why:** there is no single accepted cut-off. Below 50 is the most sensitive range, 250 indicates large ulcers. Calprotectin is strongly right-skewed.
- **Source:** Ananthakrishnan 2023 (150); Rokkas 2018 (50); D'Haens 2012 (250).
- **Note:** in Hall only 3 of the 12 controls have a value, so the comparison there is uncertain.

### 5.2 HBI categories
- **Rule:** < 5 remission, 5–7 mild, 8–16 moderate, > 16 severe. HBI is analysed in CD only; controls should have no HBI.
- **Why:** HBI is a symptom score for CD. A score of 4 or less matches a CDAI of 150 or less, the usual definition of remission.
- **Source:** Harvey & Bradshaw 1980; Vermeire 2010.

### 5.3 HBI and calprotectin are compared, not merged
- **Rule:** samples with both values are compared, but the two scores are used separately. Test: mixed model log10(calprotectin) ~ HBI + (1 | person); check: Spearman on the median per person.
- **Why:** HBI measures symptoms, and calprotectin was not related to symptom scores in CD.
- **Source:** D'Haens 2012.
- **Changes:** 2026-11-16: before, a Spearman correlation over all samples, which treated repeated samples as independent (HMP2 p < 0.001; now p = 0.10, per person p = 0.035).

---

## 6 Cohort description

### 6.1 Age per person
- **Rule:** age is summarised per person, so that people with many samples do not count more. Hall has only age categories, so they are shown separately.

---

## 7 Diversity

### 7.1 Three filters
- **Rule:** every diversity analysis runs with 0, 5 and 10 % prevalence filter (1.3). 10 % is used for the main figures.

### 7.2 Samples with very few taxa
- **Rule:** samples among the lowest 5 % in the number of taxa are flagged, not removed. Plots show whether they stand out. Own decision.

### 7.3 Alpha diversity measures
- **Rule:** richness, Shannon and Simpson are calculated; Shannon is the main measure.
- **Why:** richness counts taxa; Shannon and Simpson also use abundances, so rare taxa count less.
- **Source:** Shannon 1948; Simpson 1949.

### 7.4 Stability within a person (Hall, HMP2)
- **Rule:** SD of Shannon for people with at least 3 samples (own decision). Bray-Curtis distance of samples of the same person vs. other persons with the same diagnosis, and between consecutive samples.
- **Source:** Bray & Curtis 1957.

### 7.5 Change over time (Hall, HMP2)
- **Rule:** mixed model with time × diagnosis interaction; `REML = FALSE` so that models can be compared; a second model with a separate slope per person. The share of people with at least 5 samples (own decision) whose diversity increases is counted.

### 7.6 Possible confounders
- **Rule:** the CD effect is checked adjusted for age (not Hall, no age), without people on antibiotics (Franzosa only), and per centre (centres with at least 3 persons of both diagnoses, own decision).

### 7.7 CLR transformation
- **Rule:** CLR with zeros replaced by half the smallest non-zero value.
- **Why:** sequencing data are compositional. CLR needs positive values; the pseudocount is a convention, not a fixed rule.
- **Source:** Aitchison 1986; Gloor 2017; Lubbe 2021.

### 7.8 Distances
- **Rule:** Aitchison distance (Euclidean on CLR) is the main distance, Bray-Curtis the sensitivity analysis. Bray-Curtis is not Euclidean, so it is shown with PCoA instead of PCA.
- **Source:** Aitchison 2000; Bray & Curtis 1957; Gower 1966.

### 7.9 PERMANOVA
- **Rule:** Hall and HMP2: main test on all samples, the diagnosis is permuted between whole persons (999 permutations), so that samples of one person stay together. Check: earliest sample per person. Adjustment for centre, the dispersion test and the run without antibiotics (Franzosa) use the earliest sample.
- **Why:** samples of one person are not exchangeable, so permuting single samples would be invalid. A different spread can also make PERMANOVA significant, hence the dispersion test. Centre and all samples together are not combined, because restricted permutations with a covariate get complicated; this is a limitation.
- **Source:** Anderson 2001; Anderson 2006; Anderson & ter Braak 2003.
- **Changes:** 2026-11-16: all-samples version added as main test; before, only the earliest sample.

### 7.10 Outliers in the ordination
- **Rule:** distance to the median above median + 3 MAD. The alpha test is repeated without them.
- **Why:** the MAD is robust against the outliers themselves; the factor 3 is a convention.
- **Source:** Leys 2013.

### 7.11 UMAP for display only
- **Rule:** UMAP is not used for tests; `min_dist` is an own choice.
- **Source:** McInnes 2018.

---

## 8 Most and least abundant taxa

### 8.1 One value per person
- **Rule:** samples are averaged per person first, so that people with many samples do not count more. 10 % filter.

### 8.2 Length and content of the lists
- **Rule:** 10 taxa or samples per cohort and group (own decision). The least abundant taxa are only those that occur in the group at all. Plots show Franzosa and HMP2; Hall stays in the CSV.
- **Why:** the lists are descriptive; the least abundant taxa often occur in only one person.

### 8.3 Literature check
- **Note:** the top-10 lists agree with the literature: Escherichia (Franzosa, HMP2) and Streptococcus, Lactobacillus, Veillonella appear only in the CD lists (Gevers 2014: Enterobacteriaceae, Veillonellaceae, Streptococcus and Lactobacillales higher in CD); Faecalibacterium is lower in CD (Sokol 2008; Joossens 2011); Akkermansia, Coprococcus and Dorea appear only in the control lists (Png 2010; Gevers 2014). Dialister in the Franzosa CD list goes against Joossens 2011. Both cohorts are dominated by Bacteroides with little Prevotella, the pattern of a Western diet (Wu 2011). Least abundant genera are mostly single persons and not interpreted. Descriptive only; tested taxa are in 12.1.
- **Changes:** 2026-11-23: added (week 7 to-do).

---

## 9 Robustness

### 9.1 Variance explained
- **Rule:** marginal R² (fixed effects) and conditional R² (fixed effects + person). `unique_r2` = how much the marginal R² falls without a variable. For Franzosa this is the usual R².
- **Why:** overlapping variables share variance, so the values do not add up; small negative values mean no contribution.
- **Source:** Nakagawa & Schielzeth 2013.

### 9.2 Disease activity
- **Rule:** CD samples with low and high activity are compared with controls, using the cut-offs of 5.1 and 5.2. A group with fewer than 5 samples is skipped. Within CD, diversity is related to the scores; values below 1 are set to 1 for the log.
- **Why:** the literature is mixed: lower diversity is also found in CD in remission (Magro 2019; Pisani 2022, IBD in remission, CD and UC together), but in HMP2 it is mainly reported for active, dysbiotic periods (Lloyd-Price 2019).
- **Source:** Magro 2019; Pisani 2022; Lloyd-Price 2019.

### 9.3 Franzosa outliers
- **Rule:** boxplot rule (1.5 × IQR) plus the 10 lowest CD samples (own decision). They are checked, not removed from the main analysis. The CD effect is compared with and without them, also as Hodges-Lehmann estimate.
- **Source:** Tukey 1977; Hodges & Lehmann 1963.

### 9.4 Reads and samples per person
- **Rule:** reads enter the models as log10(reads).
- **Why:** more reads can mean more rare taxa, so richness can depend on depth; this has to be checked in each dataset. Samples of one person vary over time, so people with more samples can look different.
- **Source:** Cattonaro 2020; Vandeputte 2021.

### 9.5 Spread within a person
- **Rule:** the SD is corrected with c4 for the number of samples, and recalculated with exactly 4 random samples per person, 1000 times (own decision).
- **Why:** the SD is too small when a person has few samples.
- **Source:** Gurland & Tripathi 1971.

### 9.6 Time trend
- **Rule:** the trend is checked for richness vs. evenness, reads and centre (almost all HMP2 controls are from M), and whether the CD effect is there from the start.

### 9.7 Hall with Boston only
- **Rule:** Hall is analysed with Boston only, and a simulation estimates which CD effect could be detected with so few controls, using the variance between and within persons and the real number of samples per person.
- **Source:** Arnold 2011.

### 9.8 HMP2 vs. Franzosa
- **Rule:** the weaker CD effects in HMP2 are explained with existing results only, no new test: cohort description (6), centre (1.6, 7.6), disease activity (9.2) and the H1/H2 effect sizes (10, 11). Genus level.
- **Why:** HMP2 shows the same direction as Franzosa, but smaller and mostly not significant effects. Three possible reasons were checked: centre, activity and the taxa effects themselves.
- **Note:** HMP2 controls are people without IBD after the first colonoscopy (screening or non-specific symptoms), described as not completely healthy; 22 of the 56 Franzosa controls come from a general population study (LifeLines-DEEP), the other 34 are PRISM non-IBD controls. In HMP2, dysbiosis occurs in episodes (about a quarter of CD samples), is only weakly tied to activity and not to medication.
- **Centre within CD:** HMP2 CD patients of centre M and C are compared (mixed model `y ~ site + (1 | person)`, Wilcoxon on the median per person as check), because the centre effect cannot be estimated in controls (1 control in C). Age, calprotectin and HBI are compared per person.
- **Why fewer genera in C:** within CD, centre C and M are compared per person for sequencing depth (log10 reads), age, disease duration (age minus age at diagnosis), bowel surgery ever, ileal involvement (Montreal L1 or L3), immunosuppressants ever and antibiotics ever (raw HMP2 metadata; Wilcoxon or Fisher). The richness model within CD is then refitted with one of these variables at a time; a variable explains part of the centre effect if the centre estimate shrinks. Exploratory (30 persons), read descriptively, no correction; biologics, disease course, BMI and smoking have too many missing values.
- **Decision:** the model with centre stays the main HMP2 model; without centre and centre M only are sensitivity analyses. Batch correction is not used for tests (1.6).
- **Source:** Lloyd-Price 2019; Franzosa 2019.
- **Changes:** 2026-11-20: added.

---

## 10 H1: dysbiosis

### 10.1 MD-index as log10
- **Rule:** log10 of the ratio of taxa up in CD to taxa down in CD.
- **Why:** Gevers 2014 does not state the base; Wei 2021 write it as natural log. The base only rescales the index (log10 = ln / 2.303); tests and p-values stay the same.
- **Source:** Gevers 2014; Wei 2021.

### 10.2 MD-index taxa
- **Rule:** the 13 taxa of Table S2I of Gevers 2014. Up: Enterobacteriaceae, Pasteurellaceae, Veillonellaceae, Fusobacteriaceae, Neisseriaceae, Micrococcaceae, Gemellaceae; down: Bacteroidales, Clostridiales, Erysipelotrichaceae, Bifidobacteriaceae, Verrucomicrobiaceae, Coriobacteriaceae. Greengenes names are mapped to the taxonomy used here: Gemellaceae → genus Gemella; Verrucomicrobiaceae → Akkermansiaceae; Coriobacteriaceae → orders Coriobacteriales and Eggerthellales (split since Greengenes). Check: the 7 taxa named in the summary of Gevers 2014.
- **Why:** the summary names only 7 taxa; Table S2I lists the 13 taxa per sample, and all go in the same direction as the CD coefficients in Table S2A (Coriobacteriaceae only q = 0.09).
- **Note:** the MD-index was derived from ileal biopsies of treatment-naive children with 16S data; controls had non-inflammatory symptoms. Gevers 2014 report that the dysbiosis is weaker in stool (AUC 0.66 vs. 0.85 in ileum). Recomputed from Table S2I (log10, no antibiotics), the median CD − control difference is 0.64 in ileum, 0.53 in rectum and 0.39 in stool.
- **Source:** Gevers 2014 (Tables S2A, S2I).
- **Changes:** 2026-11-23: 13 taxa of Table S2I as main, 7 taxa of the summary as check.

### 10.3 Zero in the numerator
- **Rule:** replaced by half of the smallest non-zero numerator of the cohort (own decision). Sensitivity: a value 10 times smaller.
- **Why:** Gevers 2014 does not describe it (checked in the full text and Table S2).

### 10.4 MD-index without filter
- **Rule:** species level, no prevalence filter, so that rare numerator taxa stay.

### 10.5 Bray-Curtis dysbiosis score
- **Rule:** median Bray-Curtis distance to the reference controls, without samples of the same person, as in Lloyd-Price 2019 (species level; in HMP2 only control samples after week 20 are the reference). With the filter of this thesis (1.3). Sensitivity analyses: no filter, all controls as reference, genus level, and Aitchison distance (CLR as in 7.7).
- **Source:** Lloyd-Price 2019; Aitchison 1986.

### 10.6 Dysbiotic samples
- **Rule:** a sample is dysbiotic above the 90th percentile of the controls. Sensitivity: 80th and 95th percentile. In Franzosa the samples are compared with Fisher. In HMP2: logistic mixed model dysbiotic ~ diagnosis + site + (1 | person) on all samples; check: Wilcoxon on the share of dysbiotic time points per person.
- **Source:** Lloyd-Price 2019.
- **Changes:** 2026-11-16: HMP2 main test changed from the Wilcoxon per person to the mixed model.

### 10.7 Taxonomy for Franzosa
- **Rule:** Franzosa has only species names, so family and order come from the HMP2 taxonomy via the genus. Franzosa "genera" that are really higher levels or missing in HMP2 are listed by hand; genera with more than one family in HMP2 are checked.

### 10.8 Models
- **Rule:** Franzosa: Wilcoxon as main test, `lm(score ~ diagnosis + site)` as check (1.11). HMP2: `lmer(score ~ diagnosis + site + (1 | person))` as main test, Wilcoxon on the median per person as check. Effects with 95 % CI.
- **Changes:** 2026-11-16: Franzosa main test changed from lm to Wilcoxon; HMP2 Wilcoxon check changed from the earliest sample to the median per person.

---

## 11 H2: fibre-associated taxa

### 11.1 No prevalence filter for chosen taxa
- **Rule:** abundances and ratios of Prevotella, Bacteroides, *F. prausnitzii* and Roseburia are not prevalence-filtered.
- **Why:** the taxa were chosen in advance.

### 11.2 Prevotella/Bacteroides in three parts
- **Rule:** (1) is Prevotella detected, (2) log10 ratio in samples with both taxa, (3) Prevotella-dominated samples (Prevotella > Bacteroides).
- **Why:** Prevotella is absent from about two thirds of the samples, and a replacement value would decide the result.
- **Tests:** (1) logistic model / GLMM, Fisher as check; (2) lm / lmer, Wilcoxon as check; (3) Franzosa: Fisher; HMP2: GLMM dom ~ diagnosis + site + (1 | person), Fisher per person (at least one dominated sample) as check. The GLMM estimate in (3) is very uncertain, because a person is almost always or never Prevotella-dominated.
- **Changes:** 2026-11-16: GLMM added as HMP2 main test for (3); before, Fisher per person only.

### 11.3 H2b main method: MaAsLin2
- **Rule:** abundances are TSS-scaled over all species beforehand, so MaAsLin2 runs with normalization NONE, transform LOG, method LM, no prevalence filter. MaAsLin2 corrects over all terms including site; the hypothesis is only about diagnosis.
- **Source:** Mallick 2021.

### 11.4 H2b check: CLR models
- **Rule:** linear models on CLR values (7.7) with the 10 % genus filter (1.3). All Roseburia species become one part, because the hypothesis is about the genus.

### 11.5 Multiple testing
- **Rule:** Holm correction over the two H2b taxa, per cohort.

### 11.6 Other butyrate producers (exploratory)
- **Rule:** E. rectale, A. hadrus, E. hallii and the genus Coprococcus are tested like H2b (MaAsLin2 main, CLR models as check, person as random effect in HMP2). Holm correction over these four taxa, separately from the two pre-specified H2b taxa. Exploratory, because they were chosen after the week 8 meeting.
- **Why:** they are butyrate producers present in both cohorts (prevalence ≥ 45 %) and are linked to fibre intake. Butyricicoccus (too rare), Subdoligranulum (Franzosa only) and R. bromii (degrades resistant starch but does not produce butyrate) are not included.
- **Source:** Louis & Flint 2009; Louis & Flint 2017; Duncan 2007; Walker 2011; Rivière 2016.
- **Changes:** 2026-11-16: added. MaAsLin2 results of 11.3 and 11.6 now with a 95 % CI (normal approximation, coefficient ± 1.96 SE).

### 11.7 Plant/animal taxa ratio (David 2014)
- **Rule:** second H2a measure: log10 [(Roseburia + E. rectale + R. bromii) / (Alistipes + Bilophila + Bacteroides)] on all samples; Roseburia, Alistipes, Bilophila and Bacteroides as genus. Zeros replaced as for the MD-index (10.3), sensitivity with a value 10 times smaller. Tests as in 1.11; without antibiotics in Franzosa only (7.6).
- **Why:** Prevotella/Bacteroides can only be tested as ratio in samples with both taxa (Franzosa 54 of 144, 11.2), which costs power. The six taxa changed in opposite directions under a plant-based vs. animal-based diet in the same controlled study (David 2014), and are rarely zero (zeros in plant / animal part: Franzosa 1 / 2 of 144 samples, HMP2 34 / 23 of 571). A log-ratio between taxa does not depend on the total microbial load (Morton 2019). Prevotella/Bacteroides stays the main H2a measure: it separates long-term carbohydrate from protein/animal-fat diets (Wu 2011), and a high ratio predicts a larger response to a high-fibre diet (Hjorth 2019).
- **Not used:** Firmicutes/Bacteroidetes, because it is not diet-specific and disputed as a marker (Magne 2020; Wei 2021).
- **Note:** all six taxa are also in the list of 12.1; chosen before looking at the ratio results.
- **Source:** David 2014; Morton 2019; Hjorth 2019; Magne 2020; Wei 2021.
- **Changes:** 2026-11-23: added after the week 8 meeting (slide 9).

---

## 12 Diet- and CD-associated taxa

### 12.1 Literature list
- **Rule:** 17 taxa linked to diet in studies of healthy people (fibre/plant or animal/fat) or reported as changed in CD/IBD, present in both cohorts. A taxon can have both links. Tested like H2b (11.3; person as random effect in HMP2; no prevalence filter, 11.1); Roseburia, Coprococcus, Prevotella, Bacteroides and Alistipes as genus, as in section 11. Holm correction over the 17 taxa per cohort. Expected direction for the diet taxa: lower (fibre/plant) or higher (animal/fat) in CD, if CD patients eat less fibre (diet proxy, see 0). Two forest plots (diet, CD) show only taxa with Holm p < 0.05 in at least one cohort; a second version shows all 17 taxa (filled points: Holm p < 0.05), so that taxa against or without the expected direction are visible too. All results are in the CSV. A descriptive plot shows how common each taxon is in CD and controls: share of samples with the taxon and median relative abundance in the samples that have it (HMP2: all samples, no test). It shows whether a CD effect comes from how often a taxon is present or from how much of it there is. Exploratory, because the list was chosen after the week 8 meeting.
- **Why:** the top-6 tables show which taxa are common, not which differ, and give no test.
- **Source:** diet: Wu 2011; De Filippo 2010; David 2014; Walker 2011; Ze 2012; Ramirez-Farias 2009; Duncan 2007; Vandeputte 2017 (Anaerostipes); Benítez-Páez 2019 (E. hallii). CD: Sokol 2008; Joossens 2011; Darfeuille-Michaud 2004; Gevers 2014; Png 2010; Lloyd-Price 2019.
- **Note:** B. wadsworthia is also linked to colitis (milk fat in mice), but only listed under diet. The coefficients of the H2b taxa are the same as in 11.3 and 11.6; only the Holm values differ.
- **Note:** Gevers 2014 report R. gnavus lower in CD in 43 stool metagenomes (text; not significant in Table S2F), against the expected direction.
- **Note:** Coprococcus is listed under CD only (lower in CD biopsies and stool metagenomes, Gevers 2014), not under diet: no study found where it rises with fibre, and inulin lowered it in people with high habitual fibre intake (Healey 2018). It stays in H2b as a butyrate producer (11.6).
- **Checked:** expected directions checked in abstracts and, where open access, full texts (2026-11-23). Well supported: F. prausnitzii, Roseburia, E. rectale, R. bromii, B. adolescentis, Prevotella, Bacteroides, Alistipes, B. wadsworthia, F. nucleatum, E. hallii (Benítez-Páez 2019). Weaker: A. hadrus (fibre effect shown for the genus Anaerostipes only, Vandeputte 2017); D. invisus (Joossens 2011 only, p = 0.04); R. gnavus (Gevers 2014 report it lower in stool); E. coli and A. muciniphila are mucosal findings (Darfeuille-Michaud 2004; Png 2010), not stool.
- **Changes:** 2026-11-19: added; replaces the top-6 table of the week 8 slides. 2026-11-23: Coprococcus moved from diet to CD; diet sources for A. hadrus and E. hallii added.

---

## 13 Disease severity within CD

### 13.1 Genera vs. calprotectin
- **Rule:** CD only. Each genus (10 % filter on all samples, 1.3) vs. log10 calprotectin (values below 1 set to 1, as in 9.2) with MaAsLin2, settings as in 11.3, centre as fixed effect, person as random effect in HMP2. BH correction over all genera; q < 0.25 counts as a hit (MaAsLin2 default). Checks: HMP2 with HBI instead of calprotectin (more persons); Franzosa without antibiotics. Agreement: share of genera with the same direction in both cohorts and Spearman correlation of the coefficients.
- **Why:** answers the week 8 question which bacteria go with more severe disease. The cohorts are analysed separately and compared side by side, not pooled (0.10): the severity measures differ (Franzosa calprotectin only) and the processing differs (2.1). BH instead of Holm because all genera are tested, not a short list; exploratory.
- **Source:** Mallick 2021.
- **Changes:** 2026-11-23: added after the week 8 meeting (slide 12).

### 13.2 MD-index vs. calprotectin
- **Rule:** CD only. Franzosa `lm(md ~ log10 calprotectin + site)`; HMP2 `lmer(md ~ log10 calprotectin + site + (1 | person))`.
- **Why:** Gevers 2014 found that the MD-index rises with disease severity (PCDAI, in biopsies of children); this checks whether it holds in stool of adults.
- **Source:** Gevers 2014.

---

## Literature

- Aitchison J (1986). The statistical analysis of compositional data. doi:10.1007/978-94-009-4109-0
- Aitchison J et al. (2000). Logratio analysis and compositional distance. Math Geol 32(3):271-275. doi:10.1023/a:1007529726302
- Ananthakrishnan AN et al. (2023). AGA clinical practice guideline on the role of biomarkers for the management of Crohn's disease. Gastroenterology 165(6):1367-1399. doi:10.1053/j.gastro.2023.09.029
- Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecol 26(1):32-46. doi:10.1111/j.1442-9993.2001.01070.pp.x
- Anderson MJ, ter Braak CJF (2003). Permutation tests for multi-factorial analysis of variance. J Stat Comput Simul 73(2):85-113. (to check in full text)
- Anderson MJ (2006). Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62(1):245-253. doi:10.1111/j.1541-0420.2005.00440.x
- Arnold BF et al. (2011). Simulation methods to estimate design power: an overview for applied research. BMC Med Res Methodol 11:94. doi:10.1186/1471-2288-11-94
- Bates D et al. (2015). Fitting linear mixed-effects models using lme4. J Stat Softw 67(1). doi:10.18637/jss.v067.i01
- Beghini F et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. eLife 10:e65088. doi:10.7554/eLife.65088
- Benítez-Páez A et al. (2019). A multi-omics approach to unraveling the microbiome-mediated effects of arabinoxylan oligosaccharides in overweight humans. mSystems 4(4):e00209-19. doi:10.1128/msystems.00209-19
- Bray JR, Curtis JT (1957). An ordination of the upland forest communities of southern Wisconsin. Ecol Monogr 27(4):325-349. doi:10.2307/1942268
- Cao Q et al. (2021). Effects of rare microbiome taxa filtering on statistical analysis. Front Microbiol 11:607325. doi:10.3389/fmicb.2020.607325
- Cattonaro F et al. (2020). Do you cov me? Effect of coverage reduction on metagenome shotgun sequencing studies. F1000Res 7:1767. doi:10.12688/f1000research.16804.4
- Darfeuille-Michaud A et al. (2004). High prevalence of adherent-invasive Escherichia coli associated with ileal mucosa in Crohn's disease. Gastroenterology 127(2):412-421. doi:10.1053/j.gastro.2004.04.061
- David LA et al. (2014). Diet rapidly and reproducibly alters the human gut microbiome. Nature 505(7484):559-563. doi:10.1038/nature12820
- De Filippo C et al. (2010). Impact of diet in shaping gut microbiota revealed by a comparative study in children from Europe and rural Africa. Proc Natl Acad Sci USA 107(33):14691-14696. doi:10.1073/pnas.1005963107
- D'Haens G et al. (2012). Fecal calprotectin is a surrogate marker for endoscopic lesions in inflammatory bowel disease. Inflamm Bowel Dis 18(12):2218-2224. doi:10.1002/ibd.22917
- Duncan SH et al. (2007). Reduced dietary intake of carbohydrates by obese subjects results in decreased concentrations of butyrate and butyrate-producing bacteria in feces. Appl Environ Microbiol 73(4):1073-1078.
- Franzosa EA et al. (2019). Gut microbiome structure and metabolic activity in inflammatory bowel disease. Nat Microbiol 4(2):293-305. doi:10.1038/s41564-018-0306-4
- Gevers D et al. (2014). The treatment-naive microbiome in new-onset Crohn's disease. Cell Host Microbe 15(3):382-392. doi:10.1016/j.chom.2014.02.005
- Gloor GB et al. (2017). Microbiome datasets are compositional: and this is not optional. Front Microbiol 8:2224. doi:10.3389/fmicb.2017.02224
- Gower JC (1966). Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53(3/4):325-338. doi:10.2307/2333639
- Gurland J, Tripathi RC (1971). A simple approximation for unbiased estimation of the standard deviation. Am Stat 25(4):30-32. doi:10.1080/00031305.1971.10477279
- Hall AB et al. (2017). A novel Ruminococcus gnavus clade enriched in inflammatory bowel disease patients. Genome Med 9:103. doi:10.1186/s13073-017-0490-5
- Harvey RF, Bradshaw JM (1980). A simple index of Crohn's-disease activity. Lancet 315(8167):514. doi:10.1016/s0140-6736(80)92767-1
- Healey G et al. (2018). Habitual dietary fibre intake influences gut microbiota response to an inulin-type fructan prebiotic: a randomised, double-blind, placebo-controlled, cross-over, human intervention study. Br J Nutr 119(2):176-189. doi:10.1017/s0007114517003440
- Hillmann B et al. (2018). Evaluating the information content of shallow shotgun metagenomics. mSystems 3(6):e00069-18. doi:10.1128/msystems.00069-18
- Hjorth MF et al. (2019). Prevotella-to-Bacteroides ratio predicts body weight and fat loss success on 24-week diets varying in macronutrient composition and dietary fiber: results from a post-hoc analysis. Int J Obes 43(1):149-157. doi:10.1038/s41366-018-0093-2
- Hodges JL, Lehmann EL (1963). Estimates of location based on rank tests. Ann Math Stat 34(2):598-611. doi:10.1214/aoms/1177704172
- Johnson WE et al. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8(1):118-127. doi:10.1093/biostatistics/kxj037
- Joossens M et al. (2011). Dysbiosis of the faecal microbiota in patients with Crohn's disease and their unaffected relatives. Gut 60(5):631-637. doi:10.1136/gut.2010.223263
- Kuznetsova A et al. (2017). lmerTest package: tests in linear mixed effects models. J Stat Softw 82(13). doi:10.18637/jss.v082.i13
- Leys C et al. (2013). Detecting outliers: do not use standard deviation around the mean, use absolute deviation around the median. J Exp Soc Psychol 49(4):764-766. doi:10.1016/j.jesp.2013.03.013
- Lloyd-Price J et al. (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. Nature 569(7758):655-662. doi:10.1038/s41586-019-1237-9
- Louis P, Flint HJ (2009). Diversity, metabolism and microbial ecology of butyrate-producing bacteria from the human large intestine. FEMS Microbiol Lett 294(1):1-8. doi:10.1111/j.1574-6968.2009.01514.x
- Louis P, Flint HJ (2017). Formation of propionate and butyrate by the human colonic microbiota. Environ Microbiol 19(1):29-41. doi:10.1111/1462-2920.13589
- Lubbe S et al. (2021). Comparison of zero replacement strategies for compositional data with large numbers of zeros. Chemometr Intell Lab Syst. doi:10.1016/j.chemolab.2021.104248
- Lumley T et al. (2002). The importance of the normality assumption in large public health data sets. Annu Rev Public Health 23:151-169. (to check in full text)
- Ma S et al. (2022). Population structure discovery in meta-analyzed microbial communities and inflammatory bowel disease using MMUPHin. Genome Biol 23:208. doi:10.1186/s13059-022-02753-4
- Magne F et al. (2020). The Firmicutes/Bacteroidetes ratio: a relevant marker of gut dysbiosis in obese patients? Nutrients 12(5):1474. doi:10.3390/nu12051474
- Magro DO et al. (2019). Remission in Crohn's disease is accompanied by alterations in the gut microbiota and mucins production. Sci Rep 9:13263. doi:10.1038/s41598-019-49893-5
- Mallick H et al. (2021). Multivariable association discovery in population-scale meta-omics studies. PLoS Comput Biol 17(11):e1009442. doi:10.1371/journal.pcbi.1009442
- Mann HB, Whitney DR (1947). On a test of whether one of two random variables is stochastically larger than the other. Ann Math Stat 18(1):50-60. doi:10.1214/aoms/1177730491
- McInnes L et al. (2018). UMAP: uniform manifold approximation and projection. J Open Source Softw 3(29):861. doi:10.21105/joss.00861
- Morton JT et al. (2019). Establishing microbial composition measurement standards with reference frames. Nat Commun 10:2719. doi:10.1038/s41467-019-10656-5
- Nakagawa S, Schielzeth H (2013). A general and simple method for obtaining R2 from generalized linear mixed-effects models. Methods Ecol Evol 4(2):133-142. doi:10.1111/j.2041-210x.2012.00261.x
- Nearing JT et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. Nat Commun 13:342. doi:10.1038/s41467-022-28034-z
- Nygaard V et al. (2016). Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. Biostatistics 17(1):29-39. doi:10.1093/biostatistics/kxv027
- Pasolli E et al. (2017). Accessible, curated metagenomic data through ExperimentHub. Nat Methods 14(11):1023-1024. doi:10.1038/nmeth.4468
- Pisani A et al. (2022). Dysbiosis in the gut microbiota in patients with inflammatory bowel disease during remission. Microbiol Spectr 10(3):e0061622. doi:10.1128/spectrum.00616-22
- Png CW et al. (2010). Mucolytic bacteria with increased prevalence in IBD mucosa augment in vitro utilization of mucin by other bacteria. Am J Gastroenterol 105(11):2420-2428. doi:10.1038/ajg.2010.281
- Ramirez-Farias C et al. (2009). Effect of inulin on the human gut microbiota: stimulation of Bifidobacterium adolescentis and Faecalibacterium prausnitzii. Br J Nutr 101(4):541-550. doi:10.1017/s0007114508019880
- Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res 43(7):e47. doi:10.1093/nar/gkv007
- Rivière A et al. (2016). Bifidobacteria and butyrate-producing colon bacteria: importance and strategies for their stimulation in the human gut. Front Microbiol 7:979. doi:10.3389/fmicb.2016.00979
- Rokkas T et al. (2018). Fecal calprotectin in assessing inflammatory bowel disease endoscopic activity: a diagnostic accuracy meta-analysis. J Gastrointestin Liver Dis 27(3):299-306. doi:10.15403/jgld.2014.1121.273.pti
- Satterthwaite FE (1946). An approximate distribution of estimates of variance components. Biometrics Bull 2(6):110-114. (to check in full text)
- Schielzeth H et al. (2020). Robustness of linear mixed-effects models to violations of distributional assumptions. Methods Ecol Evol 11(9):1141-1152. (to check in full text)
- Shannon CE (1948). A mathematical theory of communication. Bell Syst Tech J 27(3):379-423. doi:10.1002/j.1538-7305.1948.tb01338.x
- Simpson EH (1949). Measurement of diversity. Nature 163(4148):688. doi:10.1038/163688a0
- Sokol H et al. (2008). Faecalibacterium prausnitzii is an anti-inflammatory commensal bacterium identified by gut microbiota analysis of Crohn disease patients. Proc Natl Acad Sci USA 105(43):16731-16736. doi:10.1073/pnas.0804812105
- Truong DT et al. (2015). MetaPhlAn2 for enhanced metagenomic taxonomic profiling. Nat Methods 12(10):902-903. doi:10.1038/nmeth.3589
- Tukey JW (1977). Exploratory data analysis. Addison-Wesley.
- Vandeputte D et al. (2017). Prebiotic inulin-type fructans induce specific changes in the human gut microbiota. Gut 66(11):1968-1974. doi:10.1136/gutjnl-2016-313271
- Vandeputte D et al. (2021). Temporal variability in quantitative human gut microbiome profiles and implications for clinical research. Nat Commun 12:6740. doi:10.1038/s41467-021-27098-7
- Vermeire S et al. (2010). Correlation between the Crohn's disease activity and Harvey-Bradshaw indices in assessing Crohn's disease severity. Clin Gastroenterol Hepatol 8(4):357-363. doi:10.1016/j.cgh.2010.01.001
- Vital M et al. (2014). Revealing the bacterial butyrate synthesis pathways by analyzing (meta)genomic data. mBio 5(2):e00889-14. doi:10.1128/mbio.00889-14
- Walker AW et al. (2011). Dominant and diet-responsive groups of bacteria within the human colonic microbiota. ISME J 5(2):220-230.
- Wei S et al. (2021). Determining gut microbial dysbiosis: a review of applied indexes for assessment of intestinal microbiota imbalances. Appl Environ Microbiol 87(11):e00395-21. doi:10.1128/AEM.00395-21
- Wilcoxon F (1945). Individual comparisons by ranking methods. Biometrics Bull 1(6):80-83. doi:10.2307/3001968
- Wong B (2011). Points of view: Color blindness. Nat Methods 8(6):441. doi:10.1038/nmeth.1618
- Wu GD et al. (2011). Linking long-term dietary patterns with gut microbial enterotypes. Science 334(6052):105-108. doi:10.1126/science.1208344
- Ze X et al. (2012). Ruminococcus bromii is a keystone species for the degradation of resistant starch in the human colon. ISME J 6(8):1535-1543. doi:10.1038/ismej.2012.4
