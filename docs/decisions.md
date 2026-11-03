# Decision log

What was decided about data and methods, why, and where it is in the code.
Add new decisions below under the matching cohort or method, with date.
All rules that apply to several scripts are variables in `R/00_setup.R`.

## Overview: from raw data to analysis

As of 2026-10-28, genus level, main analysis (10 % prevalence filter).

| Cohort | Raw data | After selection (metadata file) | In the analysis | Persons (CD / control) |
|---|---|---|---|---|
| Franzosa 2019 | 220 samples | 144 (without UC) | 144 | 88 / 56 (1 sample per person) |
| Hall 2017 | 259 samples, 32 persons | 162 (without UC) | 158 | 9 / 12 |
| HMP2 | 1638 metagenome samples, 130 persons | 578 (CD + controls, adults) | 571 | 30 / 14 |

## Rules for all cohorts

### CD and controls only, no UC
- **Rule:** samples with ulcerative colitis (UC) are removed.
- **Why:** the thesis compares Crohn's disease with healthy controls.
- **Code:** `R/02_clean_franzosa.R`, `R/02_clean_hall.R`, `R/02_clean_hmp2.R` (write the metadata files `*_meta.csv`).
- **Effect:** Franzosa −76 samples, Hall −97 samples, HMP2 −459 samples.
- **Date:** from the start. Done by hand in Excel until 2026-10-28, in code since then.

### Minimum sequencing depth: 1 million reads
- **Rule:** samples with fewer than `min_reads` = 1,000,000 reads are excluded.
- **Why:** with very few reads only the 2–3 most abundant genera are found. Example: Hall p9216_mo2 has 48,368 reads and 2 genera, the other samples of the same person have 41–49 genera. This distorts richness and Shannon.
- **Source of the read count:** Hall `number_reads` (curatedMetagenomicData), HMP2 `reads_filtered` (after QC and removal of human reads). Franzosa has no read counts, so the rule does not apply there.
- **Code:** `min_reads`, `read_depth()` and `low_read_samples()` in `R/00_setup.R`, applied in `06`, `07` and `09`.
- **Effect:** Hall −3 samples (p9216_mo2, p9281_mo6, SKST037_3_G103026), HMP2 −6 samples. No person is lost completely.
- **Date:** introduced in `06`/`07` before the clean-up, in `09` since 2026-10-28.

### Prevalence filter
- **Rule:** a genus is kept only if it occurs in at least 10 % of the samples of a cohort (`thr_10`). 0 % and 5 % run as sensitivity analysis.
- **Why:** very rare genera are mostly noise and make the CLR transformation unstable.
- **Code:** `thresholds` and `main_thr` in `R/06_diversity.R`.

### Baseline
- **Rule:** cross-sectional comparisons (alpha diversity, PERMANOVA) use the earliest sample per person. Longitudinal analyses use all samples with a mixed model (person as random effect).
- **Why:** several samples of the same person are not independent.
- **Code:** `baseline` in `R/06_diversity.R`.

### Batch effects (centre)
- **Rule:** in tests, the centre (`site`) is included as covariate. Corrected data (ComBat) are used for figures only (PCA/UMAP), and only if every centre has both diagnoses.
- **Why:** correcting with diagnosis in the model and then testing makes p-values too optimistic. In Hall, Atlanta has only controls. A correction would separate artificially there, so Hall stays uncorrected.
- **Centre:** Franzosa: PRISM or Validation. Hall: SKST = Atlanta, otherwise Boston. HMP2: first letter of the Participant ID.
- **Code:** batch section in `R/06_diversity.R` (`use_combat`).

### Randomness
- **Rule:** fixed seed `seed <- 42`. It is additionally set directly before permutation tests, UMAP and simulations.
- **Why:** same results on every run. Two runs of all scripts give identical files (checked 2026-10-28).
- **Code:** `R/00_setup.R`, `R/06_diversity.R`, `R/08_robustness.R`.

## Franzosa 2019

- **Data source:** Supplementary Datasets 1, 2, 4 and 6 of Franzosa et al. 2019 (Nature Microbiology).
- **Processing:** species and enzymes were generated with MetaPhlAn 2 / HUMAnN 2 (bioBakery 2), HMP2 and Hall with bioBakery 3 (MetaPhlAn 3 / HUMAnN 3). The reference databases differ, so absolute values (especially richness) are not directly comparable with the other cohorts. Only direction and effect size of CD vs. control within each cohort are compared (D14). Confirmed in the "Code availability" section of Franzosa et al. 2019 (KneadData, MetaPhlAn2, HUMAnN2).
- **No further exclusions.** The 3 samples with very low Shannon (boxplot rule) are flagged but kept in the analysis. `R/08_robustness.R` additionally runs the comparison without them.
- **Butyrate and propionate:** the metabolite clusters HILIC-neg_Cluster_0013 (butyrate / isobutyrate) and HILIC-neg_Cluster_0003 (propionate) from Dataset 1. Butyrate and isobutyrate are not separated. The values are peak intensities, not concentrations.

## Hall 2017

- **Data source:** curatedMetagenomicData, version 2021-03-31, `HallAB_2017.relative_abundance`. `R/01_download_cmd.R` downloads the data and saves them unchanged.
- **Controls:** in curatedMetagenomicData, controls have no `disease_subtype` (NA). `R/02_clean_hall.R` sets `control` instead.
- **SKST006_10_G102994 excluded** (`exclude_samples`): the sample consists of 100 % Streptococcus, has only 1 genus and Shannon = 0, despite 7.3 million reads. This is biologically implausible. It is the only sample of this person, so person SKST006 (control) is lost.
- **Genus level, open point:** 23 species have no genus in the curatedMetagenomicData taxonomy. They include *[Eubacterium] rectale*, several *[Clostridium]* and *[Eubacterium]* species and "Firmicutes bacterium CAG:…". They are dropped when summing to genus level but kept at species level. *E. rectale* is an important butyrate producer. **Still to decide:** should they be assigned to a genus (e.g. *Agathobacter* for *E. rectale*)?
- **Date:** until 2026-10-28 the raw file had been edited in Excel (UC removed, `control` entered). Since then it is unchanged, and the selection is in the code. The taxa values are identical.

## HMP2 (IBDMDB)

- **Data source:** ibdmdb.org, bioBakery 3 products: `hmp2_metadata_2018-08-20.csv`, `taxonomic_profiles_3.tsv`, `pathabundances_3.tsv`, `ecs_3.tsv`.
- **Metagenome samples only** (`data_type == "metagenomics"`).
- **Adults with known age only** (`hmp2_min_age <- 18`, `consent_age`):
  - **Why:** comparability; own decision of 21.09.2026, without a sensitivity analysis including all age groups (D07 in `docs/analysis_log.md`). Persons without age are dropped because it cannot be checked whether they are adults.
  - **Effect:** 43 children are dropped (Cincinnati 26, MGH Pediatrics 11, Emory 6), plus 5 persons without age (C3007, C3031, E5006, E5023, M2082; 18 samples).
  - **Code:** `R/02_clean_hmp2.R`.
  - **Date:** the age selection used to be done in Excel. The 4 persons without age were listed in the code without a documented reason. C3031 (CD, 11 samples) was included despite missing age. Since 2026-10-28 the rule applies to all 5. The main result hardly changes: Shannon CD vs. control −0.255 (p = 0.036) before, −0.248 (p = 0.044) after.
- **nonIBD is called `control`.**
- **MSM6J2MF excluded:** there are no EC data for this sample. It is removed from all feature tables so that they contain the same samples; it stays in `hmp2_meta.csv` (`R/02_clean_hmp2.R`).
- **CSM7KOOR excluded** (`hmp2_empty_samples`): the sample has not a single count (row sum 0) and cannot be normalised.
- **Butyrate function:** only the total rows are used, not the rows stratified by species:
  - pathways CENTFERM-PWY, PWY-5676 and PWY-5677;
  - enzymes EC 2.3.1.19 (phosphate butyryltransferase), 2.7.2.7 (butyrate kinase) and 2.8.3.8 (acetate CoA-transferase).
  - List in `R/02_clean_hmp2.R`.
