# Analysis & Decision Log - "From Taxa to Function"
*Luana da Silva do Aido - Stand 06.10.26 Updated continuously*

## 1. Overview
- **Topic:** Diet-related microbiome signatures in Crohn's disease (CD vs. healthy controls)
- **Target submission:** 07.01.2027 (official deadline 11.01.2027)
- **Supervisor:** Dr. Tugçe Bilgin Sonay
- **Datasets:** HMP2/IBDMDB, Franzosa/PRISM, HallAB

## 2. Phase plan

| Phase | Period | Content |
|---|---|---|
| Preparation | until 13.09. | Write exposé; define research questions, hypotheses and study aim; read background literature |
| Preprocessing | 14.09.–08.10. | Data search and download; data cleaning; EDA; PCA and UMAP; batch correction |
| Results | 12.10.–22.11. | Alpha/beta diversity; MDI, Bray-Curtis distance score; P/B ratio, comparison of *Faecalibacterium* and *Roseburia*, counting butyrate-associated bacteria; compare metabolite numbers with own results, predicted pathway analysis if feasible within time; visualization; answer questions and hypotheses |
| Writing | 23.11.–27.12. | Background; methods; results; discussion (limitations, future recommendations); introduction and abstract |
| Layout & Admin | 28.12.–07.01. | Figure/table/abbreviation lists, appendix, bibliography; own proofreading, then have someone else look over it; incorporate feedback, AI declaration, statement of authorship, final formatting check |


## 3. Decision Log
 
| ID | Date | Question | Decision | Rationale | Status |
|---|---|---|---|---|---|
| D01 | 20.08. | Scope: describe only, or compare? | Case-control comparison, CD vs. healthy | Fits the assignment; SQ1 (MDI) is a group comparison by definition | meeting 20.08. | decided |
| D02 | 20.08. | Order of work | Clarify scope/research question first, then data and analyses | Advice from the meeting: don't dive too early and too deep into the data | meeting 20.08. | decided |
| D03 | 26.08. | Timeline | Submit 07.01.2027, 4-day buffer | Buffer before official deadline 11.01. | Zeitplan | decided |
| D04 | 18.09. | Search criteria: what to relax? | (a) Paired same-stool-sample not needed, because I compare groups, not samples. (b) "Targeted metabolomics only" dropped; any platform is fine if SCFAs are well annotated. (c) No diet metadata used (not a filter, simply not used). | Untargeted data (HMP2) has butyrate/propionate annotated with high confidence; diet is handled via microbiome proxies | W1 s4 | decided |
| D05 | 18.09. | Which cohort for measured butyrate? | HMP2. Butyrate and propionate as relative peak intensities (not mM), compared within a cohort only | Only cohort with reliable butyrate (HMDB0000039) and propionate (HMDB0000237); acetate missing | W1 s8, s14 | decided |
| D06 | 18.09./21.09. | Role of Franzosa for SCFA | Propionate only. The merged peak "butyrate / isobutyrate*" at most as exploratory check | Butyrate and isobutyrate not separated, no HMDB ID; isobutyrate comes from protein fermentation, so one merged peak can hide a butyrate drop | W1 s14 | decided |
| D07 | 21.09. | Age groups | Adults (≥18) only; no sensitivity analysis with all ages | Comparability; my decision | project doc 21.09. | decided |
| D08 | 21.09. | Which candidates to exclude? | See exclusion log (section 4) | Reasons listed there | W1 s7 | decided |
| D09 | 22.09. | Regions goal (2 NA + 2 EU + 2 Asia) | Dropped | Only US/NL cohorts remained. Also on the question list for Tugce, so mention it to her | W1 s4 | decided |
| D10 | 21.–22.09. | What is the main result of SQ3? | **H3b (enzymes, EC level) = main analysis.** H3a (butyrate + propionate measured) = exploratory/validation. Enzymes and metabolites come from the same samples, so check whether they match (butyrate enzymes vs. butyrate in HMP2, propionate enzymes vs. propionate in Franzosa). A null result stays interpretable ("not detectable at this sample size") | Butyrate rests on 14 adult controls, so it cannot be the main hypothesis. Franzosa has EC tables but no pathway table, so enzyme level keeps H3b comparable across both cohorts | W1 s13, s15, s18 | proposal (Q1+Q3 for Tugce, answer not in my notes) |
| D11 | 22.09./29.09. | Process raw reads of IBDomics + Sheba myself (MetaPhlAn + HUMAnN)? | Not now. First finish the pipeline with the cohorts I have; look for more data only if time is left | Only extends SQ1/SQ2/H3b, not H3a; needs cluster and ≈1.5–2 TB. W1 s19 rated it "rather not" | W1 s19, W2 s7 | deferred |
| D12 | 22.09./29.09. | Are 14 controls enough for measured butyrate? Add 1000IBD (EGA)? | Keep butyrate exploratory, report effect size + CI. 1000IBD not pursued for now | EGA access takes weeks to months, may be too late. 1000IBD would add 238/255 and detect d = 0.25 | W1 s17, s19 | deferred |
| D13 | 22.09. | How to handle repeated samples in HMP2? | **A (main):** all samples + mixed model with random effect per person, e.g. `log10(butyrate) ~ diagnosis + age + sex + (1 \| subject)` (MaAsLin2 / lmer). **B (sensitivity):** one sample per person (first sample, or the one with a metabolome). **C (mean per person): not planned.** Franzosa: standard models (1 sample per person). Time course itself is not analysed | Samples from one person are not independent (pseudo-replication, p-values too optimistic); persons are the unit. ICC 0.58, so 141 samples ≈ 53 independent observations. I have no question about change over time | W1 s10, s17 | decided |
| D14 | 22.09. | Combine cohorts how? | Analyse each cohort separately (HMP2: mixed model; Franzosa: standard model with sub-cohort as covariate), then fixed-effect meta-analysis (weight = 1/SE²). No pooling of raw data. Measured butyrate (HMP2 only) is reported alone, not combined. **Also try pooled analysis** ("try both: separately and together"), details still to look at | Different scales/runs, different patient groups, HMP2 would dominate (≈577 vs. 144 metagenome samples) | W1 s11 + n | decided (pooled variant: open) |
| D15 | 22.09. | Marker enzymes for H3b | Butyrate: EC 2.8.3.8, 2.3.1.19, 2.7.2.7. Propionate: EC 4.1.1.41, 4.2.1.28, 4.2.1.54 | In both cohorts at enzyme level. Check that the supplement tables are not prevalence-filtered | W1 s12, s16 | decided |
| D16 | 22.09. | Which data for which question? | SQ1/SQ2: the largest dataset possible, taxa only is enough. SQ3: HMP2 + Franzosa | SQ1, SQ2, H3b need only shotgun metagenomes; only H3a needs measured metabolites | W1 s3 n | decided (from notes) |
| D17 | 29.09. | Hall 2017 (cMD) | In use as taxa-only cohort for SQ1/SQ2 (9 CD / 13 HC persons, 88 / 74 samples, longitudinal like HMP2). Same longitudinal handling as HMP2 (D13) | Species + genus available; no function or metabolome. **Note:** W1 s7 listed Hall as "no metabolomics" (excluded for H3a) and the 21.09. doc called it unsuitable, so please confirm this is a deliberate change | W2 s3–4, s7 n | decided (confirm) |
| D18 | 29.09. | Other cMD cohorts | Removed: Nielsen 2014, LifeLinesDeep 2016, VilaAV 2018 (+1 more, 4 in total) | Datasets only have "IBD", CD not separable (original papers talk about CD, tables don't separate it). VilaAV also has no healthy controls | W2 s3, s6 | decided |
| D19 | 22.–29.09. | Data handling | Download locally. Keep raw and cleaned tables as csv. Keep only CD and healthy, only adults. Structure each project into 2–3 data frames: Meta, Taxa, Pathways and/or Enzymes. Start a GitHub repo | Formats differ between projects (mostly long instead of wide, several IDs per person) | W2 s2–3, s7 | decided |
| D20 | 29.09. | Where to start? | Start with **one cohort**, build the pipeline, reuse it for the others. Simplest first: Franzosa. R or Python does not matter (I use R) | Cleaner than building three things at once | W2 s7 n | decided |
| D21 | 29.09. | Order of first analyses | Variability and diversity first (basis for the other analyses, like mean/SD). First taxa, then metabolites. Use PCA and UMAP and compare them. ML only if more power is needed at the end | Diversity is the basis for everything else | W2 s7 n | decided |

## 4. Analysis Log
 
### Preparation (until 13.09.)
- Exposé written.
- Research questions, hypotheses, study aim defined (SQ1 dysbiosis/MDI, SQ2 taxa, SQ3 function).
- Background literature read (e.g. Gevers 2014, Sokol 2008, Lloyd-Price 2019, Franzosa 2019, Wu 2011, Wirbel/Thomas 2019).
- Original 9-cohort dataset cannot be used: no metabolomics, so H3a is impossible.
- Decisions: D01–D03.

### Preprocessing (14.09.–08.10.)
 
#### 14.–18.09. – Data search
- **Goal:** shotgun metagenome + metabolome abundance tables, CD vs. HC, adults, SCFA covered, openly downloadable.
- **Approach (4 routes):** (1) literature incl. data-availability statements; (2) ENA: CD metagenome projects, does the same study have metabolomics?; (3) Metabolomics Workbench: CD studies, same cohort with shotgun?; (4) independent cross-check with Claude (Opus 5), which found no new literature beyond my list.
- **Funnel:** ≈50 candidates → ≈13 with both omics → 8 verified in detail → 2 fulfil every criterion (HMP2, Franzosa). Table: `Datensatz_Tracker_CD_Thesis_v3.xlsx`.
- **Problems:** (1) hardly any public data with both omics as abundances (often one is on request or behind controlled access); (2) where metabolome abundances exist, SCFAs are not usable (not measured, or butyrate/isobutyrate not separated).
- **Decisions:** D04–D08.

#### 21.09. – Case numbers and selection
| | persons CD / HC | samples CD / HC |
|---|---|---|
| HMP2, taxa + function | 35 / 14 | 358 / 219 (≈12–14 per person) |
| HMP2, with metabolome | 23 / 14 | 80 / 61 (≈4 per person) |
| Franzosa, PRISM (USA) | 68 / 34 | |
| Franzosa, NL validation | 20 / 22 | |
| Franzosa, total | 88 / 56 | 88 / 56 |
| Propionate measured (HMP2 + Franzosa) | 109 / 70 | |
 
- Probe run, mixed model, HMP2 adults, 141 samples: butyrate CD vs. HC factor 1.11 (95% CI 0.71–1.70), p = 0.65, ICC 0.58. Preview only, unadjusted, not normalized.
- Detectable effect at 80% power (Cohen's d): HMP2 one sample per person 0.95 · HMP2 mixed model 0.79 · Franzosa 0.48 · 1000IBD (238/255) 0.25.
- Franzosa supplement gives ready tables: Dataset 2 (metabolites), 4 (species), 6 (enzymes, EC). Only the `_Results.txt` files from Metabolomics Workbench are needed, not the raw zips (105–352 GB).

#### 22.09. – Meeting (Week 1)
- Presented search, exclusion log, two-cohort data basis, longitudinal options, cohort combination, hypothesis changes, limitations, Options A/B/C.
- **Decisions written in the slides:** D09, D13, D14, D15, D16; **proposed:** D10; open for Tugce: more data? cluster and storage? (see section 6).

#### 22.–29.09. – Cleaning and structuring (Week 2)
- Searched for more data, found little. Also checked GitHub repos: they often talk about CD but only include IBD.
- Downloaded locally (also used an R package for some cMD datasets, not all applicable).
- Cleaned: kept CD and HC only, adults only.
- Structure per cohort:
| Source | Meta | Species | Genus | EC | Pathway | Metabolite/SCFA |
|---|---|---|---|---|---|---|
| Franzosa 2019 | ✓ | ✓ | – | ✓ | – | ✓ (butyrate/propionate) |
| HMP2 (bioBakery3) | ✓ | ✓ | ✓ | ✓ | ✓ | not in table yet (ST000923 to download) |
| Hall 2017 (cMD) | ✓ | ✓ | ✓ | – | – | – |
 
- Numbers after cleaning, persons CD / HC: Franzosa 88/56 · HMP2 35/14 · Hall 9/13 · **total 132 / 83**. Samples CD / HC: Franzosa 88/56 · HMP2 358/219 · Hall 88/74 · **total 534 / 349**.
- Age (range, median, mean): Franzosa CD 21–77, 42, 43.0 · Franzosa HC 22–82, 32, 37.0 · HMP2 CD 18–76, 36, 38.4 · HMP2 HC 21–69, 47.5, 45.3. Hall: adults only, no further detail. Origin: Franzosa USA + NL, HMP2 USA multi-center, Hall USA.
- Observation: HC are younger than CD in Franzosa and older in HMP2, so age is a covariate (D13).
- **Problems:** (1) cMD cohorts declared as "IBD" (D18); (2) I don't know where to look for more abundance data; (3) cleaning takes longer than expected (long vs. wide format, multiple IDs per section).
- **Found:** something is not correct in my "clean" dataset, so look at it again before EDA.
- **Decisions:** D17, D18, D19.

#### 29.09. – Meeting (Week 2)
- **Decisions:** D20, D21 (from the notes on slide 7). D11/D12 parked until the pipeline exists.
#### Still to do in this phase
- Re-check the cleaned tables (error found, see above)
- Download HMP2 metabolomics (ST000923, HILIC-neg) and enzymes; confirm Franzosa Datasets 2, 4, 6 are complete
- EDA, then PCA and UMAP (compare), batch correction
- Plan the statistics: models, normalization, covariates, FDR correction, sensitivity analysis
- *(For every step: write down the case numbers before and after each filter.)*


### Results (from 12.10.)
*(empty)*
 
### Writing (from 23.11.)
*(empty)*
 
### Layout & Admin (from 28.12.)
*(empty)*

 ## 6. Open items and checks
 
- **Re-check cleaned data:** something is wrong (W2 s7).
- **Supplement completeness:** supplement tables sometimes show only the most common features. Check that 100% is there (W1 s12).
- **Units:** relative abundance or CPM (copies per million)? Normalize accordingly (W1 s12).
- **Answers from Tugce** to W1 s19 Q1–Q5 (butyrate with one cohort, more data, H3b as main result, cluster and storage, within-project comparison first). Only the notes above are recorded. Fill in the rest.
- **Cluster and storage:** cleaning is local for now. Decide later whether the cluster is needed (W1 s19 Q4).
- **Pooled analysis variant:** how exactly to run "separately and together" (D14).
- **Where to find abundance data** if more is needed (W2 s6/s7 question).
- **Idea:** with a gene list I could calculate pathways, and look at gene-expression/transcriptomics data. Is there any available? (W1 s12 n)
- **Idea:** look at the different SCFA products from different diet patterns and draw a conclusion (W1 s14 n). Fits the Discussion.
- **Unclear note** on W1 s18: "How big is the big data set". Please clarify what you meant.
- **Possible detail from the plan:** with longitudinal data I could also analyse changes within a person (W1 s10 n). Currently not planned (D13).
## 7. Limitations to state in the thesis (W1 s18)
 
- HMP2 and Franzosa come from the same lab context, so they are not fully independent replicates.
- No acetate in either open dataset.
- Measured butyrate rests on 14 adult controls.
- Faecal SCFA reflect fermentation output (≈90–95% absorbed), not fibre intake.
- Scarcity of open, CD-specific SCFA data is itself a finding.
- Propionate is not specific for fibre (also from protein, lactate, sugars; Reichardt 2014), so it supports butyrate and does not replace it.
- Different processing versions: Franzosa was profiled with bioBakery 2 (MetaPhlAn 2 / HUMAnN 2), HMP2 and Hall with bioBakery 3. Reference databases differ, so absolute values (especially richness) are not comparable across cohorts; only direction and effect size within each cohort are compared (D14). Confirmed in the "Code availability" section of Franzosa 2019 (KneadData, MetaPhlAn2, HUMAnN2).
 

