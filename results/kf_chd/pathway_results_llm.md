## Instructions for AI

You are helping a scientist interpret the **biological meaning** of pathway analysis results from a computational tool called BIFO (Biological Information Flow Ontology). Your role is to:

- Help the user understand which biological processes are enriched in their data
- Explain what the pathway names mean in plain language
- Connect the results to the disease, condition, or experimental system the user is studying
- Suggest follow-up questions, experiments, or analyses that might be warranted
- Answer questions about the input data (what genes were used, what the analysis did)

This session is intended for biological interpretation of the results. Unless the user asks you to evaluate the method, please focus this conversation on the biology: what the pathways mean, how they relate to the disease or experimental system, and what follow-up questions or experiments might be worth pursuing. If the user wants to discuss the method, answer their questions, then return focus to the biology.

*Note: these outputs are intended for exploratory biological interpretation. They are not part of any quantitative evaluation and are not used for scoring, ranking, or statistical inference.*

The user may be a biologist, clinician, or computational scientist. Do not assume prior knowledge of graph algorithms or bioinformatics. The user may have questions specific to their own experimental system, organism, tissue type, or clinical context, not necessarily the disease listed below.

---

## Analysis Summary: KF-CHD

- **Cohort:** KF-CHD
- **Disease / condition:** congenital heart disease
- **Cohort size:** 697 probands
- **Input genes (seeds):** 1276 (1287 in input file; 1276 resolved to graph CUIs)
- **Pathways scored:** 2130 total; 2015 with valid null distribution
- **Significant pathways (q < 0.05):** 265
- **Analysis date:** 2026-04-22
- **Source file:** pathway_scores_standard.csv

---

## What This Analysis Does

BIFO propagates signal from a set of input genes through a large biomedical knowledge graph connecting genes, proteins, pathways, diseases, and biological processes. Only biologically meaningful relationships are used; statistical associations and text-mining edges are excluded.

Signal flows from the input genes outward and accumulates at pathway nodes. Pathways strongly connected to the input genes through biological relationships receive higher scores. Scores are then compared to a statistical null model (1,000 random rewirings of gene-pathway membership) to identify pathways that score higher than expected by chance given the graph topology.

This approach is particularly useful when the input gene list is large or heterogeneous because it concentrates distributed biological signal rather than requiring strong direct overlap between input genes and pathway members.

---

## How to Read the Results Table

| Column | What it means |
|--------|---------------|
| **rank** | Position by BIFO score (1 = highest). Lower rank = more signal from input genes. |
| **pathway_name** | Name of the biological pathway or process. |
| **source** | Database the pathway comes from (MSIGDB, WikiPathways, Reactome, GO). |
| **n_members** | Number of genes in the pathway. |
| **degree_norm** | BIFO score: propagated signal at the pathway node, normalised by pathway size. |
| **null_z** | Standard deviations above expected score vs. 1,000 random rewirings. Higher = more specific enrichment. NaN = test not valid. |
| **empirical_q** | BH-corrected p-value. q < 0.05 = statistically significant. NaN = not tested. |
| **in_reference** | TRUE if in a pre-specified biologically relevant reference set. |

**rank and null_z are complementary, not interchangeable.** Rank reflects absolute propagated signal; null_z reflects whether signal is specifically attributable to this pathway's gene membership rather than generic graph topology. null_z reflects separation from a graph-specific null and is not a directly comparable effect size across different analyses.

---

## Top 50 Pathways by BIFO Score (well-calibrated pathways only)

*Full results in pathway_results_summary.tsv*

| rank | pathway_name                                                                 | source | n_members | degree_norm  | null_z | empirical_q | in_reference |
| ---- | ---------------------------------------------------------------------------- | ------ | --------- | ------------ | ------ | ----------- | ------------ |
| 43   | WP_CILIOPATHIES                                                              | MSIGDB | 170       | 8.509283e-06 | 41.187 | 0.0081      | True         |
| 45   | REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS                             | MSIGDB | 171       | 8.438848e-06 | 41.048 | 0.0081      | False        |
| 47   | Transmembrane Transport                                                      | GO     | 199       | 8.390948e-06 | 36.433 | 0.0081      | False        |
| 51   | Regulation of Cell Shape                                                     | GO     | 161       | 8.126826e-06 | 36.868 | 0.0081      | False        |
| 56   | Synaptic Transmission                                                        | GO     | 90        | 7.604403e-06 | 37.382 | 0.0081      | False        |
| 74   | HALLMARK_ADIPOGENESIS                                                        | MSIGDB | 178       | 6.703463e-06 | 32.406 | 0.0081      | False        |
| 75   | Apoptosis                                                                    | GO     | 238       | 6.676506e-06 | 32.266 | 0.0081      | False        |
| 77   | HALLMARK_INFLAMMATORY_RESPONSE                                               | MSIGDB | 196       | 6.598837e-06 | 34.862 | 0.0081      | False        |
| 79   | REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT                                | MSIGDB | 246       | 6.521001e-06 | 33.454 | 0.0081      | False        |
| 81   | REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING                          | MSIGDB | 110       | 6.460943e-06 | 34.780 | 0.0081      | False        |
| 87   | REACTOME_SENSORY_PROCESSING_OF_SOUND                                         | MSIGDB | 73        | 6.293131e-06 | 32.703 | 0.0081      | False        |
| 90   | WP_PANCREATIC_ADENOCARCINOMA_PATHWAY                                         | MSIGDB | 88        | 6.181623e-06 | 33.604 | 0.0081      | False        |
| 92   | REACTOME_TRANSPORT_OF_INORGANIC_CATIONS_ANIONS_AND_AMINO_ACIDS_OLIGOPEPTIDES | MSIGDB | 105       | 6.118499e-06 | 29.642 | 0.0081      | False        |
| 93   | HALLMARK_MYOGENESIS                                                          | MSIGDB | 195       | 6.081837e-06 | 29.837 | 0.0081      | False        |
| 94   | REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS                | MSIGDB | 98        | 6.062775e-06 | 33.105 | 0.0081      | False        |
| 95   | REACTOME_ANCHORING_OF_THE_BASAL_BODY_TO_THE_PLASMA_MEMBRANE                  | MSIGDB | 94        | 6.049225e-06 | 29.972 | 0.0081      | False        |
| 97   | REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE                                | MSIGDB | 101       | 5.969839e-06 | 32.653 | 0.0081      | False        |
| 99   | REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION                               | MSIGDB | 101       | 5.960940e-06 | 32.536 | 0.0081      | False        |
| 103  | HALLMARK_SPERMATOGENESIS                                                     | MSIGDB | 131       | 5.918830e-06 | 31.625 | 0.0081      | False        |
| 104  | Lipid Metabolism                                                             | GO     | 111       | 5.880109e-06 | 26.350 | 0.0081      | False        |
| 105  | WP_BREAST_CANCER_PATHWAY                                                     | MSIGDB | 152       | 5.861110e-06 | 29.958 | 0.0081      | False        |
| 107  | REACTOME_INTEGRATION_OF_ENERGY_METABOLISM                                    | MSIGDB | 108       | 5.789466e-06 | 29.254 | 0.0081      | False        |
| 108  | HALLMARK_FATTY_ACID_METABOLISM                                               | MSIGDB | 147       | 5.785287e-06 | 27.674 | 0.0081      | False        |
| 110  | Myogenesis                                                                   | GO     | 36        | 5.771170e-06 | 24.429 | 0.0081      | False        |
| 113  | REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA                         | MSIGDB | 113       | 5.648445e-06 | 30.442 | 0.0081      | False        |
| 114  | Receptor Mediated Endocytosis                                                | GO     | 34        | 5.613381e-06 | 25.417 | 0.0081      | False        |
| 117  | WP_IL18_SIGNALING_PATHWAY                                                    | MSIGDB | 269       | 5.482210e-06 | 26.107 | 0.0081      | False        |
| 118  | HALLMARK_COMPLEMENT                                                          | MSIGDB | 197       | 5.478526e-06 | 28.034 | 0.0081      | False        |
| 121  | Anabolism                                                                    | GO     | 61        | 5.359378e-06 | 23.000 | 0.0081      | False        |
| 122  | REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2                         | MSIGDB | 129       | 5.292003e-06 | 27.888 | 0.0081      | False        |
| 123  | Carbohydrate Metabolism                                                      | GO     | 80        | 5.283013e-06 | 21.116 | 0.0081      | False        |
| 124  | WP_PRADERWILLI_AND_ANGELMAN_SYNDROME                                         | MSIGDB | 61        | 5.279384e-06 | 27.449 | 0.0081      | False        |
| 126  | REACTOME_TRNA_PROCESSING_IN_THE_NUCLEUS                                      | MSIGDB | 56        | 5.206749e-06 | 27.254 | 0.0081      | False        |
| 127  | HALLMARK_G2M_CHECKPOINT                                                      | MSIGDB | 194       | 5.163151e-06 | 25.403 | 0.0081      | False        |
| 128  | REACTOME_METABOLISM_OF_STEROIDS                                              | MSIGDB | 151       | 5.152835e-06 | 26.556 | 0.0081      | False        |
| 130  | PID_RB_1PATHWAY                                                              | MSIGDB | 64        | 5.118894e-06 | 26.323 | 0.0081      | False        |
| 131  | WP_G1_TO_S_CELL_CYCLE_CONTROL                                                | MSIGDB | 64        | 5.118707e-06 | 27.407 | 0.0081      | False        |
| 133  | WP_HUMAN_THYROID_STIMULATING_HORMONE_TSH_SIGNALING_PATHWAY                   | MSIGDB | 66        | 5.085734e-06 | 26.145 | 0.0081      | False        |
| 134  | REACTOME_ASSEMBLY_OF_THE_PRE_REPLICATIVE_COMPLEX                             | MSIGDB | 65        | 5.063811e-06 | 27.422 | 0.0081      | False        |
| 135  | WP_DNA_DAMAGE_RESPONSE                                                       | MSIGDB | 67        | 5.024728e-06 | 26.623 | 0.0081      | False        |
| 136  | PID_P75_NTR_PATHWAY                                                          | MSIGDB | 67        | 5.016498e-06 | 29.249 | 0.0081      | False        |
| 137  | PID_TELOMERASE_PATHWAY                                                       | MSIGDB | 67        | 4.998091e-06 | 25.266 | 0.0081      | False        |
| 138  | WP_SMALL_CELL_LUNG_CANCER                                                    | MSIGDB | 95        | 4.992818e-06 | 24.770 | 0.0081      | False        |
| 141  | WP_MELANOMA                                                                  | MSIGDB | 68        | 4.971735e-06 | 26.299 | 0.0081      | False        |
| 143  | WP_GPCRS_OTHER                                                               | MSIGDB | 92        | 4.910387e-06 | 25.964 | 0.0081      | False        |
| 144  | REACTOME_TOLL_LIKE_RECEPTOR_CASCADES                                         | MSIGDB | 151       | 4.892571e-06 | 26.648 | 0.0081      | False        |
| 145  | WP_NONSMALL_CELL_LUNG_CANCER                                                 | MSIGDB | 71        | 4.867979e-06 | 26.229 | 0.0081      | False        |
| 146  | WP_MECP2_AND_ASSOCIATED_RETT_SYNDROME                                        | MSIGDB | 71        | 4.858346e-06 | 27.512 | 0.0081      | False        |
| 147  | PID_E2F_PATHWAY                                                              | MSIGDB | 72        | 4.822490e-06 | 24.640 | 0.0081      | False        |
| 148  | REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION                       | MSIGDB | 258       | 4.798234e-06 | 23.828 | 0.0081      | False        |

---

## Key Findings

**Significant pathways (q < 0.05):** 265

- WP_CILIOPATHIES (rank 43, null_z = 41.187, q = 0.0081)
- REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS (rank 45, null_z = 41.048, q = 0.0081)
- Transmembrane Transport (rank 47, null_z = 36.433, q = 0.0081)
- Regulation of Cell Shape (rank 51, null_z = 36.868, q = 0.0081)
- Synaptic Transmission (rank 56, null_z = 37.382, q = 0.0081)
- HALLMARK_ADIPOGENESIS (rank 74, null_z = 32.406, q = 0.0081)
- Apoptosis (rank 75, null_z = 32.266, q = 0.0081)
- HALLMARK_INFLAMMATORY_RESPONSE (rank 77, null_z = 34.862, q = 0.0081)
- REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT (rank 79, null_z = 33.454, q = 0.0081)
- REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING (rank 81, null_z = 34.780, q = 0.0081)


**Top 5 by null_z (strongest statistical enrichment):**


---

## Suggested Questions to Ask

- "What biological processes are enriched in congenital heart disease?"
- "What do the top-ranked pathways have in common biologically?"
- "Which pathways are both highly ranked and statistically significant?"
- "Are any of these pathways relevant to [your tissue / cell type / organism]?"
- "What genes from my input list are members of the top pathways?"
- "How might these pathway findings inform follow-up experiments?"
- "What is the biological function of [pathway name]?"
- "What were the input genes used to generate these results?"
- "What database did pathway [X] come from, and how was it defined?"
- "Are any of these pathways relevant to [specific disease or condition]?"

---

## Input Data Context

### Seed genes (input to BIFO) — 1276 resolved (1287 in input file; 1276 resolved to graph CUIs)

These genes were used as the starting point for signal propagation. They are typically genes carrying rare variants in the cohort, or genes of interest in the experimental system.

AAAS	4, AAR2	1, AARS2	3, ABCA1	1, ABCA12	2, ABCA13	5, ABCA2	1, ABCA3	1, ABCA7	4, ABCB11	1, ABCB4	1, ABCB6	1, ABCC2	1, ABCC6	2, ABCC8	1, ACADM	5, ACADS	3, ACADVL	6, ACAN	1, ACAT1	1, ACE	4, ACP2	1, ACSF3	1, ACTA1	1, ACTG2	2, ACTN3	1, ADA	1, ADA2	1, ADAD1	1, ADAM22	1, ADAMTS13	1, ADAMTS18	1, ADAMTS2	1, ADAMTS3	2, ADAT3	1, ADCY10	3, ADCY6	2, ADGRB2	1, ADGRV1	1, ADIPOQ	2, ADRA2B	1, ADSL	1, AGBL2	2, AGL	4, AGPS	2, AGT	1, AGTR1	1, AHCY	1, AHI1	1, AIRE	2, AK2	1, ALAD	1, ALB	1, ALDH1A3	1, ALDH2	2, ALDH6A1	1, ALDH7A1	3, ALG1	1, ALG3	1, ALG6	2, ALG8	1, ALMS1	1, ALPL	2, ALS2	1, AMACR	1, AMBN	1, AMER3	1, AMH	1, AMPD1	1, AMPD2	1, AMT	1, ANGPTL3	1, ANKRD11	1, ANKRD26	2, ANKRD36	1, ANO3	1, ANO6	2, ANTXR2	2, ANXA11	1, APC	4, APOB	4, APPL1	1, APRT	2, AQP4	2, ARHGAP11A	1, ARHGEF10L	1, ARHGEF11	1, ARHGEF12	1, ARL13B	1, ARMC9	1, ARPC1B	1, ARSA	4, ARV1	1, ASCC1	1, ASCC3	1, ASIC2	1, ASL	1, ASPA	4, ASPM	4, ASS1	1, ATG16L1	1, ATG9B	1, ATIC	1, ATL3	1, ATM	1, ATP11B	1, ATP13A2	1, ATP1A2	1, ATP2A1	1, ATP2B2	1, ATP6V0A1	1, ATP6V0A2	1, ATP6V0A4	1, ATP7B	5, ATP8A1	1, ATP8B3	1, ATXN2	1, AUH	1, AVIL	1, B3GLCT	1, B9D1	1, BAZ1A	1, BAZ2B	1, BBS10	1, BBS12	1, BBS4	2, BBS9	1, BCAM	1, BCKDHA	1, BCKDHB	1, BCKDK	1, BCL2L12	1, BCO1	2, BDP1	3, BLK	1, BLM	1, BMP6	1, BMPR1B	1, BMS1	2, BRAT1	1, BRCA1	1, BRCA2	4, BRIP1	1, BTBD2	1, BTD	6, C12orf57	2, C1QBP	1, C1QTNF5	1, C6	4, C7	4, C8B	1, C9	2, CABP1	1, CACNA1A	1, CACNA1B	1, CACNA1C	2, CACNA1E	1, CACNA1G	2, CACNA1H	1, CACNA1S	3, CACNA2D3	1, CACNA2D4	1, CANT1	1, CAPN10	1, CAPN3	1, CAPZA2	1, CARMIL2	1, CASQ1	1, CAT	1, CAV3	1, CC2D2A	5, CCAR1	1, CCBE1	1, CCDC107	2, CCDC39	1, CCT7	1, CD46	2, CD96	1, CDAN1	1, CDC45	1, CDH23	2, CDH4	2, CDH6	1, CDIN1	1, CDK13	1, CDKN2A	1, CDT1	2, CEACAM16	1, CEL	1, CELA2A	1, CELSR3	1, CENPJ	1, CEP104	2, CEP152	2, CEP250	2, CEP290	7, CEP41	2, CEP63	3, CERKL	4, CERS1	1
... and 1087 more (see --seeds file).

### Reference pathway set

These pathways were pre-specified as biologically relevant before the analysis was run. They are marked TRUE in the in_reference column.

MSIGDB:M211, MSIGDB:M219, MSIGDB:M27439, MSIGDB:M27440, MSIGDB:M27471, MSIGDB:M27472, MSIGDB:M27479, MSIGDB:M27480, MSIGDB:M27481, MSIGDB:M27497, MSIGDB:M39675, MSIGDB:M39706, MSIGDB:M39734, MSIGDB:M39826, MSIGDB:M39880, MSIGDB:M5919

### Analysis parameters

- Scores file: pathway_scores_standard.csv
- Total pathways scored: 2130
- Pathways with valid null distribution: 2015
- Null model: degree-preserving membership rewiring, N = 1,000 permutations
- Multiple testing correction: Benjamini-Hochberg
- Pathway size filter: minimum 8 members, maximum 300 members
- Pathway sources: MSigDB canonical collections, WikiPathways, Reactome, GO

---
*Generated by summarize_results.py on 2026-04-22*
