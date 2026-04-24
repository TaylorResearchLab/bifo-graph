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
- **Disease / condition:** Congenital Heart Disease
- **Cohort size:** a patient cohort
- **Input genes (seeds):** 1276
- **Pathways scored:** 2130 total; 2015 with valid null distribution
- **Significant pathways (q < 0.05):** 265
- **Analysis date:** 2026-04-23
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
| **contributing_seeds** | Seed genes whose PPR signal flowed through the graph and contributed to this pathway's score. Includes genes connected through mechanistic edges, not just direct pathway members. |
| **seed_members** | Seed genes that are direct annotated members of this pathway in the source database. This is what Fisher enrichment uses. When contributing_seeds >> seed_members, the pathway was recovered through graph propagation rather than direct overlap. |

**rank and null_z are complementary, not interchangeable.** Rank reflects absolute propagated signal; null_z reflects whether signal is specifically attributable to this pathway's gene membership rather than generic graph topology. null_z reflects separation from a graph-specific null and is not a directly comparable effect size across different analyses.

---

## Top 50 Pathways by BIFO Score (well-calibrated pathways only)

*Full results in pathway_results_summary.tsv*

| rank | pathway_name                                                                 | source | n_members | degree_norm  | null_z | empirical_q | in_reference | contributing_seeds             | seed_members                                                                                                                                                                                                                                                                                                                           |
| ---- | ---------------------------------------------------------------------------- | ------ | --------- | ------------ | ------ | ----------- | ------------ | ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 43   | WP_CILIOPATHIES                                                              | MSIGDB | 170       | 8.509283e-06 | 41.187 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 45   | REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS                             | MSIGDB | 171       | 8.438848e-06 | 41.048 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 47   | Transmembrane Transport                                                      | GO     | 199       | 8.390948e-06 | 36.433 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ABCA7;ABCB4;ABCB6;ABCC6;ANO3;ANO6;AQP4;ASIC2;ATP6V0A1;ATP7B;CACNA1C;CACNA1E;CACNA1G;CACNA1S;CHRNA3;CHRNE;CLCN1;CNGA3;CTNS;GJB4;KCNK18;KCNV2;NPR2;OCA2;PIEZO1;PIEZO2;PKD1L1;PRKAR2A;RYR1;SCN5A;SCN9A;SLC1A4;SLC24A4;SLC24A5;SLC25A13;SLC37A4;SLC39A4;SLC4A1;SLC52A2;SLC5A2;SLC6A19;SLC7A7;SLC7A8;SLC7A9;SLC9A9;SORT1;TCIRG1;TRPC6;TRPM6 |
| 51   | Regulation of Cell Shape                                                     | GO     | 161       | 8.126826e-06 | 36.868 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ACTA1;ACTG2;ARPC1B;CAPZA2;CLIP1;DLG1;DNAH5;DNAH9;DVL1;ESPN;FN1;KRT10;KRT14;KRT16;KRT18;KRT3;KRT6A;MAPT;MYBPC3;MYH11;MYH2;MYH6;NEB;PALS1;ROBO3;SLC4A1;SPTA1;SPTBN1;SYNE1                                                                                                                                                                |
| 56   | Synaptic Transmission                                                        | GO     | 90        | 7.604403e-06 | 37.382 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ASIC2;CACNA1E;CACNA1G;CHRNA3;CHRNE;CLCN1;CNTNAP1;DRD4;GNB5;GRM6;GRM7;SLC1A4;SNCA;TH;TRH                                                                                                                                                                                                                                                |
| 74   | HALLMARK_ADIPOGENESIS                                                        | MSIGDB | 178       | 6.703463e-06 | 32.406 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 75   | Apoptosis                                                                    | GO     | 238       | 6.676506e-06 | 32.266 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ATM;BCL2L12;CCAR1;CERKL;DIABLO;HUWE1;IL1B;NBN;PERP;PGAP2;PRF1;SNCA;TNFRSF10B;TNFRSF11B;TNFRSF13B;TNFRSF1A;UBA3;WWOX                                                                                                                                                                                                                    |
| 77   | HALLMARK_INFLAMMATORY_RESPONSE                                               | MSIGDB | 196       | 6.598837e-06 | 34.862 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 79   | REACTOME_SLC_MEDIATED_TRANSMEMBRANE_TRANSPORT                                | MSIGDB | 246       | 6.521001e-06 | 33.454 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 81   | REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING                          | MSIGDB | 110       | 6.460943e-06 | 34.780 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 87   | REACTOME_SENSORY_PROCESSING_OF_SOUND                                         | MSIGDB | 73        | 6.293131e-06 | 32.703 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 90   | WP_PANCREATIC_ADENOCARCINOMA_PATHWAY                                         | MSIGDB | 88        | 6.181623e-06 | 33.604 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 92   | REACTOME_TRANSPORT_OF_INORGANIC_CATIONS_ANIONS_AND_AMINO_ACIDS_OLIGOPEPTIDES | MSIGDB | 105       | 6.118499e-06 | 29.642 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 93   | HALLMARK_MYOGENESIS                                                          | MSIGDB | 195       | 6.081837e-06 | 29.837 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 94   | REACTOME_BINDING_AND_UPTAKE_OF_LIGANDS_BY_SCAVENGER_RECEPTORS                | MSIGDB | 98        | 6.062775e-06 | 33.105 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 95   | REACTOME_ANCHORING_OF_THE_BASAL_BODY_TO_THE_PLASMA_MEMBRANE                  | MSIGDB | 94        | 6.049225e-06 | 29.972 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 97   | REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE                                | MSIGDB | 101       | 5.969839e-06 | 32.653 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 99   | REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION                               | MSIGDB | 101       | 5.960940e-06 | 32.536 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 103  | HALLMARK_SPERMATOGENESIS                                                     | MSIGDB | 131       | 5.918830e-06 | 31.625 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 104  | Lipid Metabolism                                                             | GO     | 111       | 5.880109e-06 | 26.350 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ADIPOQ;ALG3;ALG6;BCO1;CYP2C9;GBA2;GLB1;GNE;HEXA;LEPR;LMF1;LPL;PLCB3;PPT1;PSAP                                                                                                                                                                                                                                                          |
| 105  | WP_BREAST_CANCER_PATHWAY                                                     | MSIGDB | 152       | 5.861110e-06 | 29.958 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 107  | REACTOME_INTEGRATION_OF_ENERGY_METABOLISM                                    | MSIGDB | 108       | 5.789466e-06 | 29.254 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 108  | HALLMARK_FATTY_ACID_METABOLISM                                               | MSIGDB | 147       | 5.785287e-06 | 27.674 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 110  | Myogenesis                                                                   | GO     | 36        | 5.771170e-06 | 24.429 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ACTA1;CAV3;DMPK;EVC;EVC2;FLNC;KLHL40;MYH6;NEB;PLN;RYR1                                                                                                                                                                                                                                                                                 |
| 113  | REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA                         | MSIGDB | 113       | 5.648445e-06 | 30.442 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 114  | Receptor Mediated Endocytosis                                                | GO     | 34        | 5.613381e-06 | 25.417 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | CAV3;CUBN;LDLR;LRP1;LRP5;LRPAP1;MRC1;PPT1                                                                                                                                                                                                                                                                                              |
| 117  | WP_IL18_SIGNALING_PATHWAY                                                    | MSIGDB | 269       | 5.482210e-06 | 26.107 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 118  | HALLMARK_COMPLEMENT                                                          | MSIGDB | 197       | 5.478526e-06 | 28.034 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 121  | Anabolism                                                                    | GO     | 61        | 5.359378e-06 | 23.000 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | COQ2;COQ6;COQ8B;DHFR;DHODH;DUOX2;GNE;GPHN;GSS;HEXA;MMAB;MMACHC;MOCS1;MVK;TH;TPO                                                                                                                                                                                                                                                        |
| 122  | REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2                         | MSIGDB | 129       | 5.292003e-06 | 27.888 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 123  | Carbohydrate Metabolism                                                      | GO     | 80        | 5.283013e-06 | 21.116 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 | ACE;ALDH2;B3GLCT;DHCR7;FBP1;G6PC3;GALM;GALT;GLB1;GNPTAB;GNPTG;GUSB;IDUA;INPPL1;PMM2;SI;SLC37A4;SLC3A1;SLC5A2;VCAN                                                                                                                                                                                                                      |
| 124  | WP_PRADERWILLI_AND_ANGELMAN_SYNDROME                                         | MSIGDB | 61        | 5.279384e-06 | 27.449 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 126  | REACTOME_TRNA_PROCESSING_IN_THE_NUCLEUS                                      | MSIGDB | 56        | 5.206749e-06 | 27.254 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 127  | HALLMARK_G2M_CHECKPOINT                                                      | MSIGDB | 194       | 5.163151e-06 | 25.403 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 128  | REACTOME_METABOLISM_OF_STEROIDS                                              | MSIGDB | 151       | 5.152835e-06 | 26.556 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 130  | PID_RB_1PATHWAY                                                              | MSIGDB | 64        | 5.118894e-06 | 26.323 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 131  | WP_G1_TO_S_CELL_CYCLE_CONTROL                                                | MSIGDB | 64        | 5.118707e-06 | 27.407 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 133  | WP_HUMAN_THYROID_STIMULATING_HORMONE_TSH_SIGNALING_PATHWAY                   | MSIGDB | 66        | 5.085734e-06 | 26.145 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 134  | REACTOME_ASSEMBLY_OF_THE_PRE_REPLICATIVE_COMPLEX                             | MSIGDB | 65        | 5.063811e-06 | 27.422 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 135  | WP_DNA_DAMAGE_RESPONSE                                                       | MSIGDB | 67        | 5.024728e-06 | 26.623 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 136  | PID_P75_NTR_PATHWAY                                                          | MSIGDB | 67        | 5.016498e-06 | 29.249 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 137  | PID_TELOMERASE_PATHWAY                                                       | MSIGDB | 67        | 4.998091e-06 | 25.266 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 138  | WP_SMALL_CELL_LUNG_CANCER                                                    | MSIGDB | 95        | 4.992818e-06 | 24.770 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 141  | WP_MELANOMA                                                                  | MSIGDB | 68        | 4.971735e-06 | 26.299 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 143  | WP_GPCRS_OTHER                                                               | MSIGDB | 92        | 4.910387e-06 | 25.964 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 144  | REACTOME_TOLL_LIKE_RECEPTOR_CASCADES                                         | MSIGDB | 151       | 4.892571e-06 | 26.648 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 145  | WP_NONSMALL_CELL_LUNG_CANCER                                                 | MSIGDB | 71        | 4.867979e-06 | 26.229 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 146  | WP_MECP2_AND_ASSOCIATED_RETT_SYNDROME                                        | MSIGDB | 71        | 4.858346e-06 | 27.512 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 147  | PID_E2F_PATHWAY                                                              | MSIGDB | 72        | 4.822490e-06 | 24.640 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |
| 148  | REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION                       | MSIGDB | 258       | 4.798234e-06 | 23.828 | 0.0081      | NA           | CEP290;PKD1L1;ASPM;CC2D2A;ROS1 |                                                                                                                                                                                                                                                                                                                                        |

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

- "What biological processes are enriched in Congenital Heart Disease?"
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

### Seed genes (input to BIFO) — 1276 resolved

These genes were used as the starting point for signal propagation. They are typically genes carrying rare variants in the cohort, or genes of interest in the experimental system.

C0162832, C0242988, C0376571, C0525037, C0537026, C0598034, C0694889, C0694897, C0812211, C0812248, C0812265, C0812267, C0812271, C0812278, C0812281, C0812290, C0812295, C0812304, C0812382, C0879290, C0879391, C0879392, C0919427, C0919481, C0919508, C0919524, C1325331, C1332001, C1332010, C1332029, C1332121, C1332375, C1332409, C1332419, C1332427, C1332724, C1332770, C1332774, C1332775, C1332803, C1332811, C1332829, C1333216, C1333217, C1333219, C1333225, C1333254, C1333356, C1333357, C1333358, C1333523, C1333532, C1333542, C1333544, C1333568, C1333569, C1333690, C1333705, C1333713, C1333935, C1334081, C1334082, C1334099, C1334112, C1334291, C1334292, C1334339, C1334469, C1334474, C1334481, C1334513, C1334523, C1334532, C1334538, C1334862, C1334895, C1335083, C1335088, C1335225, C1335232, C1335265, C1335270, C1335280, C1335607, C1335609, C1335620, C1335844, C1335872, C1335874, C1336563, C1336575, C1336642, C1336644, C1336645, C1336659, C1336674, C1336679, C1336688, C1336938, C1337033, C1337090, C1337098, C1337106, C1363984, C1364081, C1364508, C1366370, C1366475, C1366496, C1366499, C1366509, C1366515, C1366529, C1366623, C1366631, C1366757, C1366857, C1367472, C1367702, C1368803, C1412058, C1412061, C1412062, C1412066, C1412070, C1412071, C1412077, C1412081, C1412082, C1412107, C1412108, C1412110, C1412111, C1412114, C1412132, C1412136, C1412167, C1412179, C1412192, C1412208, C1412209, C1412225, C1412262, C1412277, C1412279, C1412284, C1412291, C1412304, C1412306, C1412332, C1412335, C1412340, C1412355, C1412364, C1412367, C1412374, C1412375, C1412384, C1412387, C1412388, C1412390, C1412401, C1412471, C1412494, C1412546, C1412553, C1412582, C1412591, C1412625, C1412629, C1412637, C1412641, C1412685, C1412686, C1412689, C1412691, C1412741, C1412744, C1412749, C1412759, C1412760, C1412770, C1412807, C1412833, C1412933, C1412946, C1412954, C1412981, C1413005, C1413014, C1413019, C1413024, C1413029, C1413036, C1413053, C1413056, C1413057, C1413059, C1413061, C1413062
... and 1076 more (see --seeds file).

### Reference pathway set

These pathways were pre-specified as biologically relevant before the analysis was run. They are marked TRUE in the in_reference column.

No reference set provided.

### Analysis parameters

- Scores file: pathway_scores_standard.csv
- Total pathways scored: 2130
- Pathways with valid null distribution: 2015
- Null model: degree-preserving membership rewiring, N = 1,000 permutations
- Multiple testing correction: Benjamini-Hochberg
- Pathway size filter: minimum 8 members, maximum 300 members
- Pathway sources: MSigDB canonical collections, WikiPathways, Reactome, GO

---
*Generated by summarize_results.py on 2026-04-23*
