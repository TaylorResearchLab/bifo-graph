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

## Analysis Summary: KF-NBL

- **Cohort:** KF-NBL
- **Disease / condition:** neuroblastoma
- **Cohort size:** 460 probands
- **Input genes (seeds):** 1395 (1406 in input file; 1395 resolved to graph CUIs)
- **Pathways scored:** 2196 total; 2196 with valid null distribution
- **Significant pathways (q < 0.05):** 296
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

| rank | pathway_name                                                                      | source | n_members | degree_norm  | null_z | empirical_q | in_reference |
| ---- | --------------------------------------------------------------------------------- | ------ | --------- | ------------ | ------ | ----------- | ------------ |
| 1    | Ionophore activity                                                                | GO     | 186       | 4.691683e-06 | 17.824 | 0.0143      | False        |
| 2    | Transmembrane Transport                                                           | GO     | 197       | 4.367070e-06 | 16.293 | 0.0143      | False        |
| 3    | WP_CILIOPATHIES                                                                   | MSIGDB | 169       | 4.240607e-06 | 18.365 | 0.0143      | True         |
| 4    | Regulation of Cell Shape                                                          | GO     | 161       | 3.636787e-06 | 13.396 | 0.0143      | False        |
| 5    | DNA Repair                                                                        | GO     | 158       | 3.557523e-06 | 11.744 | 0.0143      | False        |
| 6    | REACTOME_DISEASES_OF_METABOLISM                                                   | MSIGDB | 242       | 3.256396e-06 | 13.146 | 0.0143      | False        |
| 7    | Blood coagulation                                                                 | GO     | 31        | 3.250208e-06 | 13.749 | 0.0143      | False        |
| 8    | Skeletal Development                                                              | GO     | 52        | 2.862661e-06 | 14.158 | 0.0143      | False        |
| 9    | Myogenesis                                                                        | GO     | 37        | 2.829276e-06 | 10.409 | 0.0143      | False        |
| 10   | WP_AMINO_ACID_METABOLISM                                                          | MSIGDB | 86        | 2.706027e-06 | 11.245 | 0.0143      | False        |
| 11   | Heparin Binding                                                                   | GO     | 31        | 2.703321e-06 | 14.826 | 0.0143      | False        |
| 12   | Carbohydrate Metabolism                                                           | GO     | 80        | 2.683632e-06 | 8.971  | 0.0143      | False        |
| 13   | PID_CONE_PATHWAY                                                                  | MSIGDB | 21        | 2.672056e-06 | 13.821 | 0.0143      | False        |
| 14   | Vision                                                                            | GO     | 12        | 2.635712e-06 | 10.053 | 0.0143      | False        |
| 15   | Xenobiotic Metabolism                                                             | GO     | 41        | 2.608968e-06 | 9.215  | 0.0143      | False        |
| 16   | Nucleotide Metabolism                                                             | GO     | 102       | 2.544961e-06 | 8.459  | 0.0143      | False        |
| 17   | WP_JOUBERT_SYNDROME                                                               | MSIGDB | 75        | 2.536478e-06 | 9.949  | 0.0143      | False        |
| 18   | REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS                                  | MSIGDB | 171       | 2.466766e-06 | 9.519  | 0.0143      | False        |
| 19   | Homeostasis                                                                       | GO     | 94        | 2.440438e-06 | 8.378  | 0.0143      | False        |
| 20   | Amino Acid Biosynthesis                                                           | GO     | 26        | 2.438118e-06 | 9.355  | 0.0143      | False        |
| 21   | Eye Development                                                                   | GO     | 14        | 2.418142e-06 | 8.320  | 0.0143      | False        |
| 22   | REACTOME_DISEASES_OF_CARBOHYDRATE_METABOLISM                                      | MSIGDB | 34        | 2.417253e-06 | 11.170 | 0.0143      | False        |
| 23   | Proteolysis                                                                       | GO     | 178       | 2.416673e-06 | 7.488  | 0.0143      | False        |
| 24   | Smell Perception                                                                  | GO     | 8         | 2.354748e-06 | 10.024 | 0.0212      | False        |
| 25   | WP_THYROID_HORMONES_PRODUCTION_AND_THEIR_PERIPHERAL_DOWNSTREAM_SIGNALLING_EFFECTS | MSIGDB | 92        | 2.305824e-06 | 10.181 | 0.0143      | False        |
| 26   | REACTOME_SLC_TRANSPORTER_DISORDERS                                                | MSIGDB | 97        | 2.297809e-06 | 9.660  | 0.0143      | False        |
| 27   | HALLMARK_XENOBIOTIC_METABOLISM                                                    | MSIGDB | 195       | 2.274288e-06 | 8.069  | 0.0143      | False        |
| 28   | Organogenesis                                                                     | GO     | 105       | 2.208983e-06 | 6.547  | 0.0143      | False        |
| 29   | Synaptic Transmission                                                             | GO     | 89        | 2.177547e-06 | 7.639  | 0.0143      | False        |
| 30   | Metabolism                                                                        | GO     | 55        | 2.175265e-06 | 6.623  | 0.0143      | False        |
| 31   | WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR                    | MSIGDB | 93        | 2.154154e-06 | 8.301  | 0.0143      | True         |
| 32   | REACTOME_CILIUM_ASSEMBLY                                                          | MSIGDB | 194       | 2.153895e-06 | 6.715  | 0.0143      | True         |
| 33   | Histogenic Process                                                                | GO     | 35        | 2.153200e-06 | 8.923  | 0.0143      | False        |
| 34   | REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION                                        | MSIGDB | 19        | 2.151692e-06 | 9.386  | 0.0143      | False        |
| 35   | recombinational repair                                                            | GO     | 24        | 2.135986e-06 | 7.630  | 0.0143      | False        |
| 36   | Muscle Contraction                                                                | GO     | 44        | 2.132654e-06 | 7.625  | 0.0143      | False        |
| 37   | REACTOME_ORGANELLE_BIOGENESIS_AND_MAINTENANCE                                     | MSIGDB | 264       | 2.109733e-06 | 6.086  | 0.0143      | False        |
| 38   | REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT                                           | MSIGDB | 8         | 2.086918e-06 | 9.622  | 0.0143      | False        |
| 39   | Complement Activation                                                             | GO     | 35        | 2.077666e-06 | 9.000  | 0.0143      | False        |
| 40   | REACTOME_VISUAL_PHOTOTRANSDUCTION                                                 | MSIGDB | 96        | 2.060800e-06 | 8.553  | 0.0143      | False        |
| 41   | Drug Metabolism                                                                   | GO     | 10        | 2.046679e-06 | 4.732  | 0.0212      | False        |
| 42   | WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK                                               | MSIGDB | 120       | 2.023458e-06 | 7.118  | 0.0143      | False        |
| 43   | REACTOME_SENSORY_PROCESSING_OF_SOUND_BY_OUTER_HAIR_CELLS_OF_THE_COCHLEA           | MSIGDB | 53        | 2.022892e-06 | 9.268  | 0.0143      | False        |
| 44   | REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                                        | MSIGDB | 296       | 1.988659e-06 | 4.943  | 0.0143      | False        |
| 45   | WP_COMPLEMENT_AND_COAGULATION_CASCADES                                            | MSIGDB | 57        | 1.977216e-06 | 8.357  | 0.0143      | False        |
| 46   | REACTOME_INTRINSIC_PATHWAY_OF_FIBRIN_CLOT_FORMATION                               | MSIGDB | 22        | 1.970223e-06 | 9.329  | 0.0143      | False        |
| 47   | PID_RHODOPSIN_PATHWAY                                                             | MSIGDB | 23        | 1.966937e-06 | 9.573  | 0.0143      | False        |
| 48   | REACTOME_TRNA_AMINOACYLATION                                                      | MSIGDB | 31        | 1.958255e-06 | 7.373  | 0.0143      | False        |
| 49   | REACTOME_DNA_REPAIR                                                               | MSIGDB | 280       | 1.948880e-06 | 4.727  | 0.0143      | False        |
| 50   | REACTOME_THYROXINE_BIOSYNTHESIS                                                   | MSIGDB | 10        | 1.948367e-06 | 9.327  | 0.0143      | False        |

---

## Key Findings

**Significant pathways (q < 0.05):** 296

- Ionophore activity (rank 1, null_z = 17.824, q = 0.0143)
- Transmembrane Transport (rank 2, null_z = 16.293, q = 0.0143)
- WP_CILIOPATHIES (rank 3, null_z = 18.365, q = 0.0143)
- Regulation of Cell Shape (rank 4, null_z = 13.396, q = 0.0143)
- DNA Repair (rank 5, null_z = 11.744, q = 0.0143)
- REACTOME_DISEASES_OF_METABOLISM (rank 6, null_z = 13.146, q = 0.0143)
- Blood coagulation (rank 7, null_z = 13.749, q = 0.0143)
- Skeletal Development (rank 8, null_z = 14.158, q = 0.0143)
- Myogenesis (rank 9, null_z = 10.409, q = 0.0143)
- WP_AMINO_ACID_METABOLISM (rank 10, null_z = 11.245, q = 0.0143)


**Top 5 by null_z (strongest statistical enrichment):**

- Heparin Binding (rank 11, null_z = 14.826, q = 0.0143)

---

## Suggested Questions to Ask

- "What biological processes are enriched in neuroblastoma?"
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

### Seed genes (input to BIFO) — 1395 resolved (1406 in input file; 1395 resolved to graph CUIs)

These genes were used as the starting point for signal propagation. They are typically genes carrying rare variants in the cohort, or genes of interest in the experimental system.

AARS2	1, ABCA12	1, ABCA13	3, ABCA2	1, ABCA3	1, ABCA5	1, ABCA7	1, ABCB11	8, ABCB4	1, ABCB6	3, ABCC6	5, ABHD16A	1, ACACA	1, ACAD11	1, ACAD9	2, ACADM	1, ACADSB	1, ACADVL	1, ACD	1, ACE	2, ACE2	1, ACO2	1, ACP2	1, ACSF3	2, ACTA1	1, ACTL6B	1, ACTN2	1, ACTR3B	1, ACVR2B	2, ACVRL1	1, ACY1	2, ADA	1, ADA2	1, ADAM17	1, ADAM22	1, ADAM9	1, ADAMTS13	1, ADAMTS18	1, ADAMTS19	1, ADAMTS2	2, ADAMTS3	2, ADAMTSL4	2, ADAT3	1, ADCY6	2, ADCY9	1, ADD2	1, ADGRE2	1, ADGRV1	3, ADIPOQ	1, ADNP	1, ADPRS	1, ADRA2B	1, AFF2	1, AFF3	2, AGAP1	1, AGBL1	1, AGBL5	1, AGK	1, AGPS	1, AGT	1, AGXT	1, AHCY	3, AHI1	2, AIM2	1, AK1	1, AKAP9	1, AKT1	1, ALDH18A1	2, ALDH1A3	1, ALDH1B1	1, ALDH2	4, ALDH7A1	1, ALG12	2, ALK	1, ALMS1	1, ALPL	2, AMOT	1, ANGPTL3	1, ANGPTL4	1, ANK3	1, ANKHD1	1, ANKRD12	1, ANKRD36	1, ANLN	1, ANO3	3, ANO5	1, ANXA11	1, AP5Z1	2, APBB1	1, APC	2, APTX	1, AQP2	1, AR	1, ARAF	1, ARFGEF2	2, ARHGAP21	1, ARHGAP4	1, ARHGEF10L	1, ARHGEF12	1, ARHGEF18	1, ARID2	1, ARID3B	1, ARL13B	1, ARL2	1, ARMC5	3, ARSA	2, ASAP1	1, ASCC1	1, ASL	1, ASNS	2, ASPA	2, ASPM	3, ASPSCR1	1, ASTN1	1, ATIC	1, ATL1	1, ATM	3, ATP11A	1, ATP11C	1, ATP13A2	1, ATP1A3	9, ATP1B1	1, ATP2A1	1, ATP2B2	1, ATP6V0A2	1, ATP7B	5, ATP8B1	1, ATP8B3	1, ATR	1, ATRX	5, AUH	2, AURKC	2, AUTS2	1, B3GALNT2	2, B3GAT3	1, B9D1	1, BARD1	4, BAZ1A	1, BBS10	1, BBS4	2, BBS9	1, BBX	1, BCAM	1, BCAS3	1, BCKDHA	1, BCL10	1, BCO1	1, BCOR	1, BCS1L	1, BDP1	3, BEND4	1, BEST1	2, BIRC6	1, BLM	4, BMP1	1, BMP6	1, BMPER	1, BMS1	1, BRAF	1, BRAT1	2, BRCA2	1, BRD4	2, BRF1	2, BRIP1	2, BRWD1	2, BRWD3	1, BSND	1, BTBD11	1, BTD	4, BTK	1, C12orf57	1, C19orf12	1, C3	1, C5	1, C6	1, C7	1, C8B	1, C9	1, CA1	1, CA12	2, CACNA1A	1, CACNA1F	2, CACNA1G	1, CACNA1H	1, CACNA1I	1, CACNA2D1	1, CACNA2D2	1, CACNA2D3	1, CACNA2D4	1, CACNB2	1, CAD	2, CAMK2G	1, CAMKV	1, CANT1	2, CAPN14	2, CAPN3	3, CARMIL2	1, CARS2	1, CASR	1, CAV1	1
... and 1206 more (see --seeds file).

### Reference pathway set

These pathways were pre-specified as biologically relevant before the analysis was run. They are marked TRUE in the in_reference column.

MSIGDB:M211, MSIGDB:M219, MSIGDB:M27439, MSIGDB:M27440, MSIGDB:M27471, MSIGDB:M27472, MSIGDB:M27479, MSIGDB:M27480, MSIGDB:M27481, MSIGDB:M27497, MSIGDB:M39675, MSIGDB:M39706, MSIGDB:M39734, MSIGDB:M39826, MSIGDB:M39880, MSIGDB:M5919

### Analysis parameters

- Scores file: pathway_scores_standard.csv
- Total pathways scored: 2196
- Pathways with valid null distribution: 2196
- Null model: degree-preserving membership rewiring, N = 1,000 permutations
- Multiple testing correction: Benjamini-Hochberg
- Pathway size filter: minimum 8 members, maximum 300 members
- Pathway sources: MSigDB canonical collections, WikiPathways, Reactome, GO

---
*Generated by summarize_results.py on 2026-04-22*
