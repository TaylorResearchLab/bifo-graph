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
- **Disease / condition:** Neuroblastoma
- **Cohort size:** a patient cohort
- **Input genes (seeds):** 1395
- **Pathways scored:** 2196 total; 2196 with valid null distribution
- **Significant pathways (q < 0.05):** 296
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

| rank | pathway_name                                                                      | source | n_members | degree_norm  | null_z | empirical_q | in_reference | contributing_seeds        | seed_members                                                                                                                                                                                                                                                                                                            |
| ---- | --------------------------------------------------------------------------------- | ------ | --------- | ------------ | ------ | ----------- | ------------ | ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | Ionophore activity                                                                | GO     | 186       | 4.691683e-06 | 17.824 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ABCB11;ANO3;ATP1A3;ATP7B;BEST1;BSND;CACNA1G;CFTR;CHRNE;CLCN1;CLCN7;CLIC4;CNGA3;CTNS;FXN;GJB4;GRIN2B;GRIN2D;KCNQ1;OCA2;ORAI1;PIEZO1;PKD1L1;RYR1;SCN5A;SCN9A;SLC1A4;SLC34A2;SLC37A4;SLC39A4;SLC44A4;SLC5A2;SLC6A19;SLC7A8;SLC7A9;SLC8A1;SLC9A6;SLCO1B1;SLCO1B3;TCIRG1;TRPM3;TRPM4;TRPM6;TRPV4;TUSC3;WNK1;WNK2             |
| 2    | Transmembrane Transport                                                           | GO     | 197       | 4.367070e-06 | 16.293 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ABCA7;ABCB4;ABCB6;ABCC6;ANO3;ATP7B;BEST1;BSND;CACNA1G;CHRNE;CLCN1;CLCN7;CLIC4;CNGA3;CTNS;FXN;GJB4;GRIN2B;GRIN2D;KCNQ1;LMBRD1;OCA2;ORAI1;PIEZO1;PKD1L1;RYR1;SCN5A;SCN9A;SLC1A4;SLC25A13;SLC29A1;SLC34A2;SLC37A4;SLC39A4;SLC44A4;SLC5A2;SLC6A19;SLC6A8;SLC7A8;SLC7A9;SLC8A1;SLC9A6;SLCO1B3;TCIRG1;TRPM3;TRPM4;TRPM6;TRPV4 |
| 3    | WP_CILIOPATHIES                                                                   | MSIGDB | 169       | 4.240607e-06 | 18.365 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 4    | Regulation of Cell Shape                                                          | GO     | 161       | 3.636787e-06 | 13.396 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACTA1;ANLN;CEP164;CEP57;CLIP1;DNAH5;DNAH9;DNAI1;DVL1;DYNC1H1;EFHC1;EMD;FN1;FRYL;KRT10;KRT12;KRT13;KRT14;KRT16;KRT18;KRT81;KRT9;LMNB1;MACF1;MARCKS;MYBPC3;MYH2;NEB;NEFH;PALS1;PLXNA1;ROBO3;SCRIB;SPTA1;SYNE1;SYNE2;TNNT2;TUBG1;ZMYM3                                                                                     |
| 5    | DNA Repair                                                                        | GO     | 158       | 3.557523e-06 | 11.744 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | APTX;ATM;ATR;BARD1;BRCA2;BRIP1;CRY2;DDX11;ERCC2;ERCC3;ERCC4;ERCC5;ERCC6;ERCC8;EXO1;FAN1;FANCA;FANCD2;FANCE;FANCG;FANCM;HUWE1;INO80;MCM8;MLH1;MRE11;MSH2;MSH3;MSH6;MUTYH;PAXIP1;PMS2;POLG;POLK;RAD50;RAD51C;RAD51D;RAD54B;RAD54L;RBBP8;TP53;TP53BP1;TYMS;UNG;UVSSA                                                       |
| 6    | REACTOME_DISEASES_OF_METABOLISM                                                   | MSIGDB | 242       | 3.256396e-06 | 13.146 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 7    | Blood coagulation                                                                 | GO     | 31        | 3.250208e-06 | 13.749 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ADAMTS13;F11;F7;F8;KLKB1;PRKAR1A;PROC;SERPINC1;VKORC1;VWF                                                                                                                                                                                                                                                               |
| 8    | Skeletal Development                                                              | GO     | 52        | 2.862661e-06 | 14.158 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | BMP1;BMP6;CCN6;EVC;EVC2;PAX1;PPIB                                                                                                                                                                                                                                                                                       |
| 9    | Myogenesis                                                                        | GO     | 37        | 2.829276e-06 | 10.409 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACTA1;CAV3;DMPK;EMD;EVC;EVC2;KLHL40;MYOCD;NCOA2;NEB;RYR1;SPEG;TNNT2;ZFHX3                                                                                                                                                                                                                                               |
| 10   | WP_AMINO_ACID_METABOLISM                                                          | MSIGDB | 86        | 2.706027e-06 | 11.245 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 11   | Heparin Binding                                                                   | GO     | 31        | 2.703321e-06 | 14.826 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | CCN6;COLQ;FGF10                                                                                                                                                                                                                                                                                                         |
| 12   | Carbohydrate Metabolism                                                           | GO     | 80        | 2.683632e-06 | 8.971  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACE;ALDH2;DHCR7;DPM1;GALK1;GALT;GLB1;GNPTAB;GNPTG;HK1;IDH1;IDUA;INSR;MANBA;NEU1;NGLY1;PDX1;PMM2;PTEN;SI;SLC37A4;SLC5A2;VCAN                                                                                                                                                                                             |
| 13   | PID_CONE_PATHWAY                                                                  | MSIGDB | 21        | 2.672056e-06 | 13.821 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 14   | Vision                                                                            | GO     | 12        | 2.635712e-06 | 10.053 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | BCO1;CDH23;CNGA3;HPS1;USH1G                                                                                                                                                                                                                                                                                             |
| 15   | Xenobiotic Metabolism                                                             | GO     | 41        | 2.608968e-06 | 9.215  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | CYP1B1;CYP2C19;CYP2C8;FMO3;HNF4A;SLCO1B1;SLCO1B3                                                                                                                                                                                                                                                                        |
| 16   | Nucleotide Metabolism                                                             | GO     | 102       | 2.544961e-06 | 8.459  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ABCA12;ABCA2;ABCA3;ABCA5;ABCA7;ABCB11;ABCB4;ABCB6;ABCC6;ADA;ADA2;ADAT3;AK1;APTX;ATP1A3;ATP7B;CANT1;DHODH;ENPP1;GALK1;GMDS;PCYT1A;PDE2A;PDE8B;TYMS;XDH                                                                                                                                                                   |
| 17   | WP_JOUBERT_SYNDROME                                                               | MSIGDB | 75        | 2.536478e-06 | 9.949  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 18   | REACTOME_DISORDERS_OF_TRANSMEMBRANE_TRANSPORTERS                                  | MSIGDB | 171       | 2.466766e-06 | 9.519  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 19   | Homeostasis                                                                       | GO     | 94        | 2.440438e-06 | 8.378  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACE2;AR;ASPSCR1;CPS1;ESR2;F11;F7;F8;HNF4A;NCOA2;NR3C2;PKHD1;RARA;SLC44A4;THRB;TRPM6;TRPV4;VWF;WNK1;WNK2                                                                                                                                                                                                                 |
| 20   | Amino Acid Biosynthesis                                                           | GO     | 26        | 2.438118e-06 | 9.355  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | AHCY;ASL;ASNS;ASPA;CPS1;CTH;DHPS;MAT1A;MTRR;PSAT1                                                                                                                                                                                                                                                                       |
| 21   | Eye Development                                                                   | GO     | 14        | 2.418142e-06 | 8.320  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | COL5A1;CRB2;CYP1B1;KRT12;PROM1;USH1G                                                                                                                                                                                                                                                                                    |
| 22   | REACTOME_DISEASES_OF_CARBOHYDRATE_METABOLISM                                      | MSIGDB | 34        | 2.417253e-06 | 11.170 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 23   | Proteolysis                                                                       | GO     | 178       | 2.416673e-06 | 7.488  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACE;ACE2;ADAM17;ADAM9;ADAMTS13;ADAMTS2;BMP1;CAPN3;CLPP;CTSK;F11;F7;FURIN;HTRA1;KLKB1;LONP1;LPA;MMP1;PAPPA;PIDD1;PLG;PROC;PSMB10;RELN;USP1                                                                                                                                                                               |
| 24   | Smell Perception                                                                  | GO     | 8         | 2.354748e-06 | 10.024 | 0.0212      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | CNGA3;GJB4                                                                                                                                                                                                                                                                                                              |
| 25   | WP_THYROID_HORMONES_PRODUCTION_AND_THEIR_PERIPHERAL_DOWNSTREAM_SIGNALLING_EFFECTS | MSIGDB | 92        | 2.305824e-06 | 10.181 | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 26   | REACTOME_SLC_TRANSPORTER_DISORDERS                                                | MSIGDB | 97        | 2.297809e-06 | 9.660  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 27   | HALLMARK_XENOBIOTIC_METABOLISM                                                    | MSIGDB | 195       | 2.274288e-06 | 8.069  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 28   | Organogenesis                                                                     | GO     | 105       | 2.208983e-06 | 6.547  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | COL11A1;COL17A1;COL1A2;COL4A3;COL4A6;COL5A1;COL7A1;DACH1;FBN1;FN1;GNRHR;HCCS;HNF1B;HNF4A;HSPG2;LAMA2;LAMA3;NPHS2;PAX3;TLE1;TLE3;TNC;VCAN                                                                                                                                                                                |
| 29   | Synaptic Transmission                                                             | GO     | 89        | 2.177547e-06 | 7.639  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | CACNA1G;CHRNE;CLCN1;CNTNAP2;DRD4;GABBR2;GRIA3;GRIN2B;GRIN2D;GRM6;HTR2C;KCNQ1;RAPSN;SLC1A4;TH;TRH                                                                                                                                                                                                                        |
| 30   | Metabolism                                                                        | GO     | 55        | 2.175265e-06 | 6.623  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | AUH;CYP17A1;CYP1B1;CYP21A2;CYP24A1;CYP27A1;CYP27B1;CYP2C19;CYP2C8;MMAB                                                                                                                                                                                                                                                  |
| 31   | WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR                    | MSIGDB | 93        | 2.154154e-06 | 8.301  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 32   | REACTOME_CILIUM_ASSEMBLY                                                          | MSIGDB | 194       | 2.153895e-06 | 6.715  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 33   | Histogenic Process                                                                | GO     | 35        | 2.153200e-06 | 8.923  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | COL17A1;KRT12;KRT16;KRT81;KRT9;LAMA3;TSHB                                                                                                                                                                                                                                                                               |
| 34   | REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION                                        | MSIGDB | 19        | 2.151692e-06 | 9.386  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 35   | recombinational repair                                                            | GO     | 24        | 2.135986e-06 | 7.630  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | MCM8;MRE11;NSMCE2;RAD50;RAD51C;RAD51D;RAD54B;RAD54L;RBBP8                                                                                                                                                                                                                                                               |
| 36   | Muscle Contraction                                                                | GO     | 44        | 2.132654e-06 | 7.625  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | ACTA1;CHRNE;CLCN1;DMPK;DYSF;MYBPC3;MYH2;NEB;PGAM2;RYR1;SLC6A8                                                                                                                                                                                                                                                           |
| 37   | REACTOME_ORGANELLE_BIOGENESIS_AND_MAINTENANCE                                     | MSIGDB | 264       | 2.109733e-06 | 6.086  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 38   | REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT                                           | MSIGDB | 8         | 2.086918e-06 | 9.622  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 39   | Complement Activation                                                             | GO     | 35        | 2.077666e-06 | 9.000  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | C3;C5;C6;C8B;C9;CFHR5                                                                                                                                                                                                                                                                                                   |
| 40   | REACTOME_VISUAL_PHOTOTRANSDUCTION                                                 | MSIGDB | 96        | 2.060800e-06 | 8.553  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 41   | Drug Metabolism                                                                   | GO     | 10        | 2.046679e-06 | 4.732  | 0.0212      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 | CYP2C19;CYP2C8                                                                                                                                                                                                                                                                                                          |
| 42   | WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK                                               | MSIGDB | 120       | 2.023458e-06 | 7.118  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 43   | REACTOME_SENSORY_PROCESSING_OF_SOUND_BY_OUTER_HAIR_CELLS_OF_THE_COCHLEA           | MSIGDB | 53        | 2.022892e-06 | 9.268  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 44   | REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                                        | MSIGDB | 296       | 1.988659e-06 | 4.943  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 45   | WP_COMPLEMENT_AND_COAGULATION_CASCADES                                            | MSIGDB | 57        | 1.977216e-06 | 8.357  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 46   | REACTOME_INTRINSIC_PATHWAY_OF_FIBRIN_CLOT_FORMATION                               | MSIGDB | 22        | 1.970223e-06 | 9.329  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 47   | PID_RHODOPSIN_PATHWAY                                                             | MSIGDB | 23        | 1.966937e-06 | 9.573  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 48   | REACTOME_TRNA_AMINOACYLATION                                                      | MSIGDB | 31        | 1.958255e-06 | 7.373  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 49   | REACTOME_DNA_REPAIR                                                               | MSIGDB | 280       | 1.948880e-06 | 4.727  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |
| 50   | REACTOME_THYROXINE_BIOSYNTHESIS                                                   | MSIGDB | 10        | 1.948367e-06 | 9.327  | 0.0143      | NA           | TP53;JAK1;MAPK1;TSHR;JAK2 |                                                                                                                                                                                                                                                                                                                         |

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

- "What biological processes are enriched in Neuroblastoma?"
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

### Seed genes (input to BIFO) — 1395 resolved

These genes were used as the starting point for signal propagation. They are typically genes carrying rare variants in the cohort, or genes of interest in the experimental system.

C0079419, C0079471, C0162832, C0242988, C0525037, C0598034, C0694888, C0809246, C0812228, C0812241, C0812246, C0812257, C0812265, C0812267, C0812278, C0812281, C0812295, C0812306, C0879290, C0879389, C0879391, C0879392, C0879393, C0919429, C0919481, C0919524, C0919528, C0919550, C1325331, C1332001, C1332013, C1332025, C1332029, C1332065, C1332072, C1332080, C1332120, C1332121, C1332381, C1332389, C1332409, C1332417, C1332427, C1332448, C1332656, C1332669, C1332721, C1332724, C1332766, C1332773, C1332774, C1332775, C1332828, C1333217, C1333219, C1333242, C1333254, C1333338, C1333348, C1333356, C1333357, C1333358, C1333359, C1333365, C1333523, C1333532, C1333535, C1333584, C1333677, C1333699, C1333710, C1333713, C1333892, C1333897, C1333921, C1334084, C1334123, C1334133, C1334134, C1334136, C1334137, C1334140, C1334290, C1334291, C1334345, C1334469, C1334867, C1334868, C1334870, C1334902, C1334912, C1335088, C1335092, C1335192, C1335232, C1335237, C1335265, C1335280, C1335289, C1335580, C1335607, C1335620, C1335647, C1335824, C1335843, C1335872, C1335873, C1336567, C1336639, C1336640, C1336642, C1336644, C1336658, C1336679, C1336688, C1336938, C1337001, C1337098, C1337106, C1364081, C1364508, C1366370, C1366529, C1366530, C1366536, C1366544, C1366631, C1366757, C1366767, C1366824, C1366882, C1367449, C1367472, C1367553, C1367578, C1367671, C1367701, C1367710, C1412061, C1412062, C1412064, C1412066, C1412070, C1412071, C1412081, C1412104, C1412107, C1412109, C1412110, C1412127, C1412132, C1412136, C1412166, C1412173, C1412179, C1412186, C1412192, C1412207, C1412208, C1412209, C1412225, C1412231, C1412284, C1412290, C1412291, C1412305, C1412335, C1412338, C1412340, C1412355, C1412364, C1412401, C1412404, C1412437, C1412459, C1412492, C1412521, C1412537, C1412553, C1412582, C1412588, C1412591, C1412610, C1412611, C1412625, C1412630, C1412631, C1412637, C1412641, C1412689, C1412691, C1412697, C1412718, C1412741, C1412749, C1412759, C1412777, C1412804, C1412833, C1412841
... and 1195 more (see --seeds file).

### Reference pathway set

These pathways were pre-specified as biologically relevant before the analysis was run. They are marked TRUE in the in_reference column.

No reference set provided.

### Analysis parameters

- Scores file: pathway_scores_standard.csv
- Total pathways scored: 2196
- Pathways with valid null distribution: 2196
- Null model: degree-preserving membership rewiring, N = 1,000 permutations
- Multiple testing correction: Benjamini-Hochberg
- Pathway size filter: minimum 8 members, maximum 300 members
- Pathway sources: MSigDB canonical collections, WikiPathways, Reactome, GO

---
*Generated by summarize_results.py on 2026-04-23*
