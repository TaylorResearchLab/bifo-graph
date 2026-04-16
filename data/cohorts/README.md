# KF Cohort Seed Gene Files

## File naming convention

| File | MAF | Carrier filter | N genes | Purpose |
|------|-----|---------------|---------|---------|
| kf_chd_seeds.txt | 0.0001 | n_carriers>=3 | 56 | Original analysis (ultra-rare) |
| kf_chd_seeds_maf001.txt | 0.001 | none | 1287 | Poster-matched (recommended) |
| kf_chd_seeds_maf01.txt | 0.01 | none | 1405 | Broadest filter |
| kf_nbl_seeds.txt | 0.0001 | n_carriers>=3 | 88 | Original analysis (ultra-rare) |
| kf_nbl_seeds_maf001.txt | 0.001 | none | 1406 | Poster-matched (recommended) |
| kf_nbl_seeds_maf01.txt | 0.01 | none | 1507 | Broadest filter |

## Recommended for BIFO analysis
Use **maf001** files — these match the MAF threshold used in the U24 poster
(Stear et al. CFDE Meeting 2026) which found stable cilia enrichment signal.

## Source
All files generated from AutoGVP P/LP variants from Kids First WGS cohorts:
- KF-CHD: phs001138 (PCGC, n=697)
- KF-NBL: phs001436 (n=460)

GATK PASS, GQ>=20, DP>=10, gnomAD MAF filter as specified in filename.
Background burden genes excluded: ABCA4, USH2A, G6PD, TTN, FLG, OBSCN,
MYO15A, MYO3A, MYO7A, MYO1A, TCHH, KRT71, KRT86, PADI3, COL4A5, TBL1Y,
GJB2, BCHE, PAH, CD36.
