# Strict-filter exemplar cypher

These are early exploratory cypher queries from a stricter
filtering regime than the production manuscript analysis:

| File | Seeds | Filter |
|------|------:|--------|
| kf_chd_export_queries.cypher | 56 | MAF≤0.0001 gnomAD, n_carriers≥3 |
| kf_nbl_export_queries.cypher | 88 | MAF≤0.0001 gnomAD, n_carriers≥3 |

The production analyses reported in the manuscript use a looser filter
(MAF≤0.001, n_carriers≥1) yielding 1,276 CHD and 1,395 NBL seeds.
Production cypher is at `cypher/kf_chd_export_queries.cypher` and
`cypher/kf_nbl_export_queries.cypher` (regenerated at runtime from
`data/cohorts/{chd,nbl}/kf_*_seeds_maf001.txt` by
`scripts/run_kf_{chd,nbl}_export.sh`).

These strict-filter exemplars are retained for reference and as small,
fast-running test fixtures.
