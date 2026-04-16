#!/usr/bin/env bash
# =============================================================================
# merge_files.sh
# Stage 2.3: Merge cleaned CSV files using pandas (header-safe).
#
# Usage:
#   bash merge_files.sh [COHORT]
#   COHORT defaults to "chd"
#
# Inputs:
#   kf_{cohort}_nodes_clean.csv
#   kf_{cohort}_ncc_pathway_nodes.csv
#   kf_{cohort}_edges_raw_clean.csv
#   kf_{cohort}_pathway_membership_edges_clean.csv
#   kf_{cohort}_ncc_membership_edges.csv
#
# Outputs:
#   kf_{cohort}_nodes_extended.csv   (nodes_clean + ncc_pathway_nodes)
#   kf_{cohort}_edges_merged.csv     (edges_raw_clean + pathway_membership)
#   kf_{cohort}_edges_all.csv        (edges_merged + ncc_membership_edges)
#
# NOTE: Uses pandas pd.concat — never cat. The cat command does not skip
# duplicate CSV headers and produces malformed files.
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"

echo "============================================================"
echo "Merging CSV files for cohort: $COHORT"
echo "============================================================"

python3 - <<PYEOF
import pandas as pd
import sys

cohort = "${COHORT}"

def merge(files, out):
    dfs = [pd.read_csv(f, dtype=str) for f in files]
    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(out, index=False)
    print(f"  {out}: {len(merged):,} rows, columns={list(merged.columns)}")

merge(
    [f"kf_{cohort}_nodes_clean.csv",
     f"kf_{cohort}_ncc_pathway_nodes.csv"],
    f"kf_{cohort}_nodes_extended.csv"
)

merge(
    [f"kf_{cohort}_edges_raw_clean.csv",
     f"kf_{cohort}_pathway_membership_edges_clean.csv"],
    f"kf_{cohort}_edges_merged.csv"
)

merge(
    [f"kf_{cohort}_edges_merged.csv",
     f"kf_{cohort}_ncc_membership_edges.csv"],
    f"kf_{cohort}_edges_all.csv"
)

print("")
print("Merge complete.")
PYEOF

echo "Next step: bash run_seed_lookup.sh $COHORT"
