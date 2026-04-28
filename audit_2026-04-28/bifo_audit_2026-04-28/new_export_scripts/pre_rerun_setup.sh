#!/usr/bin/env bash
# =============================================================================
# pre_rerun_setup.sh
# Phase 0: Set up environment and capture pre-rerun state for the audit.
#
# Run this once before kicking off the actual rerun. It does:
#   1. Creates _rerun_2026-04-28/ working directory structure
#   2. Captures md5s of all current canonical inputs and outputs
#   3. Captures git status snapshot
#   4. Verifies conda environment and Python version
#   5. Verifies cypher-shell is available
#   6. Verifies pipeline_fixes/ patches have been applied
#
# Usage:
#   conda activate bifo-spectral
#   cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph
#   bash audit_2026-04-28/new_export_scripts/pre_rerun_setup.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
AUDIT_DIR="$REPO_DIR/audit_2026-04-28"
RERUN_DIR="$REPO_DIR/_rerun_2026-04-28"

echo "============================================================"
echo "BIFO-PPR Audit Pre-rerun Setup"
echo "  Repo:    $REPO_DIR"
echo "  Audit:   $AUDIT_DIR"
echo "  Rerun:   $RERUN_DIR"
echo "  Date:    $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "  Host:    $(hostname)"
echo "============================================================"

# =============================================================================
# 1. Create directory structure
# =============================================================================
echo ""
echo "[1] Creating directory structure ..."
mkdir -p "$RERUN_DIR"/{chd,nbl,tall,cbtn,rsbd,_logs}
mkdir -p "$AUDIT_DIR"/snapshots
echo "  Created: $RERUN_DIR/{chd,nbl,tall,cbtn,rsbd,_logs}"
echo "  Created: $AUDIT_DIR/snapshots"

# =============================================================================
# 2. Capture pre-rerun md5s
# =============================================================================
echo ""
echo "[2] Capturing pre-rerun md5s ..."

cd "$REPO_DIR"

SNAPSHOT="$AUDIT_DIR/snapshots/pre_rerun_md5s_$(date +%Y%m%d_%H%M%S).txt"
{
    echo "Pre-rerun md5 snapshot"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Host: $(hostname)"
    echo "Git HEAD: $(git rev-parse HEAD)"
    echo ""
    echo "--- Seed files ---"
    for cohort in chd nbl tall cbtn rsbd; do
        for seedfile in "data/cohorts/${cohort}"/*.txt; do
            if [ -f "$seedfile" ]; then
                md5sum "$seedfile" 2>/dev/null
            fi
        done
    done
    echo ""
    echo "--- Pipeline scripts ---"
    md5sum pipeline/bifo_conditioning.py \
           pipeline/score_pathways.py \
           pipeline/summarize_results.py \
           pipeline/seed_cui_lookup.py \
           pipeline/clean_cypher_output.py \
           pipeline/generate_export_cypher.py 2>/dev/null
    echo ""
    echo "--- Existing edges/nodes (canonical) ---"
    md5sum results/kf_chd/edges_all_noncc.csv.gz \
           results/kf_chd/nodes_clean_noncc.csv.gz \
           results/kf_nbl/edges_all_noncc.csv.gz \
           results/kf_nbl/nodes_clean_noncc.csv.gz 2>/dev/null
    echo ""
    echo "--- Cypher query files ---"
    md5sum cypher/*.cypher 2>/dev/null
    echo ""
    echo "--- Config ---"
    md5sum config/*.yaml 2>/dev/null || echo "  (no yamls in config/)"
} > "$SNAPSHOT"
echo "  Wrote: $SNAPSHOT"

# =============================================================================
# 3. Capture git status
# =============================================================================
echo ""
echo "[3] Capturing git status ..."
GIT_STATUS="$AUDIT_DIR/snapshots/git_status_$(date +%Y%m%d_%H%M%S).txt"
{
    echo "Git status snapshot"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo ""
    echo "--- HEAD ---"
    git log -1 --oneline
    echo ""
    echo "--- Branch ---"
    git branch -v
    echo ""
    echo "--- Status ---"
    git status
    echo ""
    echo "--- Last 20 commits ---"
    git log --oneline -20
} > "$GIT_STATUS"
echo "  Wrote: $GIT_STATUS"

# =============================================================================
# 4. Verify conda environment and Python version
# =============================================================================
echo ""
echo "[4] Verifying environment ..."

if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "  Conda env: $CONDA_DEFAULT_ENV"
    if [ "$CONDA_DEFAULT_ENV" != "bifo-spectral" ]; then
        echo "  WARNING: Expected 'bifo-spectral'."
        echo "  Run: conda activate bifo-spectral"
    fi
else
    echo "  WARNING: No conda environment activated."
    echo "  Run: conda activate bifo-spectral"
fi

PYTHON_VERSION=$(python3 --version 2>&1)
echo "  Python: $PYTHON_VERSION"
if [[ "$PYTHON_VERSION" != *"Python 3.11"* ]]; then
    echo "  WARNING: Expected Python 3.11.x"
fi

# =============================================================================
# 5. Verify cypher-shell
# =============================================================================
echo ""
echo "[5] Verifying cypher-shell ..."
if command -v cypher-shell &>/dev/null; then
    CYPHER_VERSION=$(cypher-shell --version 2>&1 | head -1 || echo "version unknown")
    echo "  Found: $(which cypher-shell)"
    echo "  Version: $CYPHER_VERSION"
else
    echo "  WARNING: cypher-shell not on PATH."
    echo "  Either install cypher-shell or set CYPHER_SHELL env var."
fi

# =============================================================================
# 6. Verify pipeline_fixes/ have been applied
# =============================================================================
echo ""
echo "[6] Verifying pipeline_fixes/ applied ..."

# Check that generate_export_cypher.py has Query 6
if grep -q "make_query6" "$REPO_DIR/pipeline/generate_export_cypher.py" 2>/dev/null; then
    echo "  generate_export_cypher.py: Query 6 PRESENT"
else
    echo "  generate_export_cypher.py: Query 6 MISSING"
    echo "  ACTION: Apply pipeline_fixes/generate_export_cypher.py:"
    echo "    cp $AUDIT_DIR/pipeline_fixes/generate_export_cypher.py $REPO_DIR/pipeline/"
fi

# =============================================================================
# 7. Verify seed files for all five cohorts
# =============================================================================
echo ""
echo "[7] Verifying seed files for all cohorts ..."
all_seeds_ok=1
for cohort in chd nbl tall cbtn rsbd; do
    seedfile="$REPO_DIR/data/cohorts/${cohort}/kf_${cohort}_seeds_maf001.txt"
    if [ -f "$seedfile" ]; then
        n_genes=$(grep -vc '^#' "$seedfile" 2>/dev/null || echo 0)
        echo "  $cohort: OK ($n_genes seeds, $seedfile)"
    else
        echo "  $cohort: MISSING ($seedfile)"
        all_seeds_ok=0
    fi
done
if [ $all_seeds_ok -eq 0 ]; then
    echo ""
    echo "  ACTION: Generate missing seed files from Ben's MAF gene tables"
    echo "  Format: one HGNC symbol per line, optionally with carrier count (tab-separated)"
fi

# =============================================================================
# Done
# =============================================================================
echo ""
echo "============================================================"
echo "Pre-rerun setup complete."
echo ""
echo "Next steps (in order):"
echo "  1. Verify warnings above are addressed."
echo "  2. Phase 1: bash $SCRIPT_DIR/run_all_cohorts.sh --cohorts \"chd nbl\""
echo "  3. After Phase 1 verification: bash $SCRIPT_DIR/run_all_cohorts.sh --cohorts \"tall cbtn rsbd\""
echo "============================================================"
