#!/usr/bin/env python3
"""
clean_cypher_output.py

Cleans cypher-shell --format plain output files.
Skips warning/error lines at the top, handles quoted fields with commas.

Usage:
    python clean_cypher_output.py kf_chd_nodes.csv kf_chd_nodes_clean.csv
    python clean_cypher_output.py kf_chd_edges_raw.csv kf_chd_edges_raw_clean.csv
    python clean_cypher_output.py kf_chd_pathway_membership_edges.csv kf_chd_pathway_membership_edges_clean.csv
"""
import csv, sys, re
from pathlib import Path

src, dst = sys.argv[1], sys.argv[2]

with open(src, encoding='utf-8-sig', errors='replace') as f:
    lines = f.readlines()

# Find first line that looks like a CSV header (starts with a quote or letter,
# contains commas, no spaces at start)
header_line = None
header_idx = 0
for i, line in enumerate(lines):
    stripped = line.strip()
    # Skip blank lines and Java warning lines
    if not stripped:
        continue
    if stripped.startswith('Failed') or stripped.startswith('java.') or stripped.startswith('WARNING'):
        continue
    # This looks like the header
    header_line = stripped
    header_idx = i
    break

if header_line is None:
    print(f"ERROR: Could not find header in {src}")
    sys.exit(1)

print(f"Header found at line {header_idx+1}: {header_line[:60]}")

# Parse header to get column count
import io
reader = csv.reader(io.StringIO(header_line))
header = next(reader)
n_cols = len(header)
print(f"Columns ({n_cols}): {header}")

# Now read all data lines from header_idx onwards
good_rows = []
skipped = 0
for line in lines[header_idx:]:
    stripped = line.strip()
    if not stripped:
        continue
    if stripped.startswith('Failed') or stripped.startswith('java.') or stripped.startswith('WARNING'):
        skipped += 1
        continue
    try:
        row = next(csv.reader(io.StringIO(stripped)))
        if len(row) == n_cols:
            good_rows.append(row)
        else:
            skipped += 1
    except Exception:
        skipped += 1

print(f"Good rows: {len(good_rows)-1} (excluding header), skipped: {skipped}")

# Strip spaces from header column names
if good_rows:
    good_rows[0] = [col.strip() for col in good_rows[0]]

# Strip surrounding double-quotes from all values (cypher-shell artifact)
def unquote(val):
    val = val.strip()
    if len(val) >= 2 and val.startswith('"') and val.endswith('"'):
        return val[1:-1]
    return val

cleaned_rows = []
for i, row in enumerate(good_rows):
    if i == 0:  # header — already stripped
        cleaned_rows.append(row)
    else:
        cleaned_rows.append([unquote(v) for v in row])

# Write clean CSV
with open(dst, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
    writer.writerows(cleaned_rows)

print(f"Written: {dst}")
