# Run manifests

Compact byte-level descriptions of known-good pipeline run directories. Used
by regression tests to detect whether code changes silently alter pipeline
output.

## What's a manifest?

JSONL — one JSON object per line. The first line is a **header** describing
the capture context (when, where, git SHA at capture time, host metadata).
Subsequent lines are **file records**: one per regular file in the run dir,
sorted by relative path.

Per-file records carry:

- `path` — relative to the run dir, forward slashes
- `size` — file size in bytes
- `sha256` — SHA-256 of file contents

The manifest does **not** include mtimes, absolute paths, or anything else
that would vary across machines or transfers. Two captures of the same run
dir on different hosts (assuming identical contents) produce byte-identical
manifests.

## Capturing a manifest

```bash
python tests/regression/run_manifests/capture_run_manifest.py \
    --run-dir /path/to/pipeline/run/output \
    --manifest tests/regression/run_manifests/<cohort>_<date>.run_manifest.jsonl \
    --repo-root /path/to/bifo-graph
```

`--repo-root` is the bifo-graph git working tree at capture time. The git
SHA, status (clean/dirty), and `git describe` output are recorded in the
manifest header. This matters for provenance: a manifest carries the SHA
of the code that produced the run it describes.

If `--repo-root` is omitted, all git fields in the header are `null` and
the manifest still works for diff purposes.

## Comparing manifests

Manifests are sorted JSONL, so plain text diff works:

```bash
diff -u old.run_manifest.jsonl new.run_manifest.jsonl
```

Or with git's color/pager support:

```bash
git diff --no-index old.run_manifest.jsonl new.run_manifest.jsonl
```

If only the header differs, the file contents are identical between runs.
If file records differ, individual file lines will show in the diff with
their `path` and `sha256` changes visible.

## When to capture a new manifest

A manifest is the byte-level "lock" for a known-good pipeline state. Update
or add manifests when:

- A new authoritative run lands on the HPC (e.g., a new MAF cut, new cohort,
  or a substrate fix that's been validated independently). Capture the
  manifest, commit alongside any code changes that produced it.
- A code change is *intentionally* expected to alter outputs (e.g., the
  `_compute_global_gene_stats` filter change). Validate the new output is
  correct by other means, then update the affected manifest.

Don't update a manifest just to make a regression test pass. If the test is
flagging a change you didn't expect, that's the bug.

## File naming convention

`<cohort>_<date>.run_manifest.jsonl` — e.g., `chd_apr28.run_manifest.jsonl`.
Date can be more specific (`chd_2026-04-28`) if multiple runs in a month.

One manifest per cohort, even if cohorts share substrate. Cross-cohort
manifests are rejected as a design choice — separate cohorts always get
separate manifests, simpler diff semantics.

## Schema versioning

The header has a `schema_version` field (currently `1`). If the schema
changes (e.g., adding a `q5_provenance` field for user-supplied null
backgrounds), bump this and update the loader. Old manifests remain readable
under their original schema.

## Future extensions

Anticipated:

- A `q5_provenance` field for runs using a user-supplied precomputed Q5
  null background, capturing where the precomputed background came from.
- Optional `excluded_paths` glob patterns in the capture command, for skipping
  ephemera (logs, debug dumps) that differ between runs but aren't outputs.

Both deferred until there's a concrete consumer requesting them.
