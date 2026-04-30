"""capture_run_manifest.py — produce a JSONL manifest of an HPC run directory.

Used by tests/regression/ to lock byte-level expected outputs of a known-good
pipeline run. The manifest is committed to git and serves as the diff target
for subsequent regression checks.

The manifest is JSONL (one JSON object per line):

  - First line: a header record describing the capture context (when, where,
    git SHA of the bifo-graph clone at capture time, etc.).
  - Subsequent lines: one record per regular file in the run dir, sorted by
    relative path.

To compare manifests later, use `git diff --no-index <old> <new>` or `diff`.

Usage:
    python tests/regression/run_manifests/capture_run_manifest.py \\
        --run-dir /path/to/pipeline/run \\
        --manifest /path/to/output.run_manifest.jsonl \\
        [--repo-root /path/to/bifo-graph]

The --repo-root argument controls which git working tree to introspect for
the manifest header. If omitted, git fields in the header are null.

This script has no third-party dependencies. It uses only the standard
library and shells out to `git` (optional).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Optional

TOOL_VERSION = "0.1.0"
SCHEMA_VERSION = 1
HASH_BUFFER_BYTES = 1024 * 1024  # 1 MiB chunks for incremental hashing


def hash_file_sha256(path: Path) -> str:
    """Compute SHA-256 of file contents in fixed-size chunks.

    Memory-bounded: never reads more than HASH_BUFFER_BYTES at once,
    so a 10 GB file hashes with the same memory footprint as a 1 KB file.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(HASH_BUFFER_BYTES)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _run_git(args: list[str], cwd: Path) -> Optional[str]:
    """Run a git subcommand in `cwd`, return stdout stripped of trailing
    whitespace. Returns None on any failure (non-zero exit, missing git,
    not a git repo). Errors warned to stderr.
    """
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=10,
        )
    except FileNotFoundError:
        print("warning: git not found in PATH; git fields will be null",
              file=sys.stderr)
        return None
    except subprocess.TimeoutExpired:
        print(f"warning: git {' '.join(args)} timed out; field will be null",
              file=sys.stderr)
        return None

    if result.returncode != 0:
        # Don't warn on every non-zero exit — `git describe` fails normally
        # on repos with no tags. Caller decides whether to warn.
        return None
    return result.stdout.strip()


def get_git_context(repo_root: Optional[Path]) -> dict:
    """Capture git SHA, working-tree cleanliness, and `git describe` for
    repo_root. All values are None if repo_root is None, not a git working
    tree, or git unavailable.

    Returns dict with keys: git_sha, git_status_clean, git_describe.
    """
    null_context = {
        "git_sha": None,
        "git_status_clean": None,
        "git_describe": None,
    }

    if repo_root is None:
        return null_context

    if not repo_root.is_dir():
        print(f"warning: --repo-root {repo_root} is not a directory; "
              "git fields will be null", file=sys.stderr)
        return null_context

    # Verify it's actually a git working tree before issuing further commands.
    inside = _run_git(["rev-parse", "--is-inside-work-tree"], repo_root)
    if inside != "true":
        print(f"warning: --repo-root {repo_root} is not a git working tree; "
              "git fields will be null", file=sys.stderr)
        return null_context

    sha = _run_git(["rev-parse", "HEAD"], repo_root)
    porcelain = _run_git(["status", "--porcelain"], repo_root)
    describe = _run_git(["describe", "--always", "--dirty"], repo_root)

    # status --porcelain returns empty string for a clean tree; non-empty
    # output means there are modified/untracked/staged files.
    status_clean = (porcelain == "") if porcelain is not None else None
    if status_clean is False:
        print(f"warning: working tree at {repo_root} is dirty; "
              "manifest header reflects uncommitted state", file=sys.stderr)

    return {
        "git_sha": sha,
        "git_status_clean": status_clean,
        "git_describe": describe,
    }


def build_header_record(run_dir: Path, repo_root: Optional[Path]) -> dict:
    """Construct the JSONL header record describing capture context."""
    git_ctx = get_git_context(repo_root)
    return {
        "record_type": "header",
        "tool": "capture_run_manifest.py",
        "tool_version": TOOL_VERSION,
        "schema_version": SCHEMA_VERSION,
        "captured_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hostname": socket.gethostname(),
        "python_version": sys.version.split()[0],
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV"),
        "git_sha": git_ctx["git_sha"],
        "git_status_clean": git_ctx["git_status_clean"],
        "git_describe": git_ctx["git_describe"],
        "run_dir_abs_path": str(run_dir.resolve()),
        "run_dir_basename": run_dir.name,
    }


def iter_file_records(run_dir: Path) -> Iterator[dict]:
    """Walk run_dir and yield one record per regular file (sorted by path).

    Sorting ensures two captures of the same directory produce byte-identical
    manifests. Symlinks are not followed and are not included as records.
    Empty directories are not represented (only files).
    """
    file_paths: list[Path] = []
    for root, dirs, files in os.walk(run_dir, followlinks=False):
        # Sort dirs in-place so traversal order is deterministic;
        # final sort below provides the actual ordering guarantee, but
        # this keeps memory peak lower for very large trees.
        dirs.sort()
        for name in files:
            p = Path(root) / name
            # Skip symlinks (don't follow, don't record).
            if p.is_symlink():
                continue
            # Skip anything that isn't a regular file (sockets, fifos).
            if not p.is_file():
                continue
            file_paths.append(p)

    file_paths.sort()

    for p in file_paths:
        rel = p.relative_to(run_dir)
        yield {
            "record_type": "file",
            "path": str(rel),
            "size": p.stat().st_size,
            "sha256": hash_file_sha256(p),
        }


def write_manifest(
    manifest_path: Path,
    header: dict,
    file_records: Iterator[dict],
) -> int:
    """Write JSONL manifest. Header first, then file records, one per line.

    Creates parent directories as needed. Returns the number of file records
    written (excludes header).
    """
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    n_files = 0
    with open(manifest_path, "w") as f:
        f.write(json.dumps(header, sort_keys=True) + "\n")
        for rec in file_records:
            f.write(json.dumps(rec, sort_keys=True) + "\n")
            n_files += 1
    return n_files


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--run-dir", type=Path, required=True,
        help="Directory whose contents will be hashed and recorded.",
    )
    parser.add_argument(
        "--manifest", type=Path, required=True,
        help="Output manifest path (JSONL). Parent dirs created if missing.",
    )
    parser.add_argument(
        "--repo-root", type=Path, default=None,
        help="Path to the bifo-graph git working tree. Used to capture "
             "git SHA, status, and describe in the manifest header. "
             "If omitted, git fields are null.",
    )
    args = parser.parse_args()

    if not args.run_dir.is_dir():
        print(f"error: --run-dir {args.run_dir} is not a directory",
              file=sys.stderr)
        return 2

    header = build_header_record(args.run_dir, args.repo_root)
    n = write_manifest(args.manifest, header, iter_file_records(args.run_dir))

    print(f"wrote {n} file records to {args.manifest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
