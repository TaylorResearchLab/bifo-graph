"""Unit tests for tests/regression/run_manifests/capture_run_manifest.py.

The capture script lives outside `pipeline/`, so we load it as a module via
importlib rather than a regular import.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# ----------------------------------------------------------------------
# Module loading
# ----------------------------------------------------------------------

_THIS_DIR = Path(__file__).resolve().parent
_CAPTURE_SCRIPT = (
    _THIS_DIR.parent / "regression" / "run_manifests" / "capture_run_manifest.py"
)

_spec = importlib.util.spec_from_file_location(
    "capture_run_manifest", _CAPTURE_SCRIPT
)
capture_run_manifest = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(capture_run_manifest)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _git_available() -> bool:
    """Whether `git` is on the PATH."""
    return shutil.which("git") is not None


def _make_dir_with_files(root: Path, files: dict[str, bytes]) -> None:
    """Create files at given relative paths under `root`. Parents created."""
    for rel, content in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(content)


def _init_git_repo(repo_root: Path, dirty: bool = False) -> str:
    """Initialise a git repo at `repo_root` with one commit. Returns the
    commit's SHA. If dirty=True, leaves an uncommitted change so
    `git status --porcelain` reports non-empty output.
    """
    repo_root.mkdir(parents=True, exist_ok=True)
    env = {
        **os.environ,
        # Pin commit identity so SHAs are reproducible across test runs.
        "GIT_AUTHOR_NAME": "test",
        "GIT_AUTHOR_EMAIL": "test@example.com",
        "GIT_COMMITTER_NAME": "test",
        "GIT_COMMITTER_EMAIL": "test@example.com",
    }
    run = lambda *args: subprocess.run(
        ["git", *args], cwd=repo_root, check=True,
        capture_output=True, text=True, env=env,
    )
    run("init", "--initial-branch=main")
    (repo_root / "README.md").write_text("test\n")
    run("add", "README.md")
    run("commit", "-m", "initial")
    sha = run("rev-parse", "HEAD").stdout.strip()

    if dirty:
        # Modify a TRACKED file. Adding an untracked file would make
        # `git status --porcelain` non-empty (status_clean=False) but
        # `git describe --dirty` would NOT flag it; describe only checks
        # changes to tracked files. We want both to flag dirty here.
        (repo_root / "README.md").write_text("modified\n")

    return sha


# ----------------------------------------------------------------------
# Tests: hash_file_sha256
# ----------------------------------------------------------------------

class TestHashFileSha256:
    def test_known_content(self, tmp_path: Path):
        """SHA-256 of a known string matches the standard digest."""
        p = tmp_path / "a.txt"
        content = b"hello world\n"
        p.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert capture_run_manifest.hash_file_sha256(p) == expected

    def test_empty_file(self, tmp_path: Path):
        """SHA-256 of empty file is the well-known empty digest."""
        p = tmp_path / "empty.txt"
        p.write_bytes(b"")
        # SHA-256 of empty input
        empty_digest = (
            "e3b0c44298fc1c149afbf4c8996fb924"
            "27ae41e4649b934ca495991b7852b855"
        )
        assert capture_run_manifest.hash_file_sha256(p) == empty_digest

    def test_large_file_streamed(self, tmp_path: Path):
        """Hashing a file larger than HASH_BUFFER_BYTES returns same digest
        as hashing the contents in one shot. Confirms chunked hashing is
        equivalent to monolithic."""
        # Create file slightly larger than the buffer to exercise multi-chunk path.
        size = capture_run_manifest.HASH_BUFFER_BYTES + 7
        content = bytes(range(256)) * (size // 256 + 1)
        content = content[:size]
        p = tmp_path / "big.bin"
        p.write_bytes(content)
        expected = hashlib.sha256(content).hexdigest()
        assert capture_run_manifest.hash_file_sha256(p) == expected


# ----------------------------------------------------------------------
# Tests: get_git_context
# ----------------------------------------------------------------------

class TestGetGitContext:
    @pytest.mark.skipif(not _git_available(), reason="git not installed")
    def test_clean_repo(self, tmp_path: Path):
        """Clean working tree: git_status_clean is True, sha matches."""
        sha = _init_git_repo(tmp_path, dirty=False)
        ctx = capture_run_manifest.get_git_context(tmp_path)
        assert ctx["git_sha"] == sha
        assert ctx["git_status_clean"] is True
        assert ctx["git_describe"] is not None

    @pytest.mark.skipif(not _git_available(), reason="git not installed")
    def test_dirty_repo(self, tmp_path: Path):
        """Dirty working tree: git_status_clean is False."""
        sha = _init_git_repo(tmp_path, dirty=True)
        ctx = capture_run_manifest.get_git_context(tmp_path)
        assert ctx["git_sha"] == sha
        assert ctx["git_status_clean"] is False
        # describe with --dirty appends "-dirty" suffix
        assert ctx["git_describe"] is not None
        assert "dirty" in ctx["git_describe"]

    def test_non_git_dir(self, tmp_path: Path):
        """Plain (non-git) directory: all git fields None."""
        ctx = capture_run_manifest.get_git_context(tmp_path)
        assert ctx == {
            "git_sha": None,
            "git_status_clean": None,
            "git_describe": None,
        }

    def test_none_repo_root(self):
        """repo_root=None: all git fields None."""
        ctx = capture_run_manifest.get_git_context(None)
        assert ctx == {
            "git_sha": None,
            "git_status_clean": None,
            "git_describe": None,
        }

    def test_nonexistent_repo_root(self, tmp_path: Path):
        """Path that doesn't exist: all git fields None, no crash."""
        bogus = tmp_path / "does_not_exist"
        ctx = capture_run_manifest.get_git_context(bogus)
        assert ctx["git_sha"] is None


# ----------------------------------------------------------------------
# Tests: iter_file_records
# ----------------------------------------------------------------------

class TestIterFileRecords:
    def test_records_sorted(self, tmp_path: Path):
        """Records yielded in path-sorted order regardless of creation order."""
        _make_dir_with_files(tmp_path, {
            "z.txt": b"z",
            "a.txt": b"a",
            "m.txt": b"m",
        })
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        paths = [r["path"] for r in records]
        assert paths == ["a.txt", "m.txt", "z.txt"]

    def test_record_fields(self, tmp_path: Path):
        """Each record has exactly the expected fields and types."""
        _make_dir_with_files(tmp_path, {"f.txt": b"hello"})
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        assert len(records) == 1
        r = records[0]
        assert set(r.keys()) == {"record_type", "path", "size", "sha256"}
        assert r["record_type"] == "file"
        assert r["path"] == "f.txt"
        assert r["size"] == 5
        assert r["sha256"] == hashlib.sha256(b"hello").hexdigest()

    def test_recursive(self, tmp_path: Path):
        """Files in subdirectories use forward-slash relative paths."""
        _make_dir_with_files(tmp_path, {
            "top.txt": b"t",
            "sub/nested.txt": b"n",
            "sub/deeper/leaf.txt": b"l",
        })
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        paths = sorted(r["path"] for r in records)
        # Note: os.path.join uses os.sep, but Path.relative_to + str gives
        # the platform separator. We test that the paths are present in
        # whichever form the OS uses.
        # On POSIX (which is where the pipeline runs), this is forward slashes.
        assert "top.txt" in paths
        assert any("sub" in p and "nested.txt" in p for p in paths)
        assert any("deeper" in p and "leaf.txt" in p for p in paths)

    def test_skips_directories(self, tmp_path: Path):
        """Empty directories produce no records."""
        (tmp_path / "empty_subdir").mkdir()
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        assert records == []

    def test_does_not_follow_symlinks(self, tmp_path: Path):
        """Symlinks not yielded (don't follow, don't record)."""
        target = tmp_path / "real.txt"
        target.write_bytes(b"real")
        link = tmp_path / "link.txt"
        link.symlink_to(target)
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        paths = [r["path"] for r in records]
        # "real.txt" should be present, "link.txt" should not.
        assert "real.txt" in paths
        assert "link.txt" not in paths

    def test_empty_dir(self, tmp_path: Path):
        """Empty run dir: no records."""
        records = list(capture_run_manifest.iter_file_records(tmp_path))
        assert records == []


# ----------------------------------------------------------------------
# Tests: end-to-end (main flow)
# ----------------------------------------------------------------------

class TestEndToEnd:
    def test_manifest_is_valid_jsonl(self, tmp_path: Path):
        """Output is parseable JSONL with header followed by file records."""
        run_dir = tmp_path / "run"
        _make_dir_with_files(run_dir, {
            "a.txt": b"hello",
            "sub/b.txt": b"world",
        })
        manifest_path = tmp_path / "out" / "test.run_manifest.jsonl"
        header = capture_run_manifest.build_header_record(run_dir, None)
        n = capture_run_manifest.write_manifest(
            manifest_path,
            header,
            capture_run_manifest.iter_file_records(run_dir),
        )
        assert n == 2

        with open(manifest_path) as f:
            lines = [json.loads(line) for line in f]
        assert len(lines) == 3  # header + 2 files
        assert lines[0]["record_type"] == "header"
        assert lines[1]["record_type"] == "file"
        assert lines[2]["record_type"] == "file"
        # File records sorted by path
        assert lines[1]["path"] < lines[2]["path"]

    def test_creates_parent_dirs(self, tmp_path: Path):
        """Manifest path with non-existent parent dirs: parents created."""
        run_dir = tmp_path / "run"
        _make_dir_with_files(run_dir, {"f.txt": b"x"})
        manifest_path = tmp_path / "deep" / "nested" / "path" / "m.jsonl"
        assert not manifest_path.parent.exists()
        header = capture_run_manifest.build_header_record(run_dir, None)
        capture_run_manifest.write_manifest(
            manifest_path,
            header,
            capture_run_manifest.iter_file_records(run_dir),
        )
        assert manifest_path.exists()
        assert manifest_path.parent.is_dir()

    def test_empty_run_dir_produces_header_only(self, tmp_path: Path):
        """Empty run dir: manifest has just one line (the header)."""
        run_dir = tmp_path / "empty_run"
        run_dir.mkdir()
        manifest_path = tmp_path / "m.jsonl"
        header = capture_run_manifest.build_header_record(run_dir, None)
        n = capture_run_manifest.write_manifest(
            manifest_path,
            header,
            capture_run_manifest.iter_file_records(run_dir),
        )
        assert n == 0

        with open(manifest_path) as f:
            lines = f.readlines()
        assert len(lines) == 1
        assert json.loads(lines[0])["record_type"] == "header"

    def test_two_captures_of_same_dir_byte_identical(self, tmp_path: Path):
        """Capturing the same directory twice (without git) produces
        byte-identical manifests aside from the captured_at timestamp.

        We strip the timestamp before comparing because that's the only
        legitimate source of difference between two captures of the same
        contents."""
        run_dir = tmp_path / "run"
        _make_dir_with_files(run_dir, {
            "a.txt": b"alpha",
            "b/c.txt": b"gamma",
            "b/d.txt": b"delta",
        })

        def capture_to(path):
            header = capture_run_manifest.build_header_record(run_dir, None)
            capture_run_manifest.write_manifest(
                path,
                header,
                capture_run_manifest.iter_file_records(run_dir),
            )
            return path.read_text().splitlines()

        m1 = capture_to(tmp_path / "m1.jsonl")
        m2 = capture_to(tmp_path / "m2.jsonl")

        # Strip captured_at field from headers before comparing
        h1 = json.loads(m1[0])
        h2 = json.loads(m2[0])
        h1.pop("captured_at")
        h2.pop("captured_at")
        assert h1 == h2

        # File records (everything after the header line) must be identical
        assert m1[1:] == m2[1:]

    @pytest.mark.skipif(not _git_available(), reason="git not installed")
    def test_header_includes_git_when_repo_root_supplied(self, tmp_path: Path):
        """With --repo-root pointing at a git repo, header carries SHA."""
        repo_root = tmp_path / "repo"
        sha = _init_git_repo(repo_root, dirty=False)

        run_dir = tmp_path / "run"
        _make_dir_with_files(run_dir, {"f.txt": b"x"})

        header = capture_run_manifest.build_header_record(run_dir, repo_root)
        assert header["git_sha"] == sha
        assert header["git_status_clean"] is True


# ----------------------------------------------------------------------
# Tests: CLI behavior
# ----------------------------------------------------------------------

class TestCLI:
    def test_cli_invocation(self, tmp_path: Path, monkeypatch):
        """Running main() with valid args produces a manifest."""
        run_dir = tmp_path / "run"
        _make_dir_with_files(run_dir, {"hello.txt": b"hi"})
        manifest = tmp_path / "out.jsonl"

        monkeypatch.setattr(sys, "argv", [
            "capture_run_manifest.py",
            "--run-dir", str(run_dir),
            "--manifest", str(manifest),
        ])
        rc = capture_run_manifest.main()
        assert rc == 0
        assert manifest.exists()

        with open(manifest) as f:
            lines = [json.loads(line) for line in f]
        assert len(lines) == 2  # header + 1 file
        assert lines[1]["path"] == "hello.txt"

    def test_cli_rejects_nonexistent_run_dir(
        self, tmp_path: Path, monkeypatch, capsys
    ):
        """--run-dir must point at an existing directory."""
        bogus = tmp_path / "does_not_exist"
        manifest = tmp_path / "out.jsonl"

        monkeypatch.setattr(sys, "argv", [
            "capture_run_manifest.py",
            "--run-dir", str(bogus),
            "--manifest", str(manifest),
        ])
        rc = capture_run_manifest.main()
        assert rc == 2
        captured = capsys.readouterr()
        assert "is not a directory" in captured.err
