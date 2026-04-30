"""
pipeline/check_configs.py

Pre-flight CLI tool that validates the bifo-graph configuration before
a pipeline run. Validates pipeline_config.yaml schema, bifo_mapping
structure, and cross-consistency between the two files. Reports
problems in human-readable form.

Designed for two audiences:
  - Pipeline users (including reviewers) running configuration spot-checks
    before committing to a long pipeline run.
  - CI / automation, via numeric exit codes.

Usage:
    python -m pipeline.check_configs                      # default paths
    python -m pipeline.check_configs --pipeline-config X
    python -m pipeline.check_configs --bifo-mapping Y
    python -m pipeline.check_configs --quiet              # only print on failure

Exit codes:
    0  All checks passed.
    1  Configuration errors found (cross-validation, schema, structure).
    2  Misuse or setup failure (file not found, can't read, bad CLI args).
"""
import argparse
import difflib
import importlib.util
import sys
from pathlib import Path

# Load the pipeline_config module from a sibling file. This lets the CLI
# work without bifo-graph being installed as a package.
_HERE = Path(__file__).resolve().parent
_PIPELINE_CONFIG = _HERE / "pipeline_config.py"
_spec = importlib.util.spec_from_file_location("pipeline_config", _PIPELINE_CONFIG)
pipeline_config = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pipeline_config)

DEFAULT_PIPELINE_CONFIG = "config/pipeline_config.yaml"
DEFAULT_BIFO_MAPPING = "config/bifo_mapping_ddkg.yaml"

EXIT_OK = 0
EXIT_CONFIG_ERROR = 1
EXIT_SETUP_ERROR = 2

# Markers used in human-readable output. Plain ASCII only so output is
# safe across terminal/log/CI environments.
MARK_OK = "[OK]"
MARK_FAIL = "[FAIL]"


def _print_header(
    pipeline_config_path: Path,
    bifo_mapping_path: Path,
    out,
) -> None:
    print("Validating bifo-graph configuration...", file=out)
    print(file=out)
    print(f"  Pipeline config: {pipeline_config_path}", file=out)
    print(f"  KG mapping:      {bifo_mapping_path}", file=out)
    print(file=out)


def _suggest_predicate_alternative(
    typo: str,
    known_predicates: list[str],
) -> str | None:
    """Find the closest match in known_predicates for a likely typo.

    Returns the suggestion or None if nothing close enough exists.
    Uses difflib.get_close_matches (stdlib).
    """
    matches = difflib.get_close_matches(typo, known_predicates, n=1, cutoff=0.6)
    return matches[0] if matches else None


def _suggest_inverse_alternative(
    forward_predicate: str,
    known_predicates: list[str],
) -> str | None:
    """If a forward-direction predicate is listed where a g_to_pw inverse
    is needed, suggest the inverse_<name> form if it exists.

    DDKG convention: forward edges flow pathway -> gene; inverse_ forms
    flow gene -> pathway. The user wants the inverse form for membership.
    """
    candidate = f"inverse_{forward_predicate}"
    if candidate in known_predicates:
        return candidate
    return None


def _format_cross_validation_errors(
    error_message: str,
    config: "pipeline_config.PipelineConfig",
    bifo_predicates: dict,
    out,
) -> None:
    """Pretty-print a ValueError raised by cross-validation.

    The error text from the loader is structured as:
        "pipeline_config is inconsistent with <path>:
            - source 'X': <error 1>
            - source 'X': <error 2>
        "

    This function unpacks that, groups errors by source SAB, and adds
    'Did you mean' hints based on the loaded bifo_mapping predicates.
    """
    print(f"{MARK_FAIL} Cross-validation found problems:", file=out)
    print(file=out)

    # Parse out the per-source errors from the message body.
    lines = error_message.splitlines()
    error_lines = [line.lstrip(" -").rstrip() for line in lines if line.lstrip().startswith("-")]

    known_pred_names = list(bifo_predicates.keys())

    # Group by source SAB if present in the line.
    by_source: dict[str, list[str]] = {}
    for line in error_lines:
        # Look for the "source 'X':" pattern.
        sab = None
        if "source '" in line:
            try:
                sab = line.split("source '", 1)[1].split("'", 1)[0]
            except IndexError:
                sab = None
        key = sab if sab else "(unattributed)"
        by_source.setdefault(key, []).append(line)

    for sab, errors in by_source.items():
        print(f"  source {sab}:", file=out)
        for i, line in enumerate(errors, start=1):
            # Strip the redundant "source 'X':" prefix for readability
            # within this grouping.
            display = line
            prefix = f"source '{sab}': "
            if display.startswith(prefix):
                display = display[len(prefix):]
            print(f"    [{i}] {display}", file=out)

            # Add "Did you mean" hint where applicable.
            hint = _hint_for_error(line, sab, known_pred_names)
            if hint:
                print(f"        --> {hint}", file=out)
        print(file=out)

    print("To resolve:", file=out)
    print("  - Edit the affected files (see paths above).", file=out)
    print("  - Re-run this command to confirm the fix.", file=out)
    print(file=out)
    print(
        "For BIFO predicate direction conventions and the classification "
        "taxonomy,\nsee: docs/configuring_pathway_sources.md",
        file=out,
    )


def _hint_for_error(
    error_line: str,
    sab: str,
    known_predicates: list[str],
) -> str | None:
    """Return a 'Did you mean ...' style hint for an error line, or None.

    Two heuristics:
      1. Predicate not in bifo_mapping -> suggest closest spelling.
      2. Predicate has wrong direction -> suggest the inverse_ form
         if it exists.
    """
    # Extract the predicate name from the error line. The loader writes
    # errors in the form: "predicate 'X' is not defined ..." or
    # "predicate 'X' has direction ..."
    if "predicate '" not in error_line:
        return None
    try:
        pred = error_line.split("predicate '", 1)[1].split("'", 1)[0]
    except IndexError:
        return None

    # Case 1: typo (predicate missing from bifo_mapping)
    if "is not defined" in error_line:
        suggestion = _suggest_predicate_alternative(pred, known_predicates)
        if suggestion:
            return (
                f"Did you mean '{suggestion}'? "
                f"(closest known predicate in bifo_mapping)"
            )

    # Case 2: wrong direction (forward used where inverse needed)
    if "has direction 'target_to_source'" in error_line:
        suggestion = _suggest_inverse_alternative(pred, known_predicates)
        if suggestion:
            return (
                f"Did you mean '{suggestion}'? "
                f"(forward edges flow pathway -> gene; inverse_ forms "
                f"flow gene -> pathway, which is what BIFO-PPR uses)"
            )

    return None


def _print_success(
    config: "pipeline_config.PipelineConfig",
    bifo_predicates: dict,
    pipeline_config_path: Path,
    bifo_mapping_path: Path,
    out,
) -> None:
    """Pretty-print the all-checks-passed report."""
    print(
        f"{MARK_OK} pipeline_config.yaml schema valid "
        f"(schema_version: {config.schema_version})",
        file=out,
    )
    print(
        f"{MARK_OK} bifo_mapping_ddkg.yaml structure valid "
        f"({len(bifo_predicates)} predicates loaded)",
        file=out,
    )
    print(f"{MARK_OK} Cross-validation passed:", file=out)
    print(file=out)

    enabled = config.enabled_sources()
    if not enabled:
        print(
            "       (no pathway sources are enabled; nothing to "
            "cross-validate)",
            file=out,
        )
    for source in enabled:
        n = len(source.propagating_predicates_g_to_pw)
        print(
            f"       {source.sab} (enabled): {n} predicate(s), all "
            f"propagating in bifo_mapping",
            file=out,
        )
        for pred in source.propagating_predicates_g_to_pw:
            entry = bifo_predicates[pred]
            print(f"         - {pred}", file=out)
            print(
                f"           classification: {entry.classification}",
                file=out,
            )

    # List disabled sources for completeness so the user sees them.
    disabled = [
        s for s in config.pathway_sources.values() if not s.enabled
    ]
    if disabled:
        for source in disabled:
            print(
                f"       {source.sab} (disabled): not validated", file=out
            )

    print(file=out)
    print(f"{MARK_OK} Gene SABs: {', '.join(config.gene_sabs)}", file=out)
    n_excl = len(config.burden_control_exclusions)
    print(
        f"{MARK_OK} Burden control exclusions: {n_excl} entries", file=out
    )
    print(file=out)
    print(
        "All checks passed. Configuration is ready for pipeline run.",
        file=out,
    )


def check(
    pipeline_config_path: Path,
    bifo_mapping_path: Path,
    quiet: bool,
    out,
    err,
) -> int:
    """Run all checks. Returns an exit code.

    Args:
        pipeline_config_path: path to pipeline_config.yaml.
        bifo_mapping_path: path to bifo_mapping_ddkg.yaml.
        quiet: if True, suppress success-path output.
        out: stream for normal output.
        err: stream for error output.

    Returns:
        EXIT_OK on success; EXIT_CONFIG_ERROR on validation failure;
        EXIT_SETUP_ERROR on missing files / unreadable input.
    """
    # Setup checks: do the files exist?
    if not pipeline_config_path.exists():
        print(
            f"{MARK_FAIL} pipeline_config not found: {pipeline_config_path}",
            file=err,
        )
        return EXIT_SETUP_ERROR
    if not bifo_mapping_path.exists():
        print(
            f"{MARK_FAIL} bifo_mapping not found: {bifo_mapping_path}",
            file=err,
        )
        return EXIT_SETUP_ERROR

    if not quiet:
        _print_header(pipeline_config_path, bifo_mapping_path, out)

    # Run validation. If anything fails, the loader raises ValueError
    # with a descriptive message.
    try:
        # Load bifo_mapping first so we have the predicate inventory
        # for "did you mean" hints when we present cross-validation errors.
        bifo_predicates = pipeline_config._load_bifo_mapping_predicates(
            bifo_mapping_path
        )
    except ValueError as e:
        print(f"{MARK_FAIL} {e}", file=err)
        return EXIT_CONFIG_ERROR
    except FileNotFoundError as e:
        print(f"{MARK_FAIL} {e}", file=err)
        return EXIT_SETUP_ERROR

    try:
        config = pipeline_config.load_config(
            pipeline_config_path,
            bifo_mapping_path=bifo_mapping_path,
        )
    except ValueError as e:
        # Distinguish between schema errors (no cross-validation
        # happened yet) and cross-validation errors (we have predicates
        # loaded and can offer hints).
        msg = str(e)
        if "is inconsistent with" in msg:
            _format_cross_validation_errors(
                msg, config=None, bifo_predicates=bifo_predicates, out=err
            )
        else:
            print(f"{MARK_FAIL} {msg}", file=err)
        return EXIT_CONFIG_ERROR
    except FileNotFoundError as e:
        print(f"{MARK_FAIL} {e}", file=err)
        return EXIT_SETUP_ERROR

    # Success.
    if not quiet:
        _print_success(
            config,
            bifo_predicates,
            pipeline_config_path,
            bifo_mapping_path,
            out,
        )
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate bifo-graph configuration before pipeline run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--pipeline-config",
        type=Path,
        default=Path(DEFAULT_PIPELINE_CONFIG),
        help=f"path to pipeline_config.yaml (default: {DEFAULT_PIPELINE_CONFIG})",
    )
    parser.add_argument(
        "--bifo-mapping",
        type=Path,
        default=Path(DEFAULT_BIFO_MAPPING),
        help=f"path to bifo_mapping_ddkg.yaml (default: {DEFAULT_BIFO_MAPPING})",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress success-path output (errors still printed)",
    )
    args = parser.parse_args(argv)

    return check(
        pipeline_config_path=args.pipeline_config,
        bifo_mapping_path=args.bifo_mapping,
        quiet=args.quiet,
        out=sys.stdout,
        err=sys.stderr,
    )


if __name__ == "__main__":
    sys.exit(main())
