#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from logging_handler import Logger as CustomLogger


VERSION = "1.2.0"
FASTA_SUFFIXES = frozenset({".fasta", ".fa", ".fna"})

class DuplicateAccessionError(Exception):
    """Raised when more than one FASTA file matches the same accession."""


def parse_targets(targets_file: Path) -> list[str]:
    """
    Parse target accessions from a tab-delimited file.

    The first column is used and everything after the first dot is removed.

    Examples:
        GCF_123456.1      -> GCF_123456
        GNARF123.contig   -> GNARF123

    Empty lines are ignored.

    Args:
        targets_file: File containing target accessions.

    Returns:
        Ordered list of unique accessions.

    Raises:
        FileNotFoundError: If the target file does not exist.
        OSError: If the file cannot be read.
    """
    if not targets_file.is_file():
        raise FileNotFoundError(f"Target file not found: {targets_file}")

    with targets_file.open("r", encoding="utf-8") as handle:
        targets = [
            line.split("\t", 1)[0].strip().split(".", 1)[0]
            for line in handle
            if line.strip()
        ]

    # Deduplicate while preserving order.
    return list(dict.fromkeys(targets))


def is_fasta_file(file_path: Path) -> bool:
    """Return True if the path is a regular file with a supported FASTA suffix."""
    return file_path.is_file() and file_path.suffix.lower() in FASTA_SUFFIXES


def extract_accession(file_path: Path) -> str:
    """
    Extract the accession from a FASTA filename.

    The accession is defined as the part of the filename stem before the first
    dot.

    Examples:
        GCF_123456.1.fna  -> GCF_123456
        GNARF12.contig.fa -> GNARF12
        sample.fasta      -> sample
    """
    return file_path.stem.split(".", 1)[0]


def iter_fasta_files(source_dir: Path):
    """Yield FASTA files in the source directory (non-recursive)."""
    for file_path in source_dir.iterdir():
        if is_fasta_file(file_path):
            yield file_path


def build_accession_index(source_dir: Path) -> dict[str, list[Path]]:
    """
    Build a mapping from accession to matching FASTA files.

    Args:
        source_dir: Directory containing assembly files.

    Returns:
        Dictionary mapping accession -> list of FASTA files.
    """
    index: defaultdict[str, list[Path]] = defaultdict(list)

    for file_path in iter_fasta_files(source_dir):
        index[extract_accession(file_path)].append(file_path)

    return dict(index)


def copy_file(source: Path, destination_dir: Path) -> None:
    """Copy a file into the destination directory, preserving metadata."""
    shutil.copy2(source, destination_dir)


def resolve_single_match(
    accession: str,
    matches: list[Path],
    logger: logging.Logger,
) -> Path | None:
    """
    Resolve a single FASTA match for an accession.

    Returns:
        The matched file path, or None if no match exists.

    Raises:
        DuplicateAccessionError: If multiple files match the accession.
    """
    if not matches:
        logger.warning(
            "Accession '%s' not found. Are you in the folder with the "
            "assemblies? Please double-check!",
            accession,
        )
        return None

    if len(matches) > 1:
        raise DuplicateAccessionError(
            "Multiple assemblies found for accession "
            f"'{accession}': {[str(path) for path in matches]}. "
            "Please keep only one assembly (preferably the highest quality)."
        )

    return matches[0]


def copy_target_files(
    targets: list[str],
    accession_index: dict[str, list[Path]],
    target_dir: Path,
    logger: logging.Logger,
) -> None:
    """Copy all matched target FASTA files to the target directory."""
    for accession in targets:
        matches = accession_index.get(accession, [])

        logger.debug(
            "Matches for '%s': %s",
            accession,
            [str(path) for path in matches],
        )

        matched_file = resolve_single_match(accession, matches, logger)
        if matched_file is None:
            continue

        logger.debug("Copying %s to %s", matched_file, target_dir)
        copy_file(matched_file, target_dir)


def copy_neighbour_files(
    source_dir: Path,
    target_dir: Path,
    neighbour_dir: Path,
    logger: logging.Logger,
) -> None:
    """
    Copy all FASTA files not already present in the target directory into the
    neighbour directory.
    """
    for file_path in iter_fasta_files(source_dir):
        if (target_dir / file_path.name).exists():
            logger.debug(
                "File %s already exists in %s; skipping.",
                file_path.name,
                target_dir,
            )
            continue

        copy_file(file_path, neighbour_dir)
        logger.debug("Copied %s to %s", file_path.name, neighbour_dir)


def create_output_directories(base_dir: Path) -> tuple[Path, Path]:
    """
    Create and return the target and neighbour directories.

    Returns:
        Tuple of (target_dir, neighbour_dir).
    """
    target_dir = base_dir / "FUR.target"
    neighbour_dir = base_dir / "FUR.neighbour"

    target_dir.mkdir(parents=True, exist_ok=True)
    neighbour_dir.mkdir(parents=True, exist_ok=True)

    return target_dir, neighbour_dir


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="Dipper2",
        description="Target and Neighbour Folder Sorting",
        epilog=(
            "Bugs, suggestions, criticism, cake, capybaras and praise to "
            "t.wacker2@exeter.ac.uk"
        ),
    )

    parser.add_argument(
        "-t",
        "--target",
        type=Path,
        required=True,
        help="File containing target accessions.",
    )

    parser.add_argument(
        "-f",
        "--folder",
        type=Path,
        required=True,
        help=(
            "Results folder that will contain the "
            "FUR.target and FUR.neighbour subfolders."
        ),
    )

    parser.add_argument(
        "-V",
        "--verbose",
        action="store_true",
        help="Increase logging verbosity.",
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}",
    )

    return parser.parse_args()


def configure_logger(log_dir: Path, verbose: bool) -> logging.Logger:
    """
    Configure and return the application logger.

    The custom Logger wrapper is expected to expose a ``get_logger()`` method
    that returns an instance of ``logging.Logger``.
    """
    module_name = Path(__file__).name
    return CustomLogger(module_name, log_dir, verbose).get_logger()


def validate_paths(source_dir: Path, output_dir: Path) -> None:
    """
    Validate source and output paths.

    Raises:
        NotADirectoryError: If the source directory does not exist.
        ValueError: If the output directory is the same as the source directory.
    """
    if not source_dir.is_dir():
        raise NotADirectoryError(f"Source directory not found: {source_dir}")

    # Prevent accidental writes directly into the source directory.
    if output_dir.resolve() == source_dir.resolve():
        raise ValueError(
            "The output folder must not be the same as the source directory. "
            "Please specify a separate destination folder."
        )


def run(args: argparse.Namespace, logger: logging.Logger) -> None:
    """Execute the application workflow."""
    source_dir = Path.cwd()

    validate_paths(source_dir, args.folder)

    targets = parse_targets(args.target)
    target_dir, neighbour_dir = create_output_directories(args.folder)

    accession_index = build_accession_index(source_dir)

    copy_target_files(
        targets=targets,
        accession_index=accession_index,
        target_dir=target_dir,
        logger=logger,
    )

    copy_neighbour_files(
        source_dir=source_dir,
        target_dir=target_dir,
        neighbour_dir=neighbour_dir,
        logger=logger,
    )


def main() -> int:
    """Application entry point."""
    args = parse_arguments()
    logger = configure_logger(Path.cwd(), args.verbose)

    try:
        run(args, logger)
    except Exception as exc:
        logger.exception("Processing failed: %s", exc)
        return 1

    logger.info("Processing completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
