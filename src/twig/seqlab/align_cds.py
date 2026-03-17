#!/usr/bin/env python

"""Run MACSE codon-aware CDS alignment."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from loguru import logger

from twig.utils.exceptions import TwigError
from twig.utils.logger_setup import set_log_level


BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")


def _parse_fasta_records(inpath: Path) -> list[tuple[str, str]]:
    """Return FASTA records in input order."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq_chunks: list[str] = []

    with inpath.open("r") as datain:
        for raw_line in datain:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_chunks)))
                name = line[1:]
                seq_chunks = []
                continue
            if name is None:
                raise TwigError(
                    "malformed FASTA file: "
                    f"sequence data appears before a header in {inpath}"
                )
            seq_chunks.append(line)

    if name is not None:
        records.append((name, "".join(seq_chunks)))
    return records


def _load_seq_lr_selectors(seq_lr_args: Sequence[str]) -> list[str]:
    """Resolve less-reliable selectors from inline args or a selector file."""
    if len(seq_lr_args) == 1:
        candidate = Path(seq_lr_args[0])
        if candidate.exists() and candidate.is_file():
            selectors: list[str] = []
            with candidate.open("r") as datain:
                for raw_line in datain:
                    line = raw_line.strip()
                    if not line or line.startswith("#"):
                        continue
                    selectors.append(line)
            return selectors
    return list(seq_lr_args)


def _split_records_by_selectors(
    records: Sequence[tuple[str, str]],
    selectors: Sequence[str],
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Partition FASTA records into reliable and less-reliable sets."""
    try:
        patterns = [re.compile(pattern) for pattern in selectors]
    except re.error as exc:
        raise TwigError(f"invalid --seq-lr regex pattern: {exc}") from exc

    reliable_records: list[tuple[str, str]] = []
    less_reliable_records: list[tuple[str, str]] = []
    for name, seq in records:
        if any(pattern.search(name) for pattern in patterns):
            less_reliable_records.append((name, seq))
        else:
            reliable_records.append((name, seq))
    return reliable_records, less_reliable_records


def _write_fasta_records(records: Sequence[tuple[str, str]], outpath: Path) -> Path:
    """Write FASTA records to disk in their current order."""
    with outpath.open("w") as out:
        for name, seq in records:
            out.write(f">{name}\n{seq}\n")
    return outpath


def _get_tmp_aa_path(outpath: Path) -> Path:
    """Return the temporary amino-acid path MACSE writes alongside the CDS output."""
    return Path(f"{outpath}.tmp.aa.fa")


def build_macse_align_cmd(
    cds_fasta: Path,
    outpath: Path,
    max_iter: int,
    seq_lr_fasta: Path | None = None,
) -> list[str]:
    """Build the MACSE `alignSequences` command.

    Args:
        cds_fasta: FASTA file passed to MACSE as the reliable `-seq` input.
        outpath: Output path for the aligned nucleotide FASTA.
        max_iter: Maximum number of MACSE refinement iterations.
        seq_lr_fasta: Optional FASTA file passed to MACSE as `-seq_lr`.

    Returns:
        The command tokens to execute with ``subprocess.run``.
    """
    cmd = [
        BIN_MACSE,
        "-prog",
        "alignSequences",
        "-seq",
        str(cds_fasta),
    ]
    if seq_lr_fasta is not None:
        cmd.extend(["-seq_lr", str(seq_lr_fasta)])
    cmd.extend(
        [
            "-out_NT",
            str(outpath),
            "-out_AA",
            str(_get_tmp_aa_path(outpath)),
            "-max_refine_iter",
            str(max_iter),
        ]
    )
    return cmd


def call_macse_align(
    cds_fasta: Path,
    outpath: Path,
    max_iter: int,
    verbose: bool,
    seq_lr_fasta: Path | None = None,
) -> int:
    """Run MACSE `alignSequences` on reliable and optional less-reliable CDS inputs.

    Args:
        cds_fasta: FASTA file containing reliable CDS sequences.
        outpath: Path where the aligned nucleotide FASTA will be written.
        max_iter: Maximum number of MACSE refinement iterations.
        verbose: If ``True``, stream MACSE stderr to the terminal.
        seq_lr_fasta: Optional FASTA file containing less-reliable sequences for
            MACSE's ``-seq_lr`` input.

    Returns:
        The MACSE subprocess return code.

    Raises:
        TwigError: If the MACSE subprocess exits with a non-zero return code.
    """
    cmd = build_macse_align_cmd(
        cds_fasta=cds_fasta,
        outpath=outpath,
        max_iter=max_iter,
        seq_lr_fasta=seq_lr_fasta,
    )
    logger.debug(" ".join(cmd))
    if verbose:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=sys.stderr,
            text=True,
            check=False,
        )
    else:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if proc.returncode:
        raise TwigError((proc.stderr or "") + (proc.stdout or ""))
    _get_tmp_aa_path(outpath).unlink(missing_ok=True)
    return proc.returncode


def run_align_cds(args) -> int:
    """Run the `twig align-cds` workflow.

    Args:
        args: Parsed CLI arguments from ``get_parser_align_cds``.

    Returns:
        ``0`` when the alignment is skipped because output already exists or when
        MACSE completes successfully.

    Raises:
        IOError: If required input paths are missing or invalid.
        TwigError: If user-provided ``--seq-lr`` selectors are invalid or select
            every sequence in the input FASTA.
    """
    set_log_level(args.log_level)

    if not Path(BIN_MACSE).exists():
        raise TwigError(f"macse binary not found. Checked: {BIN_MACSE}")

    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    if args.outpath.is_dir():
        raise IOError("outpath must be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    if args.outpath.exists() and not args.force:
        logger.warning(
            f"[skipping] {args.outpath} already exists. Use --force to overwrite"
        )
        return 0

    if args.seq_lr:
        records = _parse_fasta_records(args.input)
        selectors = _load_seq_lr_selectors(args.seq_lr)
        reliable_records, less_reliable_records = _split_records_by_selectors(
            records,
            selectors,
        )

        if not less_reliable_records:
            logger.warning(
                "--seq-lr did not match any input records; continuing without "
                "-seq_lr"
            )
        elif not reliable_records:
            raise TwigError(
                "--seq-lr matched every input record; at least one reliable "
                "sequence must remain"
            )
        else:
            logger.info(
                f"[{args.outpath.name}] split {len(records)} records into "
                f"{len(reliable_records)} reliable and "
                f"{len(less_reliable_records)} less-reliable sequences"
            )
            with TemporaryDirectory(
                prefix=f"{args.outpath.name}.seq_lr.",
                dir=args.outpath.parent,
            ) as temp_dir:
                temp_dir_path = Path(temp_dir)
                reliable_fasta = _write_fasta_records(
                    reliable_records,
                    temp_dir_path / "reliable.nt.fa",
                )
                less_reliable_fasta = _write_fasta_records(
                    less_reliable_records,
                    temp_dir_path / "less_reliable.nt.fa",
                )
                call_macse_align(
                    reliable_fasta,
                    args.outpath,
                    args.max_refine_iter,
                    args.verbose,
                    seq_lr_fasta=less_reliable_fasta,
                )
                logger.info(
                    f"[{args.outpath.name}] alignment written to {args.outpath}"
                )
                return 0

    call_macse_align(args.input, args.outpath, args.max_refine_iter, args.verbose)
    logger.info(f"[{args.outpath.name}] alignment written to {args.outpath}")
    return 0


def main() -> int:
    """Run `align-cds` as a module entrypoint.

    Returns:
        The exit code returned by :func:`run_align_cds`.
    """
    from ..cli.subcommands import get_parser_align_cds

    parser = get_parser_align_cds()
    args = parser.parse_args()
    return run_align_cds(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
