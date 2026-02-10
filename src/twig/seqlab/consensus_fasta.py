#!/usr/bin/env python

"""Write a consensus sequence from an aligned FASTA."""

from collections import Counter
from pathlib import Path
from typing import Dict, List

from loguru import logger

from twig.utils import TwigError, set_log_level


def parse_fasta(path: Path) -> Dict[str, str]:
    """Parse FASTA file into {name: sequence}."""
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path, "r") as data:
        for line in data:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                if not name:
                    raise TwigError("encountered empty FASTA header")
                if name in seqs:
                    raise TwigError(f"duplicate FASTA header: {name}")
                seqs[name] = []
            else:
                if name is None:
                    raise TwigError("FASTA starts with sequence before any header")
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def infer_ambiguous_char(seqs: Dict[str, str]) -> str:
    """Infer default ambiguous char based on sequence alphabet."""
    chars = set("".join(seqs.values()))
    dna_chars = set("ACGTUNacgtun-?.")
    if chars and chars.issubset(dna_chars):
        return "N"
    return "X"


def build_consensus(
    seqs: Dict[str, str],
    ignore_chars: set[str],
    ambiguous_char: str,
    min_majority: float,
) -> str:
    """Build consensus sequence across aligned input sequences."""
    if not seqs:
        raise TwigError("input fasta is empty")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise TwigError("input contains variable sequence lengths (not aligned)")
    seqlen = lengths.pop()
    out = []

    values = list(seqs.values())
    for idx in range(seqlen):
        col = [s[idx] for s in values]
        votes = [c for c in col if c not in ignore_chars]
        if not votes:
            out.append(ambiguous_char)
            continue

        counts = Counter(votes)
        top = counts.most_common(2)
        best_char, best_count = top[0]
        frac = best_count / len(votes)

        tie = len(top) > 1 and top[1][1] == best_count
        if tie or frac < min_majority:
            out.append(ambiguous_char)
        else:
            out.append(best_char.upper())
    return "".join(out)


def run_consensus_fasta(args) -> None:
    """Run consensus-fasta subcommand."""
    set_log_level(args.log_level)

    if not (args.input.exists() and args.input.is_file()):
        raise TwigError(f"{args.input} not found")
    if args.output.is_dir():
        raise TwigError("output must be a file path, not a directory")
    args.output.parent.mkdir(exist_ok=True, parents=True)

    if args.output.exists() and not args.force:
        logger.warning(f"[skipping] {args.output} already exists. Use --force to overwrite")
        return

    if len(args.name.strip()) == 0:
        raise TwigError("--name cannot be empty")
    if len(args.ambiguous_char) > 1 if args.ambiguous_char else False:
        raise TwigError("--ambiguous-char must be a single character")
    if not (0.0 < args.min_majority <= 1.0):
        raise TwigError("--min-majority must be in (0, 1]")

    seqs = parse_fasta(args.input)
    ambig = args.ambiguous_char if args.ambiguous_char else infer_ambiguous_char(seqs)
    consensus = build_consensus(
        seqs=seqs,
        ignore_chars=set(args.ignore_chars),
        ambiguous_char=ambig,
        min_majority=args.min_majority,
    )

    with open(args.output, "w") as out:
        out.write(f">{args.name}\n{consensus}\n")
    logger.info(f"consensus written to {args.output}")


def main() -> None:
    from twig.cli.subcommands import get_parser_consensus_fasta

    parser = get_parser_consensus_fasta()
    args = parser.parse_args()
    run_consensus_fasta(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
