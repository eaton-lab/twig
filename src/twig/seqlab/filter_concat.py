#!/usr/bin/env python

"""Filter aligned loci and concatenate a supermatrix."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import re

from loguru import logger

from twig.utils import TwigError, set_log_level


NT_CHARS = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv-?.")


@dataclass
class LocusData:
    name: str
    length: int
    seq_by_taxon: Dict[str, str]
    missing_char: str


def parse_fasta_alignment(path: Path) -> Dict[str, str]:
    """Return aligned fasta records as {header: sequence}."""
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
                    raise TwigError(f"{path}: empty FASTA header")
                if name in seqs:
                    raise TwigError(f"{path}: duplicate FASTA header: {name}")
                seqs[name] = []
            else:
                if name is None:
                    raise TwigError(f"{path}: FASTA starts with sequence before any header")
                seqs[name].append(line)
    records = {k: "".join(v) for k, v in seqs.items()}
    if not records:
        raise TwigError(f"{path}: no FASTA records")
    lengths = {len(s) for s in records.values()}
    if len(lengths) != 1:
        preview = ",".join(str(i) for i in sorted(lengths)[:6])
        raise TwigError(f"{path}: alignment has variable sequence lengths ({preview})")
    return records


def parse_imap(path: Path) -> Tuple[List[str], Dict[str, str]]:
    """Parse IMAP TSV rows: <name> <pop>."""
    order: List[str] = []
    popmap: Dict[str, str] = {}
    with open(path, "r") as data:
        for idx, line in enumerate(data, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise TwigError(f"{path}:{idx}: expected exactly 2 columns (<name> <pop>)")
            name, pop = parts[0], parts[1]
            if name in popmap:
                raise TwigError(f"{path}:{idx}: duplicate IMAP name: {name}")
            order.append(name)
            popmap[name] = pop
    if not order:
        raise TwigError(f"{path}: no IMAP rows found")
    return order, popmap


def parse_minmap(path: Path) -> Dict[str, int]:
    """Parse MINMAP TSV rows: <imap-pop> <min-samples>."""
    minmap: Dict[str, int] = {}
    with open(path, "r") as data:
        for idx, line in enumerate(data, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise TwigError(f"{path}:{idx}: expected exactly 2 columns (<imap-pop> <min-samples>)")
            pop, mincov = parts
            try:
                mincov_i = int(mincov)
            except ValueError as exc:
                raise TwigError(f"{path}:{idx}: min-samples must be an integer") from exc
            if mincov_i < 0:
                raise TwigError(f"{path}:{idx}: min-samples must be >= 0")
            if pop in minmap:
                raise TwigError(f"{path}:{idx}: duplicate population in minmap: {pop}")
            minmap[pop] = mincov_i
    if not minmap:
        raise TwigError(f"{path}: no MINMAP rows found")
    return minmap


def parse_taxon(name: str, delim: str | None, delim_idxs: List[int], delim_join: str) -> str:
    """Parse taxon label from a sequence header."""
    if delim is None:
        return name
    parts = name.split(delim)
    if not parts:
        return name
    try:
        kept = [parts[i] for i in delim_idxs]
    except IndexError as exc:
        raise TwigError(
            f"name '{name}' cannot be parsed by --delim '{delim}' and --delim-idxs {delim_idxs}"
        ) from exc
    parsed = delim_join.join(kept)
    return parsed if parsed else name


def scan_loci_taxa(paths: List[Path], delim: str | None, delim_idxs: List[int], delim_join: str) -> Tuple[List[str], int]:
    """Fast first pass through headers only to gather taxa and max per-locus richness."""
    taxa_ordered: "OrderedDict[str, None]" = OrderedDict()
    max_taxa = 0
    for path in paths:
        seen = set()
        with open(path, "r") as data:
            for line in data:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                taxon = parse_taxon(header, delim, delim_idxs, delim_join)
                seen.add(taxon)
                taxa_ordered.setdefault(taxon, None)
        max_taxa = max(max_taxa, len(seen))
    return list(taxa_ordered.keys()), max_taxa


def is_nucleotide_like(seqs: List[str]) -> bool:
    chars = set("".join(seqs))
    return chars.issubset(NT_CHARS)


def sanitize_for_phylip(name: str) -> str:
    s = re.sub(r"\s+", "_", name.strip())
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "taxon"


def sanitize_partition_name(name: str) -> str:
    s = re.sub(r"\s+", "_", name.strip())
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "locus"


def build_partitions(loci: List[LocusData]) -> List[Tuple[str, int, int]]:
    """Return ordered list of (partition_name, start, end) in 1-based coords."""
    parts: List[Tuple[str, int, int]] = []
    start = 1
    used: Dict[str, int] = {}
    for locus in loci:
        pname = sanitize_partition_name(locus.name)
        if pname in used:
            used[pname] += 1
            pname = f"{pname}__{used[pname]}"
        else:
            used[pname] = 0
        end = start + locus.length - 1
        parts.append((pname, start, end))
        start = end + 1
    return parts


def write_fasta(path: Path, matrix: Dict[str, str], taxa_order: List[str]) -> None:
    with open(path, "w") as out:
        for taxon in taxa_order:
            out.write(f">{taxon}\n{matrix[taxon]}\n")


def write_phylip(path: Path, matrix: Dict[str, str], taxa_order: List[str]) -> None:
    ntaxa = len(taxa_order)
    nchar = len(matrix[taxa_order[0]]) if taxa_order else 0
    used: Dict[str, int] = {}
    labels: List[str] = []
    for taxon in taxa_order:
        label = sanitize_for_phylip(taxon)
        if label in used:
            used[label] += 1
            label = f"{label}__{used[label]}"
        else:
            used[label] = 0
        labels.append(label)
    with open(path, "w") as out:
        out.write(f"{ntaxa} {nchar}\n")
        for label, taxon in zip(labels, taxa_order):
            out.write(f"{label} {matrix[taxon]}\n")


def write_nexus(path: Path, matrix: Dict[str, str], taxa_order: List[str]) -> None:
    ntaxa = len(taxa_order)
    nchar = len(matrix[taxa_order[0]]) if taxa_order else 0
    concat = [matrix[t] for t in taxa_order]
    datatype = "DNA" if is_nucleotide_like(concat) else "PROTEIN"
    with open(path, "w") as out:
        out.write("#NEXUS\n")
        out.write("BEGIN DATA;\n")
        out.write(f"  DIMENSIONS NTAX={ntaxa} NCHAR={nchar};\n")
        out.write(f"  FORMAT DATATYPE={datatype} MISSING=? GAP=-;\n")
        out.write("  MATRIX\n")
        for taxon in taxa_order:
            out.write(f"  {taxon} {matrix[taxon]}\n")
        out.write("  ;\n")
        out.write("END;\n")


def write_partition_file(path: Path, partitions: List[Tuple[str, int, int]]) -> None:
    """Write a simple per-locus partition range file."""
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "w") as out:
        for pname, start, end in partitions:
            out.write(f"{pname} = {start}-{end}\n")


def run_filter_concat(args) -> None:
    """Run filter-concat subcommand."""
    set_log_level(args.log_level)

    if args.imap and not args.imap.exists():
        raise TwigError(f"{args.imap} not found")
    if args.minmap and not args.minmap.exists():
        raise TwigError(f"{args.minmap} not found")
    for path in args.input:
        if not path.exists():
            raise TwigError(f"{path} not found")
    if args.output.is_dir():
        raise TwigError("output must be a file path, not a directory")
    args.output.parent.mkdir(exist_ok=True, parents=True)
    if args.output.exists() and not args.force:
        raise TwigError(f"{args.output} exists; use --force to overwrite")
    if args.partition_file and args.partition_file.exists() and not args.force:
        raise TwigError(f"{args.partition_file} exists; use --force to overwrite")
    if args.min_present < 1:
        raise TwigError("--min-present must be >= 1")
    if args.minmap and not args.imap:
        raise TwigError("--minmap requires --imap")

    if args.imap:
        taxa_order, imap = parse_imap(args.imap)
        max_locus_taxa = None
    else:
        imap = None
        taxa_order, max_locus_taxa = scan_loci_taxa(args.input, args.delim, args.delim_idxs, args.delim_join)
        if max_locus_taxa > 1000:
            logger.warning(
                f"found locus with {max_locus_taxa} unique taxa; this may indicate a taxon-name delimiting error"
            )

    pop_minmap = parse_minmap(args.minmap) if args.minmap else None
    if pop_minmap and imap:
        pops_in_imap = set(imap.values())
        diff = set(pop_minmap) - pops_in_imap
        if diff:
            raise TwigError(f"populations in --minmap not present in --imap: {sorted(diff)}")

    loci: List[LocusData] = []
    excluded_multi = 0
    excluded_min = 0
    total = 0

    for locus_path in args.input:
        total += 1
        records = parse_fasta_alignment(locus_path)
        length = len(next(iter(records.values())))

        parsed: Dict[str, str] = {}
        duplicate = False
        for header, seq in records.items():
            taxon = parse_taxon(header, args.delim, args.delim_idxs, args.delim_join)
            if taxon in parsed:
                duplicate = True
                break
            parsed[taxon] = seq
        if duplicate:
            excluded_multi += 1
            continue

        if imap:
            present = {tax: seq for tax, seq in parsed.items() if tax in imap}
            if pop_minmap:
                counts: Dict[str, int] = {}
                for tax in present:
                    pop = imap[tax]
                    counts[pop] = counts.get(pop, 0) + 1
                failed = any(counts.get(pop, 0) < mincov for pop, mincov in pop_minmap.items())
                if failed:
                    excluded_min += 1
                    continue
            else:
                groups_present = {imap[tax] for tax in present}
                if len(groups_present) < args.min_present:
                    excluded_min += 1
                    continue
        else:
            present = parsed
            if len(present) < args.min_present:
                excluded_min += 1
                continue

        missing_char = "N" if is_nucleotide_like(list(present.values())) else "X"
        loci.append(
            LocusData(
                name=locus_path.stem,
                length=length,
                seq_by_taxon=present,
                missing_char=missing_char,
            )
        )

    if not loci:
        raise TwigError(
            f"no loci retained (total={total}, excluded_multi_copy={excluded_multi}, "
            f"excluded_min_imap_samples={excluded_min})"
        )

    matrix_parts: Dict[str, List[str]] = {tax: [] for tax in taxa_order}
    for locus in loci:
        for taxon in taxa_order:
            seq = locus.seq_by_taxon.get(taxon)
            if seq is None:
                matrix_parts[taxon].append(locus.missing_char * locus.length)
            else:
                matrix_parts[taxon].append(seq)
    matrix = {k: "".join(v) for k, v in matrix_parts.items()}
    partitions = build_partitions(loci)

    lengths = {len(v) for v in matrix.values()}
    if len(lengths) != 1:
        raise TwigError("internal error: concatenated matrix has variable row lengths")
    nsites = lengths.pop()
    ntaxa = len(matrix)

    if args.format == "fasta":
        write_fasta(args.output, matrix, taxa_order)
    elif args.format == "phylip":
        write_phylip(args.output, matrix, taxa_order)
    elif args.format == "nexus":
        write_nexus(args.output, matrix, taxa_order)
    else:
        raise TwigError(f"unsupported format: {args.format}")
    if args.partition_file:
        write_partition_file(args.partition_file, partitions)

    logger.info(
        f"retained {len(loci)}/{total} loci (excluded_multi_copy={excluded_multi}, "
        f"excluded_min_threshold={excluded_min})"
    )
    logger.info(f"supermatrix dimensions: ntaxa={ntaxa}, nsites={nsites}")
    logger.info(f"concatenated alignment written to {args.output}")
    if args.partition_file:
        logger.info(f"partition ranges written to {args.partition_file}")


def main() -> None:
    from twig.cli.subcommands import get_parser_filter_concat

    parser = get_parser_filter_concat()
    args = parser.parse_args()
    run_filter_concat(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
