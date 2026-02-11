#!/usr/bin/env python

"""Run DIAMOND BLASTP with optional filtering and annotation joins."""

import gzip
import shutil
import sys
from collections import defaultdict
from subprocess import run, PIPE
from pathlib import Path
from loguru import logger
from ..utils.logger_setup import set_log_level
from .diamond_bl_annotate import detect_annotation_format, load_annotations, expand_rows_with_annotations

DEFAULT_EXTENDED_FIELDS = [
    "qseqid", "sseqid",
    "pident", "length",
    "mismatch", "gapopen",
    "qstart", "qend",
    "sstart", "send",
    "evalue", "bitscore",
    "corrected_bitscore", "qlen", "slen",
]

CONSENSUS_GO_FIELDS = [
    "consensus_qseqid",
    "go_id",
    "support_n_queries",
    "n_queries_total",
    "support_frac_queries",
    "passes_min_support",
]


def split_raw_kwargs(kwargs: list[str] | None) -> list[str]:
    """Split raw CLI kwargs into tokenized arg lists."""
    if not kwargs:
        return []
    out = []
    for arg in kwargs:
        if " " in arg:
            out.extend(arg.split())
        else:
            out.append(arg)
    return out


def resolve_diamond_binary(binary: Path | None) -> Path:
    """Resolve DIAMOND binary from user arg, PATH, or env prefix fallback."""
    if binary is not None:
        if binary.exists():
            return binary
        raise ValueError(f"DIAMOND binary not found: {binary}")

    from_path = shutil.which("diamond")
    if from_path:
        return Path(from_path)

    fallback = Path(sys.prefix) / "bin" / "diamond"
    if fallback.exists():
        return fallback
    raise ValueError("DIAMOND binary not found; set --binary or install diamond in $PATH")


def parse_float(value: str, field: str) -> float:
    """Parse float fields while allowing sparse TSV rows."""
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"cannot parse numeric field '{field}' from value: {value!r}") from None


def is_db_target(path: Path, target_is_db: bool) -> bool:
    """Infer whether target is already a DIAMOND DB."""
    return target_is_db or path.suffix in {".dmnd", ".db"}


def make_db_prefix(target_fasta: Path, tmpdir: Path) -> Path:
    """Set deterministic DB prefix for auto-created DBs."""
    return tmpdir / f"{target_fasta.name}.diamond"


def call_diamond_makedb(
    diamond_bin: Path,
    fasta: Path,
    db_prefix: Path,
    threads: int,
    kwargs: list[str] | None,
) -> Path:
    """Create a DIAMOND database file from target FASTA and return .dmnd path."""
    cmd = [
        str(diamond_bin), "makedb",
        *split_raw_kwargs(kwargs),
        "--ignore-warnings",
        "--in", str(fasta),
        "--db", str(db_prefix),
        "--threads", str(threads),
    ]

    logger.debug(" ".join(cmd))
    proc = run(cmd, stdout=PIPE, stderr=PIPE, text=True)
    if proc.returncode:
        raise ValueError(proc.stderr.strip() or proc.stdout.strip() or "diamond makedb failed")

    db_path = Path(str(db_prefix) + ".dmnd")
    if not db_path.exists():
        raise ValueError(f"diamond makedb completed but DB file was not found: {db_path}")
    return db_path


def run_blastp(
    diamond_bin: Path,
    query: Path,
    db_path: Path,
    threads: int,
    evalue: float,
    bitscore: float | None,
    max_target_seqs: int,
    ultra_sensitive: bool,
    internal_fields: list[str],
    kwargs: list[str] | None,
) -> str:
    """Run DIAMOND BLASTP and return TSV text."""
    cmd = [
        str(diamond_bin), "blastp",
        *split_raw_kwargs(kwargs),
        "--ignore-warnings",
        "--quiet",
        "--query", str(query),
        "--db", str(db_path),
        "--threads", str(threads),
        "--max-target-seqs", str(max_target_seqs),
        "--outfmt", "6", *internal_fields,
    ]

    if bitscore is not None:
        cmd += ["--min-score", str(bitscore)]
    else:
        cmd += ["--evalue", str(evalue)]

    if ultra_sensitive:
        cmd += ["--ultra-sensitive"]
    logger.debug(" ".join(cmd))

    proc = run(cmd, stdout=PIPE, stderr=PIPE, text=True)
    if proc.returncode != 0:
        err = proc.stderr.strip() or proc.stdout.strip() or "diamond blastp failed"
        raise ValueError(err)
    return proc.stdout


def filter_and_rank_rows(
    rows: list[dict[str, str]],
    min_qcov: float | None,
    min_scov: float | None,
    best_hit_only: bool,
    top_hits: int | None,
) -> list[dict[str, str]]:
    """Apply optional coverage and per-query hit ranking filters."""
    if min_qcov is not None:
        # DIAMOND reports qcovhsp as percent; CLI threshold is a 0-1 fraction.
        min_qcov_pct = min_qcov * 100.0
        rows = [row for row in rows if parse_float(row["qcovhsp"], "qcovhsp") >= min_qcov_pct]
    if min_scov is not None:
        # DIAMOND reports scovhsp as percent; CLI threshold is a 0-1 fraction.
        min_scov_pct = min_scov * 100.0
        rows = [row for row in rows if parse_float(row["scovhsp"], "scovhsp") >= min_scov_pct]

    if not best_hit_only and not top_hits:
        return rows

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for idx, row in enumerate(rows):
        # Keep input order as a deterministic tie-breaker after sorting.
        row["_row_idx"] = str(idx)
        grouped[row["qseqid"]].append(row)

    keep_n = 1 if best_hit_only else int(top_hits)
    out: list[dict[str, str]] = []
    for qid in grouped:
        candidates = grouped[qid]
        # Rank within each query by score, then evalue, then original row order.
        candidates.sort(
            key=lambda r: (
                -parse_float(r["bitscore"], "bitscore"),
                parse_float(r["evalue"], "evalue"),
                int(r["_row_idx"]),
            )
        )
        out.extend(candidates[:keep_n])
    out.sort(key=lambda r: int(r["_row_idx"]))
    for row in out:
        row.pop("_row_idx", None)
    return out


def write_output(
    rows: list[dict[str, str]],
    output_fields: list[str],
    outpath: Path | None,
    compress: bool,
    header: bool,
    annot_cols: list[str] | None,
) -> None:
    """Write TSV output to stdout or file, optionally gzip compressed."""
    columns = list(output_fields)
    if annot_cols:
        columns.extend(annot_cols)

    lines = []
    if header:
        lines.append("\t".join(columns))
    for row in rows:
        values = [row.get(col, "") for col in output_fields]
        if annot_cols:
            values.extend(row.get(col, "") for col in annot_cols)
        lines.append("\t".join(values))
    text = ("\n".join(lines) + "\n") if lines else ""

    if compress and outpath is None:
        raise ValueError("--compress requires --out path")

    if outpath is None:
        sys.stdout.write(text)
        return

    outpath.parent.mkdir(parents=True, exist_ok=True)
    if compress:
        with gzip.open(outpath, "wt") as out:
            out.write(text)
    else:
        outpath.write_text(text)


def build_go_consensus_rows(
    rows: list[dict[str, str]],
    consensus_id: str,
    n_queries_total: int,
    min_support: float,
    include_failed: bool,
) -> list[dict[str, str]]:
    """Aggregate GO-term support across distinct query IDs."""
    votes: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        go_id = row.get("go_id", "")
        qseqid = row.get("qseqid", "")
        if go_id and qseqid:
            # Count support per distinct query ID to avoid over-weighting duplicated rows.
            votes[go_id].add(qseqid)

    out: list[dict[str, str]] = []
    for go_id, qids in votes.items():
        support_n = len(qids)
        support_frac = (support_n / n_queries_total) if n_queries_total else 0.0
        passes = support_frac >= min_support
        if not include_failed and not passes:
            continue
        out.append(
            {
                "consensus_qseqid": consensus_id,
                "go_id": go_id,
                "support_n_queries": str(support_n),
                "n_queries_total": str(n_queries_total),
                "support_frac_queries": f"{support_frac:.6f}",
                "passes_min_support": "true" if passes else "false",
            }
        )

    out.sort(
        key=lambda r: (
            -int(r["support_n_queries"]),
            -float(r["support_frac_queries"]),
            r["go_id"],
        )
    )
    return out


def run_diamond_bl(args):
    """Runs diamond blast with args parsed from sys or CLI."""
    set_log_level(args.log_level)

    if not args.query.exists():
        raise ValueError(f"query FASTA not found: {args.query}")
    if not args.target.exists():
        raise ValueError(f"target path not found: {args.target}")

    if args.min_qcov is not None and not (0.0 <= args.min_qcov <= 1.0):
        raise ValueError("--min-qcov must be between 0 and 1")
    if args.min_scov is not None and not (0.0 <= args.min_scov <= 1.0):
        raise ValueError("--min-scov must be between 0 and 1")
    if args.top_hits is not None and args.top_hits < 1:
        raise ValueError("--top-hits must be >= 1")
    if args.consensus_go_min_support is not None and not (0.0 <= args.consensus_go_min_support <= 1.0):
        raise ValueError("--consensus-go-min-support must be between 0 and 1")
    if args.consensus_go_out is not None and not args.consensus_go_id:
        raise ValueError("--consensus-go-id is required when --consensus-go-out is set")

    diamond_bin = resolve_diamond_binary(args.binary)
    args.tmpdir.mkdir(parents=True, exist_ok=True)

    if args.outfmt is None:
        output_fields = list(DEFAULT_EXTENDED_FIELDS)
    elif len(args.outfmt) == 0:
        output_fields = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        ]
    else:
        output_fields = list(args.outfmt)
    internal_fields = list(output_fields)

    # Ensure fields needed for ranking/filtering are always present internally.
    for needed in ("qseqid", "sseqid", "bitscore", "evalue"):
        if needed not in internal_fields:
            internal_fields.append(needed)
    if args.min_qcov is not None and "qcovhsp" not in internal_fields:
        internal_fields.append("qcovhsp")
    if args.min_scov is not None and "scovhsp" not in internal_fields:
        internal_fields.append("scovhsp")

    dbfile: Path
    built_db = False
    if is_db_target(args.target, args.target_is_db):
        dbfile = args.target
        if not dbfile.exists():
            raise ValueError(f"DIAMOND database not found: {dbfile}")
    else:
        db_prefix = make_db_prefix(args.target, args.tmpdir)
        dbfile = Path(str(db_prefix) + ".dmnd")
        if dbfile.exists():
            if args.force:
                dbfile.unlink()
            else:
                raise ValueError(f"database already exists: {dbfile}. Use --force to overwrite.")
        dbfile = call_diamond_makedb(diamond_bin, args.target, db_prefix, args.threads, args.kwargs_makedb)
        built_db = True

    rows: list[dict[str, str]] = []
    ranked_rows: list[dict[str, str]] = []
    consensus_log_msg: str | None = None
    try:
        blast_text = run_blastp(
            diamond_bin=diamond_bin,
            query=args.query,
            db_path=dbfile,
            threads=args.threads,
            evalue=args.evalue,
            bitscore=args.bitscore,
            max_target_seqs=args.max_target_seqs,
            ultra_sensitive=args.ultra_sensitive,
            internal_fields=internal_fields,
            kwargs=args.kwargs_blastp,
        )

        if blast_text:
            for line in blast_text.rstrip("\n").splitlines():
                vals = line.split("\t")
                if len(vals) != len(internal_fields):
                    raise ValueError(f"unexpected DIAMOND output column count ({len(vals)}), expected {len(internal_fields)}")
                rows.append(dict(zip(internal_fields, vals)))

        # Coverage and per-query ranking happen before annotation expansion.
        rows = filter_and_rank_rows(rows, args.min_qcov, args.min_scov, args.best_hit_only, args.top_hits)
        # Keep pre-annotation rows for consensus denominator and query counting.
        ranked_rows = [dict(row) for row in rows]

        annotations = None
        annotation_cols = None
        annotation_match_mode = "exact"
        if args.target_annotations is not None:
            if not args.target_annotations.exists():
                raise ValueError(f"annotation file not found: {args.target_annotations}")
            annotation_format = detect_annotation_format(args.target_annotations, args.annotation_format)
            if annotation_format in {"gaf", "gpa"}:
                # DIAMOND sseqid may include isoform suffix (e.g., ".1"), while GO keys may be gene-level.
                annotation_match_mode = "exact_then_gene_fallback"
            annotations, annotation_cols = load_annotations(
                path=args.target_annotations,
                annotation_format=args.annotation_format,
                key_col=args.annotation_key,
                wanted_cols=args.annotation_cols,
            )

        rows = expand_rows_with_annotations(
            rows,
            annotations,
            annotation_cols,
            match_mode=annotation_match_mode,
        )

        write_output(
            rows=rows,
            output_fields=output_fields,
            outpath=args.out,
            compress=args.compress,
            header=args.header,
            annot_cols=annotation_cols,
        )

        if args.consensus_go_out is not None:
            if args.target_annotations is None:
                raise ValueError("--consensus-go-out requires --target-annotations with GO annotations")
            if not annotation_cols or "go_id" not in annotation_cols:
                raise ValueError("--consensus-go-out requires GO annotations with 'go_id' available")

            n_queries_total = len({row.get("qseqid", "") for row in ranked_rows if row.get("qseqid", "")})
            consensus_rows = build_go_consensus_rows(
                rows=rows,
                consensus_id=args.consensus_go_id,
                n_queries_total=n_queries_total,
                min_support=args.consensus_go_min_support,
                include_failed=args.consensus_go_include_failed,
            )
            # Consensus output is intentionally separate from the per-hit output TSV.
            write_output(
                rows=consensus_rows,
                output_fields=CONSENSUS_GO_FIELDS,
                outpath=args.consensus_go_out,
                compress=args.consensus_go_out.suffix == ".gz",
                header=True,
                annot_cols=None,
            )
            consensus_log_msg = (
                f"wrote {len(consensus_rows)} consensus GO rows for {n_queries_total} query IDs to {args.consensus_go_out}"
            )
    finally:
        if built_db and dbfile.exists() and (not args.save_db):
            dbfile.unlink()

    qids = {row.get("qseqid", "") for row in rows}
    sids = {row.get("sseqid", "") for row in rows}
    logger.info(f"wrote {len(rows)} hits ({len(qids)} query IDs, {len(sids)} subject IDs)")
    if consensus_log_msg:
        logger.info(consensus_log_msg)
    logger.complete()


def main():
    from ..cli.subcommands import get_parser_diamond_bl
    parser = get_parser_diamond_bl()
    args = parser.parse_args()
    run_diamond_bl(args)


if __name__ == "__main__":
    main()
