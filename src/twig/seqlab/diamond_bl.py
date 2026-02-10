#!/usr/bin/env python

"""Run DIAMOND BLASTP with optional filtering and annotation joins."""

import csv
import gzip
import shutil
import sys
from collections import defaultdict
from subprocess import run, PIPE
from pathlib import Path
from loguru import logger
from ..utils.logger_setup import set_log_level

DEFAULT_EXTENDED_FIELDS = [
    "qseqid", "sseqid",
    "pident", "length",
    "mismatch", "gapopen",
    "qstart", "qend",
    "sstart", "send",
    "evalue", "bitscore",
    "corrected_bitscore", "qlen", "slen",
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


def load_annotations(path: Path, key_col: str, wanted_cols: list[str] | None) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Load subject annotation TSV keyed by key_col."""
    table: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or key_col not in reader.fieldnames:
            raise ValueError(f"annotation file must contain a '{key_col}' column: {path}")
        if wanted_cols:
            for col in wanted_cols:
                if col not in reader.fieldnames:
                    raise ValueError(f"annotation column '{col}' not found in {path}")
            cols = wanted_cols
        else:
            cols = [i for i in reader.fieldnames if i != key_col]
        for row in reader:
            key = row.get(key_col)
            if key:
                table[key] = {col: row.get(col, "") for col in cols}
    return table, cols


def filter_and_rank_rows(
    rows: list[dict[str, str]],
    min_qcov: float | None,
    min_scov: float | None,
    best_hit_only: bool,
    top_hits: int | None,
) -> list[dict[str, str]]:
    """Apply optional coverage and per-query hit ranking filters."""
    if min_qcov is not None:
        rows = [row for row in rows if parse_float(row["qcovhsp"], "qcovhsp") >= min_qcov]
    if min_scov is not None:
        rows = [row for row in rows if parse_float(row["scovhsp"], "scovhsp") >= min_scov]

    if not best_hit_only and not top_hits:
        return rows

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for idx, row in enumerate(rows):
        row["_row_idx"] = str(idx)
        grouped[row["qseqid"]].append(row)

    keep_n = 1 if best_hit_only else int(top_hits)
    out: list[dict[str, str]] = []
    for qid in grouped:
        candidates = grouped[qid]
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
    annots: dict[str, dict[str, str]] | None,
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
        if annots and annot_cols:
            match = annots.get(row.get("sseqid", ""), {})
            values.extend(match.get(col, "") for col in annot_cols)
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


def run_diamond_bl(args):
    """Runs diamond blast with args parsed from sys or CLI."""
    set_log_level(args.log_level)

    if not args.query.exists():
        raise ValueError(f"query FASTA not found: {args.query}")
    if not args.target.exists():
        raise ValueError(f"target path not found: {args.target}")

    if args.min_qcov is not None and not (0.0 <= args.min_qcov <= 100.0):
        raise ValueError("--min-qcov must be between 0 and 100")
    if args.min_scov is not None and not (0.0 <= args.min_scov <= 100.0):
        raise ValueError("--min-scov must be between 0 and 100")
    if args.top_hits is not None and args.top_hits < 1:
        raise ValueError("--top-hits must be >= 1")

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

        rows = filter_and_rank_rows(rows, args.min_qcov, args.min_scov, args.best_hit_only, args.top_hits)

        annotations = None
        annotation_cols = None
        if args.subject_annotations is not None:
            if not args.subject_annotations.exists():
                raise ValueError(f"annotation TSV not found: {args.subject_annotations}")
            annotations, annotation_cols = load_annotations(args.subject_annotations, args.annotation_key, args.annotation_cols)

        write_output(
            rows=rows,
            output_fields=output_fields,
            outpath=args.out,
            compress=args.compress,
            header=args.header,
            annots=annotations,
            annot_cols=annotation_cols,
        )
    finally:
        if built_db and dbfile.exists() and (not args.save_db):
            dbfile.unlink()

    qids = {row.get("qseqid", "") for row in rows}
    sids = {row.get("sseqid", "") for row in rows}
    logger.info(f"wrote {len(rows)} hits ({len(qids)} query IDs, {len(sids)} subject IDs)")
    logger.complete()


def main():
    from ..cli.subcommands import get_parser_diamond_bl
    parser = get_parser_diamond_bl()
    args = parser.parse_args()
    run_diamond_bl(args)


if __name__ == "__main__":
    main()
