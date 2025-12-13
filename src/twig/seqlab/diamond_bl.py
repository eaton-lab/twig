#!/usr/bin/env python

"""...

orthogroup synteny
evo-analysis-synteny-orthology
synortho

$ twig diamond-blast A.faa B.faa -e 1e-5 -k 1 > A_B.tsv

TODO
----
--comment | option to add comment line with diamond command parameters

"""

import sys
from subprocess import Popen, STDOUT, PIPE
from pathlib import Path
from loguru import logger


DIAMOND_BIN = Path(sys.prefix) / "bin" / "diamond"


def call_diamond_makedb(fasta: Path, tmpdir: Path, threads: int, kwargs: list[str]) -> Path:
    """Create a diamond database file (.dmnd) from a sample fasta file.
    """
    # set sample fasta_db to diamond suffix
    assert DIAMOND_BIN.exists(), f"diamond binary not found at {DIAMOND_BIN}"
    assert fasta.exists(), f"query sequence not found: {fasta}"

    # make database file for this sample (in tmp folder?)
    db_path = tmpdir / (fasta.name + ".db")
    cmd = [
        str(DIAMOND_BIN), "makedb",
        "--ignore-warnings",
        "--in", str(fasta),
        "--db", str(db_path),
        "--threads", str(threads),
    ]

    # additional kwargs of any kind
    if kwargs:
        for argstr in kwargs:
            cmd += [argstr]

    # log command and run
    logger.debug(" ".join(cmd))
    with Popen(cmd, stderr=STDOUT, stdout=PIPE) as proc:
        out = proc.communicate()
    if proc.returncode:
        out = out[0].decode()
        raise ValueError(out)
    logger.complete()
    return db_path


def call_diamond_blastp(
    query: Path,
    database: Path,
    threads: int,
    evalue: float,
    bitscore: float | None,
    max_target_seqs: int,
    compress: bool,
    tmpdir: Path,
    ultra_sensitive: bool,
    outfmt: list[str],
    kwargs: list[str],
) -> str:
    """Create a diamond blastp result file for two samples.

    Notes
    -----
    - repeat masking is applied to query and DB unless --masking 0
    - memory limits can be reduced by lowering the block size -b. (ultasensitive mode benefits little from increasing this parameter)
    - diamond is optimized for millions of proteins.

    The following ordered data columns are contained in the blastp TSV:
    0. qseqid = query id
    1. sseqid = subject id
    2. pident = percentage of identical matches
    3. length = alignment length
    4. mismatch = number of mismatches
    5. gapopen = number of gap openings
    6. qstart = start of alignment in query
    7. qend = end of alignment in query
    8. sstart = start of alignment in subject
    9. send = end of alignment in subject
    10. evalue = Expect value
    11. bitscore = Bit score
    12. corrected_bitscore = Bit score corrected for edge effects
    13. qlen = length of the query sequence
    14. slen = length of the subject sequence
    """
    # logger.warning(outfmt)
    # standard outfmt 6 + [corrected_bitscore, qlen, slen]
    FIELDS = [
        "qseqid", "sseqid",
        "pident", "length",
        "mismatch", "gapopen",
        "qstart", "qend",
        "sstart", "send",
        "evalue", "bitscore",
        "corrected_bitscore", "qlen", "slen",  # these are extra
    ]

    # make database file for this sample (in tmp folder?)
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "--ignore-warnings",
        "--quiet",
        "--compress", str(int(compress)),
        "--query", str(query),
        "--db", str(database),
        "--threads", str(threads),
        "--max-target-seqs", str(max_target_seqs),
    ]

    # set outfmt. default is outfmt6 + [corrected_bitscore, qlen, slen]
    if outfmt is None:
        cmd += ["--outfmt", "6"] + FIELDS
    elif not outfmt:
        cmd += ["--outfmt", "6"]
    else:
        cmd += ["--outfmt", "6"] + outfmt

    # score cutoff options
    if bitscore is not None:
        cmd += ["--min-score", str(bitscore)]
    else:
        cmd += ["--evalue", str(evalue)]

    # add sensitivity arg
    if ultra_sensitive:
        cmd += ["--ultra-sensitive"]

    # additional kwargs of any kind
    if kwargs:
        for argstr in kwargs:
            if ' ' in argstr:
                cmd += argstr.split()
            else:
                cmd += [argstr]

    # TODO: optionally add header
    # if header:
        # pass

    # log command, run, and catch errors
    logger.debug(" ".join(cmd))
    # raise SystemExit(0)
    block_size = 1024 * 1024 * 2  # 2Mb
    process = Popen(cmd, stderr=PIPE, stdout=PIPE)
    while True:
        block = process.stdout.read(block_size)
        if not block:
            break
        sys.stdout.buffer.write(block)
        sys.stdout.flush()

    # capture errors
    _, err = process.communicate()
    process.stdout.close()
    if process.returncode != 0:
        logger.error(err.decode())
    process.wait()
    logger.complete()


def run_diamond_bl(args):
    """Runs diamond blast with args parsed from sys or CLI."""
    # make database if target is a fasta
    if args.database or args.target.suffix == ".db":
        dbfile = Path(args.target)
    else:
        dbfile = call_diamond_makedb(args.target, args.tmpdir, args.threads, args.kwargs_makedb)

    # search query against db > stdout
    call_diamond_blastp(args.query, dbfile, args.threads, args.evalue, args.bitscore, args.max_target_seqs, args.compress, args.tmpdir, args.ultra_sensitive, args.outfmt, args.kwargs_blastp)

    # cleanup database file if it was created by twig, and exists, and not asked to be kept
    if dbfile.exists() and (not args.save_db) and (not args.database):
        dbfile.unlink()

def main():
    from ..cli.subcommands import get_parser_diamond_bl
    parser = get_parser_diamond_bl()
    args = parser.parse_args()
    run_diamond_bl(args)


if __name__ == "__main__":
    main()
