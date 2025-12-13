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
import textwrap
from subprocess import Popen, STDOUT, PIPE
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
from tempfile import gettempdir
from loguru import logger

DIAMOND_BIN = Path(sys.prefix) / "bin" / "diamond"
TMPDIR = gettempdir()


KWARGS = dict(
    prog="diamond-bl",
    usage="diamond-bl [options]",
    help="search blastp hits of one protein fasta to another",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | diamond-bl: write .tsv of diamond blastp hits to stdout         |
        -------------------------------------------------------------------
        | There are many more options available in the `diamond blastp`   |
        | tool itself. This is intended only as a simple wrapper to make  |
        | common workflows easily available in twig. This can create an   |
        | indexed database, search against a database, or create and      |
        | search against a database in one command. Results are written to|
        | a TSV file with default outfmt-6 + corrected_bitscore qlen slen.|
        | Use --outfmt to set to standard outfmt-6, or to list custom     |
        | fields to append. This checks for diamond installation in PATH  |
        | or -B path, and raises an error if not found.                   |
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        # find hits of A to B
        $ diamond-bl -q A.fa -d B.fa > A_B.tsv

        # find hits of A to B and save and reuse B database
        $ diamond-bl -q A.fa -d B.fa -s > A_B.tsv
        $ diamond-bl -q C.fa -D B.db > C_B.tsv

        # set options for parallel, tmp dir, and sensitivity
        $ diamond-bl -q A.fa -d B.fa -j 10 -d /tmp -e 1e-5 > A_B.tsv

        # set to standard outfmt-6
        $ diamond-bl -q A.fa -d B.fa --outfmt > A_B.tsv

        # set to custom outfmt (metrics to append to outfmt-6)
        $ diamond-bl -q A.fa -d B.fa --outfmt qseqid sseqid evalue > A_B.tsv

        # set any additional arbitrary kwargs allowed by diamond blast
        $ diamond-bl -q A.fa -D B.db --kwargs-blastp '--taxonlist 33090' ... > A_B.tsv

        # set path to diamond binary
        $ diamond-bl -q A.fa -t B.fa -B /opt/diamond > A_B.tsv
    """)
)

def get_parser_diamond_blastp(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # add arguments
    parser.add_argument("-q", "--query", type=Path, metavar="path", required=True, help="fasta protein sequence to align to the target")
    parser.add_argument("-t", "--target", type=Path, metavar="path", required=True, help="fasta protein sequence or diamond db file as target")
    parser.add_argument("-e", "--evalue", type=float, metavar="float", default=0.001, help="max evalue to report an alignment; default=0.001")
    parser.add_argument("-b", "--bitscore", type=float, metavar="float", default=None, help="min bit-score to report an alignment (overrides -e); default=None")
    parser.add_argument("-k", "--max-target-seqs", type=int, metavar="int", default=25, help="max number of target sequences to report alignments for; default=25")
    parser.add_argument("-j", "--threads", type=int, metavar="int", default=4, help="number of threads per diamond job; default=4")
    parser.add_argument("-d", "--tmpdir", type=Path, metavar="path", default=TMPDIR, help=f"directory for tmp files; default={TMPDIR}")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite .db files if they exist in tmpdir")
    parser.add_argument("-c", "--compress", action="store_true", help="gzip compress output")
    parser.add_argument("-u", "--ultra-sensitive", action="store_true", help="enable ultra sensitive mode")
    parser.add_argument("-D", "--database", action="store_true", help="specify that target file is already a diamond indexed database")
    parser.add_argument("-s", "--save-db", action="store_true", help="do not remove diamond indexed .db in tmpdir, if created")
    parser.add_argument("-H", "--header", action="store_true", help="add header to TSV")
    parser.add_argument("-B", "--binary", type=Path, metavar="path", help="designate path to diamond binary if not in $PATH")
    parser.add_argument("--outfmt", type=str, metavar="str", nargs="*", default=None, help="set to standard outfmt-6 (no additional args) or set custom fields")
    # parser.add_argument("--blastx", action="store_true", help="search DNA against a protein database")
    # parser.add_argument("-m", "--map", type=Path, help="map file to translate sample labels.")
    parser.add_argument("--kwargs-makedb", type=str, nargs="*", default=None, help="additional makedb args")
    parser.add_argument("--kwargs-blastp", type=str, nargs="*", default=None, help="additional blastp args")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")

    return parser


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


def run_diamond_blastp(args):
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
    parser = get_parser_diamond_blastp()
    args = parser.parse_args()
    run_diamond_blastp(args)


if __name__ == "__main__":
    main()
