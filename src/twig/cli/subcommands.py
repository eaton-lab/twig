#!/usr/bin/env python

import re
from textwrap import dedent
from pathlib import Path
from tempfile import gettempdir
from argparse import ArgumentParser, RawDescriptionHelpFormatter

TMPDIR = gettempdir()
ISOFORM_REGEX_DEFAULT = r"^([^|]+)\|.*?__(.+?)_i\d+"
# group 1 (shared part up until |)
# group 2 (after __ and up until first _i)



def get_parser_csubst(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
    prog="csubst",
    usage="twig csubst -i MSA -t TREE -g FG [options]",
    help="run csubst in wrapper",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=dedent("""
        -------------------------------------------------------------------
        |
        -------------------------------------------------------------------
        |
        |
        |
        -------------------------------------------------------------------
    """),
    epilog=dedent("""
        Examples
        --------
        $ twig csubst -a MSA -t NWK -g FG -o OUT/PRE   # OUT/PRE.cb...
        $ ...
    """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-a", "--alignment", type=Path, metavar="path", required=True, help="input CDS alignment")
    parser.add_argument("-t", "--tree", type=Path, metavar="path", required=True, help="input rooted tree file")
    parser.add_argument("-g", "--foreground", type=Path, metavar="path", help="foreground file")
    # parser.add_argument("-o", "--out", type=Path, metavar="path", help="out prefix; default is input path [{input}]")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", default="./csubst", help="Out directory will be created at this path [%(default)s]")
    # parser.add_argument("-p", "--prefix", type=str, metavar="str", help="optional outfile prefix. If None the cds filename is used")

    parser.add_argument("-m", "--max-arity", type=int, metavar="int", default=2, help="max combinatorial number of branches (K) [%(default)s]")
    parser.add_argument("-u", "--exhaustive-until", type=int, metavar="int", default=1, help="perform exhaustive (non-heuristic) search up N branch combs [%(default)s]")
    parser.add_argument("-c", "--cutoff-stat", type=str, metavar="str", default="OCNany2spe,2.0|omegaCany2spe,5.0", help="Cutoff stats for searching higher-order branch combs [%(default)s]")

    parser.add_argument("-b", "--foreground-table", action="store_true", help="foreground file is a table (fg_format=2)")
    parser.add_argument("-w", "--fg-exclude-wg", action="store_true", help="sets 'yes' to exclude within-group comparisons")
    parser.add_argument("-s", "--fg-stem-only", action="store_true", help="sets 'yes' to only stem comparisons")

    # others
    parser.add_argument("-e", "--env", type=Path, metavar="path", default="csubst", help="conda env name where 'csubst' in installed [csubst]")
    parser.add_argument("-j", "--threads", type=int, metavar="int", default=1, help="number of threads")
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_diamond_bl(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="diamond-bl",
        usage="diamond-bl [options]",
        help="search blastp hits of one protein fasta to another",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
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
        epilog=dedent("""
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


def get_parser_diamond_pw(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="diamond-pw",
        usage="%(prog)s paths [args]",
        help="write .tsv blast hits for all pairs of fasta sequence files",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | diamond-pw: write .tsv's of all-by-all pairwise diamond blastp  |
            -------------------------------------------------------------------
            | Perform all-by-all blastp on a set of fasta protein sequence
            | files. There are many more options available from the `diamond
            | blast` tool itself. This is intended only as a simple wrapper.  |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 --no-self
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 -k 1 -j 10 -t 4

            # Example result dir (blast/) contains
            # d1_d1.tsv  d2_d2.tsv  d3_d3.tsv
            # d1_d2.tsv  d1_d3.tsv  d2_d3.tsv
            # d2_d1.tsv  d3_d1.tsv  d3_d2.tsv
            # sample_map.tsv
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("-i", "--input", type=Path, metavar="path", help="fasta sequence files for all-by-all blastp")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", default="./blast", help="directory to write results (created if necessary) [%(default)s]")
    parser.add_argument("-e", "--evalue", type=float, metavar="float", default=1e-5, help="max evalue to report an alignment [%(default)s]")
    parser.add_argument("-m", "--min-bitscore", type=float, metavar="float", default=None, help="min bit-score to report an alignment (overrides -e) [%(default)s]")
    parser.add_argument("-k", "--max-target-seqs", type=int, metavar="int", default=25, help="max n target seqs to report hits for [%(default)s]")
    parser.add_argument("-d", "--dbdir", type=Path, metavar="path", default=TMPDIR, help="directory for temp .db files [%(default)s]")
    parser.add_argument("-t", "--threads", type=int, default=4, help="number of threads per diamond job [%(default)s]")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="number of diamond jobs to run in parallel [%(default)s]")
    parser.add_argument("-n", "--no-self", action="store_true", help="do not perform self-self search")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite .db files if they exist in dbdir")
    # parser.add_argument("-r", "--relabel-headers", type=int, default=1, help="...")
    # parser.add_argument("-r", "--relabel-taxa", type=int, default=1, help="...")
    # parser.add_argument("-I", "--imap", type=Path, help="map file to translate sample labels.")
    return parser


def get_parser_format_fasta(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for format-fasta tool.
    """
    KWARGS = dict(
        prog="format-fasta",
        usage="twig format-fasta FASTA [options]",
        help="format fasta to sort/filter/subset/relabel/modify sequences",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | format-fasta: sort/filter/subset/relabel/modify headers or seqs |
            -------------------------------------------------------------------
            | This tool is used to modify a fasta file by performing one or   |
            | more of the following operations in order: (1) sorting; (2)     |
            | filtering; (3) sub-selecting; (4) relabeling; and (5) modifying |
            | headers and sequences. A modified fasta is written to STDOUT.   |
            | A subset of sequences/scaffolds can be returned in any order    |
            | with --subset-idx or --subset-names but note that --sort-alpha  |
            | or --sort-len override the ordering. This tool can be useful    |
            | for modifying genome fasta files or CDS or PEP sequence files.  |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # operations on genome scaffolds
            $ twig format-fasta GENOME --subset 3 --relabel-list A B C > genome.fa
            $ twig format-fasta GENOME --subset-idx {0..3} --relabel-prefix-idx A_ > subgenomeA.fa
            $ twig format-fasta GENOME --subset-names Chr1 Chr2 --relabel-list chromosome_1 chromosome_2 > genome.fa
            $ twig format-fasta GENOME --sort-len --subset 8 --max-width 80 > chromosomes.fa

            # operations on sequences (e.g., cds, pep)
            $ twig format-fasta FA --relabel-prefix-idx spp --map NAME.map.tsv > NAME.fa
            $ twig format-fasta FA --sort-len --min-len 200 --subset 50 > NAME.fa
            $ twig format-fasta FA --max-width 80 --unmask > NAME.fa

            # complex operations: select sequences from a subset of scaffolds
            $ ids=$(twig format-gff GFF --subset-names Chr1 --feat mRNA --export-ids)
            $ twig format-fasta FA --subset-names $ids > Chr1.pep.fa
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # add arguments
    parser.add_argument("fasta", type=Path, help="a fasta sequence file (can be .gz)")

    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort chromosomes by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort chromosomes by length (longest to shortest)")

    relabel_options = parser.add_mutually_exclusive_group()
    relabel_options.add_argument("--relabel-prefix", type=str, metavar="str", default="", help="append prefix label to start of existing labels [%(default)s]")
    relabel_options.add_argument("--relabel-prefix-idx", type=str, metavar="str", default="", help="overwrite labels with {{prefix}}{{counter}} using incrementing counter [%(default)s]")
    relabel_options.add_argument("--relabel-list", type=str, metavar="str", nargs="+", help="overwrite labels with new labels provided as a list [%(default)s]")
    parser.add_argument("--map", type=Path, metavar="Path", default=None, help="path to write optional relabel translation map (TSV format) [%(default)s]")

    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N sequences after optional sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more sequences by 1-based index after optional sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more sequences by name before optional relabeling")

    parser.add_argument("--min-len", type=int, metavar="int", default=None, help="sequences shorter than this length are excluded (default=None)")
    parser.add_argument("--max-width", type=int, metavar="int", default=None, help="write sequences with max width textwrapping (e.g., 80)")
    parser.add_argument("--unmask", action="store_true", help="convert any soft-masked sequences to upper case")
    parser.add_argument("--revcomp", action="store_true", help="reverse complement DNA sequences, or reverse AA sequences")

    # compress?
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_format_gff(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for format-gff tool.
    """
    KWARGS = dict(
        prog="format-gff",
        # usage="%(prog)s fasta [options]",
        usage="twig format-gff GFF [options]",
        help="format gff to filter and relabel genes and scaffolds",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | format-gff: filter/relabel genome annotation table              |
            -------------------------------------------------------------------
            | This tool is used to ...
            -------------------------------------------------------------------
        """),
        epilog=dedent(r"""
            Examples
            --------
            # select a subset of scaffolds by name
            $ format-gff GFF --subset-names Chr1 > Chr1.gff

            # sort scaffolds by length and select first 10
            $ format-gff GFF --sort-len --subset 10 > Chr1-10.gff

            # relabel sorted scaffold names using str subst, save map of old names to new
            $ format-gff GFF --sort-len --relabel {idx}_Chr{sidx}.{gidx} --relabel-map scaff-map.tsv > relabeled.gff

            # filter to include only gene features
            $ format-gff GFF --features gene > genes.gff

            # subset to genes on Chr1, keep only ID attrs, and relabel IDs with spp prefix
            $ format-gff GFF --subset-names Chr1 --feat gene --attrs ID --relabel-ids spp_{id} > relabeled.gff

            # export GFF data to bed format and remove overlapping/duplicate features
            $ format-gff GFF --subset-names Chr1 --feat mRNA --to bed --filter-duplicates > Chr1-transcripts.bed

            # get all sequence IDs of mRNAs on Chr1
            $ format-gff GFF --subset-names Chr1 --feat mRNA --attrs ID --to ids > Chr1-ids.txt
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # add arguments
    parser.add_argument("GFF", type=Path, help="a gff/gtf annotation table (can be .gz)")

    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort scaffolds by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort scaffolds by length (longest to shortest)")

    # relabel scaffolds. If you relabel in GFF you need to relabel same way in FA file using format-fasta --relabel-map MAP
    relabel_options = parser.add_mutually_exclusive_group()
    relabel_options.add_argument("--relabel", type=str, metavar="str", default="", help=r"relabel scaffs using str subst w/ {scaff},{sidx}. Applies after sorting/subset")
    relabel_options.add_argument("--relabel-names", type=str, metavar="str", default="", help=r"relabel scaffs with list of new names. Applies after sorting/subset")
    relabel_options.add_argument("--relabel-map", type=str, metavar="str", default="", help=r"optional path to write tsv mapping old scaff names to new")

    # select a subset of scaffolds
    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N scaffolds. Applies after sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more scaffolds by 1-based index. Applies after sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more scaffolds by name. Applies before relabeling/sorting")

    parser.add_argument("--features", type=str, metavar="str", nargs="+", help="subselect one or more features (e.g., gene) to keep (default=all)")
    parser.add_argument("--attrs", type=str, metavar="str", nargs="+", help="subselect one or more attributes (e.g., ID) to keep (default=all)")

    # relabel IDs. If you relabel in GFF you need to relabel same way in PEP, CDS files using format-fasta --relabel-map MAP
    parser.add_argument("--relabel-ids", type=str, metavar="str", default="", help=r"relabel attr IDs using str subst w/ {scaff},{sidx},{gidx},{idx},{type}")
    parser.add_argument("--relabel-ids-map", type=str, metavar="str", default="", help=r"optional path to write tsv mapping old IDS to new")

    parser.add_argument("--exclude-comments", action="store_true", help="exclude comments from the output")
    parser.add_argument("--add-to-comments", action="store_true", help="add this cmd to comments section")
    parser.add_argument("--filter-duplicates", action="store_true", help="keep only the highest score, or first, of duplicates coordinate features")

    parser.add_argument("--to", type=str, metavar="{gff, bed, ids}", default="gff", help="output format (default=gff)")
    # export_options.add_argument("--to-gff2", action="store_true", help="write BED format selected features")
    # export_options.add_argument("--to-gff3", action="store_true", help="write BED format selected features")
    # export_options.add_argument("--to-gtf", action="store_true", help="write BED format selected features")

    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_genome_table(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="genome-table",
        usage="twig genome-table FA GFF [options]",
        help="write table with chrom/scaff [names lengths ngenes]",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | genome-table: get table of chrom/scaff [names, lengths, ngenes] |
            -------------------------------------------------------------------
            | Parses the genome fasta (and annotation) files to extract       |
            | scaffold names, lengths, and gene numbers and write as a tsv    |
            | to stdout. This table is useful for downstream analysis such as |
            | dotplots and other synteny analyses. If an annotation file is   |
            | not provided then only names and lengths are extracted. Options |
            | can be used to sort and/or subselect scaffolds for inclusion.   |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ genome-table FA > genome.info.tsv
            $ genome-table FA GFF > genome.info.tsv
            $ genome-table FA GFF --sort-len --subset 8 > genome.info.tsv
            $ genome-table FA GFF --subset-names chr1 chr2 > genome.info.tsv
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("genome", type=Path, help="fasta genome sequence (can be .gz)")
    parser.add_argument("annotation", type=Path, nargs="?", default=None, help="gff annotation file (optional; can be .gz)")
    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort chromosomes by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort chromosomes by length")
    # parser.add_argument("--subset", type=int, metavar="int", default=None, help="write only the first N scaffolds after sorting")
    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N scaffolds after optional sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more scaffolds by 1-based index after optional sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more scaffolds by name")

    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_macse_prep(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="macse-prep",
        usage="macse-prep -i CDS -o OUTDIR [options]",
        help="prepare CDS for macse alignment (trim,filter,iso-collapse)",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | macse-prep: prepare CDS for alignment (trim, filter, iso-collapse)
            -------------------------------------------------------------------
            | Macse ...
            | This implements `macse -prog trimNonHomologousFragments` and
            | additional steps to ...
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig macse-prep -i CDS -o OUT.nt.fa
            $ twig macse-prep -i CDS -mh 0.1 -mi 0.5 -ti 50 -te 50 -mc 15 -c 20
            $ twig macse-prep -i CDS -mh 0.3 -mi 0.8 -ti 25 -te 25 -mc 15 -c 20
            $ twig macse-prep -i CDS -mh 0.5 -mc 10 -k -xa -ml 200 -e '^sppA.*'
            $ twig macse-prep -i CDS -s '^sppA.*'

            # run parallel jobs on many cds files
            $ parallel -j 10 "twig macse-prep -i {} ..."  ::: CDS/*.fa

            # full pipeline
            $ twig macse-prep -i CDS -o TRIM     # {TRIM}
            $ twig macse-align -i TRIM -o MSA    # {MSA}
            $ twig macse-refine -i MSA -o MSA    # {MSA}.nt.fa, {MSA}.aa.fa
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS (aligned or unaligned)")
    parser.add_argument("-o", "--outpath", type=Path, metavar="path", required=True, help="path to write nt fasta result")
    parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
    parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
    # options
    parser.add_argument("-ml", "--min-length", type=int, metavar="int", default=0, help="min nt sequence length after trimming [%(default)s]")
    parser.add_argument("-hf", "--min-homology-full", type=float, metavar="float", default=0.1, help="min homology required w/ >=mc others across full sequence [%(default)s]")
    parser.add_argument("-hi", "--min-homology-internal", type=float, metavar="float", default=0.5, help="min homology required w/ >=mc others in the internal sequence [%(default)s]")
    parser.add_argument("-hc", "--min-homology-coverage", type=int, metavar="int", default=3, help="min samples a seq must share homology with at >= mh and mi [%(default)s]")
    parser.add_argument("-ti", "--min-trim-length-homology-internal", type=int, metavar="int", default=50, help="trim fragments w/ len <ti and homology <mi w/ <mc sequences [%(default)s]")
    parser.add_argument("-te", "--min-trim-length-homology-external", type=int, metavar="int", default=50, help="trim fragments w/ len <tx and homology <mh w/ <mc sequences [%(default)s]")
    parser.add_argument("-mm", "--mem-length", type=int, metavar="int", default=6, help="homology is the prop of aa Maximum Exact Matches of this length [%(default)s]")
    parser.add_argument("-is", "--isoform-regex", type=re.compile, metavar="str", default=ISOFORM_REGEX_DEFAULT, help="regex used to group isoform sequences ['%(default)s']")
    # parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
    # parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment using a sliding 'half_window_size' defined as ... [%(default)s]")
    # parser.add_argument("-xa", "--skip-alignment", action="store_true", help="skip alignment step and only trim/filter sequences")

    # choose one or more
    # parser.add_argument("-xt", "--skip-trim-and-filter", action="store_true", help="skip trim and filter step")
    parser.add_argument("-xi", "--skip-isoform-collapse", action="store_true", help="skip isoform collapse step")

    # others
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_macse_align(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="macse-align",
        usage="macse-align -i CDS -o OUTDIR [options]",
        help="run macse alignment (CDS/AA)",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | macse-align: ...
            -------------------------------------------------------------------
            | Macse ...
            | This implements `macse -prog AlignSequences` and
            | additional steps to ...
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig macse-align -i CDS -o /OUT/CDS.msa

            # run parallel jobs on many cds files
            $ parallel -j 10 "twig macse-align -i {} ..." ::: CDS/*.nt.fa

            # full pipeline
            $ twig macse-prep -i CDS -o TRIM     # {TRIM}
            $ twig macse-align -i TRIM -o MSA    # {MSA}
            $ twig macse-refine -i MSA -o MSA    # {MSA}.nt.fa, {MSA}.aa.fa
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS (aligned or unaligned)")
    parser.add_argument("-o", "--outpath", type=Path, metavar="path", help="path to write aligned nt fasta")
    # others
    parser.add_argument("-m", "--max-refine-iter", type=int, metavar="int", default=-1, help="max refinement iterations during optimizing [default -1 = no limit]")
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    # parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_macse_refine(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="macse-refine",
        usage="macse-refine -i CDS [options]",
        help="...",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | macse-refine: CDS/AA jointly and filter low homology seqs
            -------------------------------------------------------------------
            | Macse codon-aware alignment of CDS and AA sequences. This method
            | also includes options to filter and trim sequences to remove low
            | homology fragments or sequences; sequences that are too short;
            | and extra isoforms. If a result with -p PREFIX name exists in -o
            | OUTDIR it will be skipped, allowing for easy checkpointing to
            | run in parallel on a large set of sequences and restart if
            | interrupted.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig macse-refine -i CDS -o PATH/PRE
            $ twig macse-refine -i CDS -mh 0.1 -mi 0.5 -ti 66 -te 66 -mc 15 -c 20
            $ twig macse-refine -i CDS -mh 0.3 -mi 0.8 -ti 33 -te 33 -mc 15 -c 20
            $ twig macse-refine -i CDS -mh 0.5 -mc 10 -k -xa -ml 200 -e '^sppA.*'
            $ twig macse-refine -i CDS -s '^sppA.*'

            # run parallel jobs on many cds files
            $ parallel -j 10 "twig macse-refine -i {} ..."  ::: CDS/*.msa.nt.fa

            # full pipeline
            $ twig macse-prep -i CDS                     # {CDS}.nt.fa, ...
            $ twig macse-align -i CDS.nt.fa -o CDS       # {CDS}.msa.nt.fa, ...
            $ twig macse-refine -i CDS.msa.nt.fa -o CDS  # {CDS}.msa.refined.nt.fa, ...
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input aligned CDS")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="out prefix; default is input path [{input}]")
    parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
    parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
    parser.add_argument("-t", "--tree", type=Path, metavar="path", help="optional newick file to subsample genes present in tree")
    # options
    parser.add_argument("-ml", "--min-length", type=int, metavar="int", default=0, help="min length of non-missing sequence in a sample [%(default)s]")
    parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
    parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment using a sliding 'half_window_size' defined as ... [%(default)s]")
    parser.add_argument("-if", "--codon-int-fs", type=str, metavar="str", default="NNN", help="codon to sub for internal frame shift [NNN]")
    parser.add_argument("-ef", "--codon-ext-fs", type=str, metavar="str", default="NNN", help="codon to sub for external frame shift [NNN]")
    parser.add_argument("-fs", "--codon-final-stop", type=str, metavar="str", default="NNN", help="codon to sub for final stop [NNN]")
    parser.add_argument("-is", "--codon-int-stop", type=str, metavar="str", default="NNN", help="codon to sub for internal stop [NNN]")

    parser.add_argument("-r", "--refine-alignment", action="store_true", help="refine alignment")
    parser.add_argument("-R", "--refine-alignment-if", action="store_true", help="refine alignment only if >=1 sequences are filtered out")
    parser.add_argument("-ri", "--max-iter-refine-alignment", type=int, metavar="int", default=-1, help="max iterations in refine alignment [%(default)s]")

    # others
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_tree_filter(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="tree-filter",
        usage="tree-filter [options]",
        help="filter and relabel a set of trees for downstream analyses",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
        description=dedent("""
            -------------------------------------------------------------------
            | tree-filter: Filter/Relabel a set of trees for downstream analyses
            -------------------------------------------------------------------
            | The order of [optional] operations is (1) parse names from delim
            | split tip labels; (2) subsample and/or relabel names by imap;
            | (3) exclude outlier edges; (4) collapse outgroups; (5) exclude
            | trees w/ < min-tips; (6) exclude trees w/ > min-copies in any name.
            | Note that if names are parsed by delim, then the imap should contain
            | the parsed names. The labels on output trees will be the same as in
            | the inputs unless --relabel-delim or --relabel-imap is used. But
            | parsing shorter names can be useful to match names to imap even
            | if keeping full names. Filtering stats are written to stderr, which
            | can be redirected to a log file with `2> log.txt`.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # relabel by split-select-join on tip names
            $ twig tree-filter -i NWK -d '-' -di 0 2 -dj '-rd' > relabeled-trees.nwk

            # get only single-copy gene trees
            $ twig tree-filter -i NWK -c 1 > single-copy-trees.nwk

            # get only trees w/ >20 tips
            $ twig tree-filter -i NWK -m 20 > min20-trees.nwk

            # exclude outlier edges
            $ twig tree-filter -i NWK --exclude -ei 4 > cleaned-trees.nwk

            Examples using IMAP/MINMAP
            --------------------------
            # subsample to names in imap
            $ twig tree-filter -i NWK -I IMAP > subsample-trees.nwk

            # subsample and relabel to pop names in imap
            $ twig tree-filter -i NWK -I IMAP -ri > relabeled-trees.nwk

            # subsample, relabel, and keep only one outgroup (best to root before)
            $ twig tree-filter -i NWK -I IMAP -ri --collapse > final-trees.nwk

            Example on a many trees in different paths
            -------------------------------------------
            $ parallel "twig tree-filter -i {} ... > {}.filtered" ::: TREES/*.nwk

            Example IMAP                       Example MINMAP
            --------------------               ----------------
            sampleA     pop1                   pop1       1
            sampleB     pop1                   pop2       1
            sampleC     pop2                   outgroup   1
            sampleD     pop2
            sampleE     outgroup
            sampleF     outgroup
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # filter on names
    # parser.add_argument("-i", "--include", type=str, metavar="str", nargs="+", help="subset of names to keep")
    parser.add_argument("-I", "--imap", type=Path, metavar="path", help=r"filepath listing names to keep, or assigning (name population)")
    parser.add_argument("-M", "--minmap", type=Path, metavar="path", help=r"filepath listing (population mincov) for filters")

    # filter on tree data
    parser.add_argument("-m", "--min-tips", type=int, metavar="int", default=4, help="min tips after pruning [4]")
    parser.add_argument("-c", "--max-copies", type=int, metavar="int", help="filter trees with >c gene copies from a taxon [None]")
    parser.add_argument("-eo", "--edge-outlier-outgroup", type=float, metavar="float", default=10, help="exclude 'outgroup' population edges if >eo stdev from mean [10]")
    parser.add_argument("-ei", "--edge-outlier-ingroup", type=float, metavar="float", default=5, help="exclude non 'outgroup' population edges if >ei stdev from mean [5]")

    # actions
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    parser.add_argument("-ri", "--relabel-imap", action="store_true", help="relabel tips to their imap mapped names")
    # parser.add_argument("-x", "--nexus", action="store_true", help="export in NEXUS format. Retains path/tree names")
    parser.add_argument("--subsample", action="store_true", help="subsample to include only tips in imap")
    parser.add_argument("--exclude-outliers", action="store_true", help="exclude tips with outlier edge lengths (>ei or >eo)")
    parser.add_argument("--require-outgroups", action="store_true", help="require at least one 'outgroup' sample")
    parser.add_argument("--collapse-outgroups", action="store_true", help="keep only the most distant 'outgroup' (assumes rooted trees)")

    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_tree_rooter(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="tree-rooter",
        usage="tree-rooter [options]",
        help="Re-root gene trees by outgroup, species tree, and/or MAD",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
        description=dedent("""
            -------------------------------------------------------------------
            | tree-rooter: ...
            -------------------------------------------------------------------
            | Root trees by designating an outgroup, or a set of potential
            | outgroups in the case that a set of trees may or may not include
            | a single outgroup. The ...
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # root trees on a specific outgroup
            $ twig tree-rooter -t NWK -r A > rooted-trees.nwk

            # root trees on an ordered set of options of outgroup clades
            $ twig tree-rooter -t NWK -r A B A,B > rooted-trees.nwk

            # root trees using minimal ancestor deviation (MAD)
            $ twig tree-rooter -t NWK --mad > rooted-trees.nwk

            # root trees on outgroups if present and mad root the rest
            $ twig tree-rooter -t NWK -r A B --mad > rooted-trees.nwk

            # root trees on outgroups in order of a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE > rooted-trees.nwk

            # root trees on a restricted set of outgroups from a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE -r Z Y > rooted-trees.nwk

            # root trees on outgroups defined in IMAP given a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE -I IMAP > rooted-trees.nwk

            # return trees that could not be rooted given the options
            $ twig tree-rooter -t NWK -s SPTREE -x > unrooted-trees.nwk

            # parallel process many files individually
            $ parallel -j 10 "twig tree-rooter -i {} -s SPTREE > ROOTED/{}" ::: TREES/*.nwk
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # rooting options
    parser.add_argument("-s", "--sptree", type=Path, metavar="path", help="rooted species tree")
    parser.add_argument("-r", "--outgroups", type=str, metavar="str", nargs="+", help="list outgroup tip labels")
    parser.add_argument("-I", "--imap", type=Path, metavar="path", help="get outgroup assignment from tabular imap file listing 'sample\tpopulation'")

    # return option
    parser.add_argument("-m", "--mad", action="store_true", help="use minimal ancestor deviation to estimate root of remaining unrooted trees")
    parser.add_argument("-x", "--not-rooted", action="store_true", help="return the trees that could not be rooted given the options")
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    # logging
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="CRITICAL", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_partition_cds(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="partition-cds",
        usage="partition-cds [options]",
        help="write a partition file for phylogenetic inference",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | partition-cds: write a partition file for phylogenetic analysis |
            -------------------------------------------------------------------
            | Write a partition file for a CDS alignment to specify substitution
            | models for phylogenetic inference based on codon position and the
            | length of the CDS. Writes to stdout.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig partition-cds -i ALN.fa -n 2 -m GTR+G GTR+G > ALN.partition
            $ twig partition-cds -i ALN.fa -n 3 -m GTR+G GTR+G GTR+G > ALN.partition
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="fasta CDS sequence")
    parser.add_argument("-n", "--number", choices=[2, 3], type=int, metavar="int", default=3, help="number of partitions (2 or 3) [3]")
    parser.add_argument("-m", "--model", type=str, metavar="str", nargs="+", default=["GTR+G"], help="subst model for each position (1+2,3) or (1,2,3) for 2 or 3 partitions")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser




# ...

