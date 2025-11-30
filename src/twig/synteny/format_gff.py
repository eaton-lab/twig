#!/usr/bin/env python

"""Utilities for working with GFF files.

Example
-------
$ python format_gff.py GFF --relabel-prefix --feat gene 
$ python format_gff.py GFF --subset-names Chr01 --feat mRNA --export-bed 

TODO
-----
- allow i/o gzip
- relabel scaffolds with rich text
- relabel ids with rich text
- update parent ids option
- filter overlapping features?

rich text options: 
  {scaff} : scaffold name
  {sidx}  : scaffold index  
  {gene}  : gene name
  {type}  : feature type
  {ctype} : feature type only of children
  {gidx}  : gene index in scaffold
  {idx}   : index across all scaffolds
  {cidx}  : child type count

g{idx}.t{cidx}.{ctype}{cidx}
  gene = g1
  mRNA = g1.t1


WGDI format (wowww?)
------------
A1  ad200s1g03486   106893663   106894130   -   3486    aradu.V14167.gnm1.ann1.Aradu.1NR1F.1
A1  ad200s1g03487   106901310   106905125   +   3487    aradu.V14167.gnm1.ann1.Aradu.P2M19.1
A1  ad200s1g03488   106904773   106908513   -   3488    aradu.V14167.gnm1.ann1.Aradu.W8AID.1
A1  ad200s1g03489   106970538   106976241   -   3489    aradu.V14167.gnm1.ann1.Aradu.B5T68.1
A1  ad200s1g03490   106981155   106983777   -   3490    aradu.V14167.gnm1.ann1.Aradu.F0BNN.1
A1  ad200s1g03491   106984200   106985794   +   3491    aradu.V14167.gnm1.ann1.Aradu.116JI.1
A1  ad200s1g03492   107017590   107021758   +   3492    aradu.V14167.gnm1.ann1.Aradu.NR4MV.1
A1  ad200s1g03493   107021216   107025068   -   3493    aradu.V14167.gnm1.ann1.Aradu.BU5Q7.1
A10 ad200s10g00001  7040    7419    -   1   aradu.V14167.gnm1.ann1.Aradu.V68SN.1
A10 ad200s10g00002  7447    10856   -   2   aradu.V14167.gnm1.ann1.Aradu.6QW61.1
A10 ad200s10g00003  13132   14655   +   3   aradu.V14167.gnm1.ann1.Aradu.65XKD.1
A10 ad200s10g00004  14820   15933   +   4   aradu.V14167.gnm1.ann1.Aradu.KB0UY.1

My format
---------
scaff  source   type    start   end   score   orient  phase  attributes
Chr1     .      gene     100    500     .        +      0     ID={prefix}s1g1
Chr1     .      gene     600    900     .        +      0     ID={prefix}s1g2
...

# what is not subset features to genes only?
Chr1     .      gene     100    500     .        +      0     ID={prefix}s1g1
Chr1     .      mRNA     100    500     .        +      0     ID={prefix}s1g1.{type}{count};Parent={ID}
Chr1     .      gene     600    900     .        +      0     ID={prefix}s1g2
...



"""

import re
import sys
import difflib
import textwrap
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from loguru import logger


KWARGS = dict(
    prog="format-gff",
    # usage="%(prog)s fasta [options]",
    usage="twig format-gff GFF [options]",
    help="format gff to filter and relabel genes and scaffolds",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),    
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | format-gff: filter/relabel genome annotation table              |
        -------------------------------------------------------------------
        | This tool is used to ...
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent(r"""
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


def get_parser_format_gff(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for format-gff tool.
    """
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


def fuzzy_choice_match(value: str, choices: list[str]) -> str:
    """Return a str matching an allowed value as a substring ignoring case
    """
    # convert choices to lowercase for case-insensitive matching
    choices_lower = [choice.lower() for choice in choices]
    value_lower = value.lower()

    # Check if user_input matches the start of any choices (partial matching)
    matches = [choice for choice in choices_lower if choice.startswith(value_lower)]

    # if only 1 possible match return it lower case
    if len(matches) == 1:
        return choices[choices_lower.index(matches[0])]

    # if multiple possible match log warning and exit
    if len(matches) > 1:
        logger.error(f"value '{value}' matches multiple possible options [{matches}], you must disambiguate")
        sys.exit(1)

    # if no partial matches, suggest a close match using difflib
    close_matches = difflib.get_close_matches(value_lower, choices_lower, n=1, cutoff=0.6)
    if close_matches:
        logger.error(f"value '{value}' not recognized. Did you mean {close_matches}?")
    sys.exit(1)


def get_gff_table(gff: Path) -> tuple[pd.DataFrame, str]:
    """Return GFF as a pandas dataframe.
    """
    return pd.read_table(gff, header=None)


def extract_comments(gff: Path, comment_char: str = "#") -> str:
    """Return the comment header lines as a str"""
    # read comment lines until a non-comment line is encountered.
    with open(gff, 'r') as file:
        comments = []
        for line in file:
            if line.startswith(comment_char):
                comments.append(line)
            else:
                break
    return "#".join(comments)


def sort_gff_table(table: pd.DataFrame, sort_alpha: bool, sort_len: bool) -> pd.DataFrame:
    """... """
    # sort table
    if sort_alpha:
        table.sort_values(by=[0, 3, 4], ascending=[True, True, False], inplace=True)
    if sort_len:
        # get length of each scaffold and then sort by list
        chroms = table.iloc[:, 0].unique()
        chroms_to_lens = {i: table.loc[table[0] == i, 4].max() for i in chroms}
        custom_order = sorted(chroms_to_lens, key=lambda x: chroms_to_lens[x], reverse=True)
        cat_type = CategoricalDtype(categories=custom_order, ordered=True)
        table[0] = table[0].astype(cat_type)
        table.sort_values(by=[0, 3, 4], ascending=[True, True, False], inplace=True)
    return table


def subset_gff_table_features(table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    """..."""
    if features:
        table = table[table[2].isin(features)]
    return table


def subset_gff_table_scaffolds(table: pd.DataFrame, subset: int, subset_idx: list[int], subset_names: list[str]) -> pd.DataFrame:
    """..."""
    if subset:
        subset_names = list(table[0].unique())[:subset]
        table = table[table[0].isin(subset_names)]        
    if subset_idx:
        subset_names = [j for (i, j) in enumerate(table[0].unique()) if i + 1 in subset_idx]
        table = table[table[0].isin(subset_names)]        
    if subset_names:
        table = table[table[0].isin(subset_names)]
    return table


def subset_gff_table_attrs(table: pd.DataFrame, attrs: list[str], attr_column: int = 8) -> pd.DataFrame:
    """Return table subselected key=value; pairs from attributes column"""
    def filter_pairs(cell):
        # Split by semicolon to get each key=value pair
        pairs = cell.split(';')
        
        # Filter pairs to keep only those with keys in keys_to_keep
        filtered_pairs = [pair for pair in pairs if pair.split('=')[0] in attrs]
        
        # Join filtered pairs back into a semicolon-separated string
        return ';'.join(filtered_pairs) + ";"
    if attrs:
        table[attr_column] = table[attr_column].apply(filter_pairs)
    return table


def filter_duplicates(table: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicate records with identical start-stop.

    >>> Chr01   AUGUSTUS    gene    4924    12067   .   -   .   ID=g1;
    >>> Chr01   AUGUSTUS    mRNA    4924    12067   1   -   .   ID=g1.t1;Parent=g1;
    >>> Chr01   AUGUSTUS    mRNA    4924    12067   1   -   .   ID=g1.t2;Parent=g1;
    >>> Chr01   AUGUSTUS    CDS 4924    5065    1   -   1   ID=g1.t1.CDS1;Parent=g1.t1;
    >>> Chr01   AUGUSTUS    exon    4924    5065    .   -   .   ID=g1.t1.exon1;Parent=g1.t1;
    >>> Chr01   AUGUSTUS    CDS 4924    5065    1   -   1   ID=g1.t2.CDS1;Parent=g1.t2;
    >>> Chr01   AUGUSTUS    exon    4924    5065    .   -   .   ID=g1.t2.exon1;Parent=g1.t2;
    >>> ...

    >>> Chr01   AUGUSTUS    gene    4924    12067   .   -   .   ID=g1;
    >>> Chr01   AUGUSTUS    mRNA    4924    12067   1   -   .   ID=g1.t1;Parent=g1;
    >>> Chr01   AUGUSTUS    CDS 4924    5065    1   -   1   ID=g1.t1.CDS1;Parent=g1.t1;
    >>> Chr01   AUGUSTUS    exon    4924    5065    .   -   .   ID=g1.t1.exon1;Parent=g1.t1;
    >>> ...
    """
    # to retain row order
    table = table.copy()
    table["_original_index"] = table.index

    # convert score column to numeric, setting errors='coerce' to convert '.' to NaN
    table[5] = pd.to_numeric(table[5], errors='coerce')

    # fill NaN scores (previously '.') with a very low value (e.g., -inf) to deprioritize them
    table[5] = table[5].fillna(-np.inf)

    # filter to select just one or >1 identical start/stop records
    table = table.loc[table.groupby([3, 4, 2])[5].idxmax()]

    # sort by the original index and drop the old index column
    table = table.sort_values("_original_index")
    table = table.drop(columns=["_original_index"])

    # reset -inf scores back to '.'
    table[5] = table[5].replace(-np.inf, '.')
    return table


def relabel_ids_subst(table: pd.DataFrame, text: str, mapfile: Path | None) -> pd.DataFrame:
    """Replace ID strings in attribute column with str substituted labels.

    Note that if IDs are relabeled in the GFF you will likely also want
    to relabel the headers in the PEP and CDS files. Those files do not
    have the meta information of scaff names, so the new names here
    should be written out as a map file and used if `format-fasta` to
    relabel those files with the option --relabel-map MAP. 

    {id}    : current id
    {scaff} : scaffold name
    {sidx}  : scaffold index after optional sorting
    {gidx}  : gene index in scaffold
    {type}  : feature type (e.g., gene, mRNA, CDS, ...)
    {idx}   : feature index across all scaffolds
    # these ones require tracking hierarchical relations and are not yet implemented
    # {gene}  : gene name
    # {mRNA}  : mRNA/transcript name
    # {midx}  : mRNA/transcript index in scaffold
    # {ctype} : feature type only of children
    # {cidx}  : child type count
    """
    df = table.copy()
    df = df[[0, 3, 8]]
    df.columns = ["scaff", "type", "attrs"]
    df["id"] = subset_gff_table_attrs(df, ["ID"], "attrs")["attrs"].apply(lambda x: x[3:-1])
    df["sidx"] = pd.factorize(df["scaff"])[0] + 1
    df["gidx"] = df.groupby("sidx").cumcount() + 1
    df["idx"] = range(1, df.shape[0] + 1)

    # apply str substitution
    df["label"] = df.apply(lambda row: text.format(**row), axis=1)

    # regular expression to match ID part (e.g., "ID=g1;")
    def replace_id_with_label(row):
        match = re.match(r"(ID=)(.*?)(;.*)", row["attrs"])
        if match:
            return match.group(1) + row["label"] + match.group(3)
        return row["attrs"]

    # apply the function to each row of the DataFrame
    df["attrs"] = df.apply(replace_id_with_label, axis=1)

    # write MAP file
    if mapfile:
        tmp_table = subset_gff_table_attrs(table, ["ID"], 8)
        old_ids = tmp_table[8].apply(lambda x: x[3:-1]).tolist()
        tmp_table = subset_gff_table_attrs(df, ["ID"], "attrs")
        new_ids = tmp_table["attrs"].apply(lambda x: x[3:-1]).tolist()
        id_map = pd.DataFrame({"old": old_ids, "new": new_ids})
        id_map.to_csv(mapfile, sep="\t", index=False, header=None)

    # replace attrs column with 
    table[8] = df["attrs"]
    return table


def relabel_gff_table_scaffolds_subst(table: pd.DataFrame, text: str | None, relabel_names: list[str] | None, mapfile: Path | None) -> pd.DataFrame:
    """..."""
    if text:
        table["scaff"] = table[0]
        table["sidx"] = pd.factorize(table["scaff"])[0] + 1
        old_scaffs = table["scaff"].tolist()
        table[0] = table.apply(lambda row: text.format(**row), axis=1)
        table = table.drop(columns=["scaff", "sidx"])
        new_scaffs = table[0].tolist()

        # write map of old to new names
        if mapfile:
            tmp_table = subset_gff_table_attrs(table, ["ID"], 8)
            old_ids = tmp_table[8].apply(lambda x: x[3:-1]).tolist()
            tmp_table = subset_gff_table_attrs(df, ["ID"], "attrs")
            new_ids = tmp_table["attrs"].apply(lambda x: x[3:-1]).tolist()
            id_map = pd.DataFrame({"old": old_ids, "new": new_ids})
            id_map.to_csv(mapfile, sep="\t", index=False, header=None)


    table[0] = table[0].apply(lambda x: ...)




def run_format_gff(args):
    """..."""

    # allow fuzzy matching of args with options
    to_choices = ["gff3", "bed", "ids"]
    args.to = fuzzy_choice_match(args.to, to_choices)

    # TODO: option to keep the original case, and to allow multiple matches
    # to_choices = ["gene", "mRNA", "CDS", "..."]
    # args.to = fuzzy_choice_match(args.to, to_choices)

    # wrap to allow graceful exit
    try:
        # extract comments
        comments = extract_comments(args.GFF)

        # process whole file to sort chroms
        table = get_gff_table(args.GFF)
        table = sort_gff_table(table, args.sort_alpha, args.sort_len)
        table = subset_gff_table_scaffolds(table, args.subset, args.subset_idx, args.subset_names)
        table = relabel_gff_table_scaffolds_subst(table, args.relabel, args.relabel_names, args.relabel_map)

        # TODO: faster subset chroms w/o sorting by streaming and not reading whole file
        # ...

        # subset features
        table = subset_gff_table_features(table, args.features)

        # takes into account feature types, but requires post-sorting
        if args.filter_duplicates:
            table = filter_duplicates(table)
            # table.sort_values(by=[0, 3, 4], ascending=[True, True, False], inplace=True)

        # relabel ids
        if args.relabel_ids:
            table = relabel_ids_subst(table, args.relabel_ids, args.relabel_ids_map)

        # always subselect "ID" tag for ID output
        if args == "ids":
            table = subset_gff_table_attrs(table, ["ID"])
            ids = table[8].apply(lambda x: x[3:-1]).tolist()
            sys.stdout.write(" ".join(ids) + "\n")
            sys.exit(0)

        # always subselect "ID" tag for BED records
        if args == "bed":
            table = table.loc[:, [0, 3, 4, 8, 5, 6]]
            table = subset_gff_table_attrs(table, ["ID"], 8)
            table.loc[:, 8] = table.loc[:, 8].apply(lambda x: x[3:-1])#.to_list()
            table.to_csv(sys.stdout, sep="\t", index=False, header=None)
            sys.exit(0)            

        # subselect for GFF output
        if args.attrs:
            table = subset_gff_table_attrs(table, args.attrs)            

        # TODO: get comments and append back to output here.
        if not args.exclude_comments:
            if args.add_to_comments:
                cmd = ["# CMD: twig"] + sys.argv[1:]
                comments += " ".join(cmd) + "\n"
            sys.stdout.write(comments)

        # write to TSV
        table.to_csv(sys.stdout, sep='\t', index=False, header=None)

    # allow graceful exit if output is piped but pipe is broken
    except BrokenPipeError:
        sys.exit(0)
    except KeyboardInterrupt:
        logger.error("KeyboardInterrupt by user")
        sys.exit(0)


def main():
    """module-level main cli"""
    parser = get_parser_format_gff()
    args = parser.parse_args()
    run_format_gff(args)



if __name__ == "__main__":
    main()
