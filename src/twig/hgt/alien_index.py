#!/usr/bin/env python

"""Calculate the alien index using the NCBI NR protein database.

Methods
-------
...

Examples
--------

# blastp all samples to each other and to NR
$ twig diamond-blastp PEP1 DB -k 5 --kwargs-blastp '--taxonlist 33090' > blast_1_nr.tsv
$ twig diamond-blastp PEP1 PEP2 > blast_1_2.tsv

# get alien index for all genes in each
$ alien-index -o blast_1_nr.tsv -i blast_1_2.tsv -a blast_1_2.tsv > alien_index_1.tsv

```
"""

import sys
import gzip
import textwrap
from pathlib import Path
from functools import lru_cache
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from subprocess import Popen, PIPE, CalledProcessError
import numpy as np
import pandas as pd
from loguru import logger
from twig.utils.make_wide import make_wide

logger = logger.bind(name="twig")

NAMES = [
    "query", "subject", "hit_pct", "len", "mismatch", 
    "gapopen", "qstart", "qend", "sstart", "send", 
    "evalue", "bitscore", "cbitscore", "qlen", "slen",
]

INDEX = [
    "query", "qlen", "alien", 
    "taxon_o", "gene_o", "hit_pct_o", "bit_o",
    "taxon_i", "gene_i", "hit_pct_i", "bit_i",
    # "taxon_o", "gene_a", "hit_pcts_a", "bits_a",
    # "og", "gene_o", "taxon_o",
]


KWARGS = dict(
    prog="alien-index",
    usage="alien-index --ingroup A.gz B.gz --outgroup C.gz [options]",
    help="metric to detect HGT from comparative blastp scores",
    formatter_class=make_wide(RawDescriptionHelpFormatter),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | alien-index: write .tsv of comparative bit-scores in groups     |
        -------------------------------------------------------------------
        | Alien index is measured by max-bit-score to an outgroup sample  |
        | minus max-bit-score to an ingroup sample. Negative values are   |
        | of interest as potential HGT events. It is often advisable to   |
        | use a blastp search against (a subset of) the NR database as    |
        | as the outgroup table. The ingroup table should be a blastp     |
        | search result against a smaller subset of the NR table, option- |
        | ally supplemented by searches against additional ingroups. You  |
        | can use the --additional flag to include bit scores to other    |
        | samples in the output table for further downstream comparisons. |
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        $ alien-index -o A_NR.tsv.gz -i A_IN.tsv.gz > AI.tsv
        $ alien-index -o A_NR.tsv A_B.tsv A_C.tsv -i A_IN.tsv > AI.tsv
        $ alien-index -o A_NR.tsv -i A_IN.tsv -a A_B.tsv A_C.tsv > AI.tsv
        $ alien-index -o A_NR.tsv -i A_[B,C].tsv -a A_D*.tsv > AI.tsv        
    """)
)


def get_parser_alien_index(parser: ArgumentParser | None = None) -> ArgumentParser:
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
    parser.add_argument("-i", "--ingroup", type=Path, metavar="path", nargs="*", required=True, default=None, help="one or more TSV blastp results of query to an ingroup (gzip OK, regex OK)")
    parser.add_argument("-o", "--outgroup", type=Path, metavar="path", nargs="*", required=True, default=None, help="one or more TSV blastp results of query to an outgroup (gzip OK, regex OK)")
    parser.add_argument("-a", "--additional", type=Path, metavar="path", nargs="*", default=None, help="one or more TSV blastp results of query to other samples (gzip OK, regex OK)")
    parser.add_argument("-e", "--exclude-prefix", type=str, metavar="str", nargs="*", default=None, help="exclude hits to a genes that start with one or more listed prefixes")
    # parser.add_argument("-e", "--efetch-nr-info", action="store_true", help="optional: efetch gene and taxon names for NR hits")
    # parser.add_argument("-b", "--efetch-binary", type=Path, metavar="path", default=None, help="optional: efetch binary path (if not in PATH)")
    parser.add_argument("-p", "--positive", action="store_true", help="include positive AI values in output (otherwise excluded)")
    parser.add_argument("-c", "--compress", action="store_true", help="gzip compress output")
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")

    return parser


######################################################################


def _expand_path_group(paths: list[Path]) -> list[Path]:
    """Return a list of Paths expanded from one or more paths or regex.
    """
    epaths = []
    for path in paths:
        expanded_paths = list(path.parent.glob(path.name))
        # TODO: could allow empty paths and just print warning.
        assert expanded_paths, f"no files matched {path}"
        for epath in expanded_paths:
            if not epath.exists():
                logger.warning(f"no files matched path {epath}")
            else:
                epaths.append(epath)
    return epaths


def _paths_to_path_dicts(paths: list[Path]) -> dict[str, Path]:
    """Return {name: path} dict where name strips suffix and 'blastp' from filenames.
    """
    names = [i.name.split(".")[0].lstrip("blastp_") for i in paths]
    names = ["_".join(i.split("_")[1:]) for i in names]
    return dict(zip(names, paths))  # {name: i for i in paths}


def _path_dict_to_dataframes(path_dict: dict[str, Path], label: str) -> list[pd.DataFrame]:
    """Return concatenated dataframe of labeled rows from all Paths.
    """
    # iterate over the paths
    data = []
    for name, path in path_dict.items():

        # check for header
        header = False
        xopen = gzip.open if path.suffix == ".gz" else open
        with xopen(path, 'r') as indata:
            if all(i in indata.readline().split() for i in ["query", "subject", "bitscore"]):
                header = True

        # parse the dataframe 
        if header:
            df = pd.read_csv(path, sep="\t", header=NAMES)
        else:
            df = pd.read_csv(path, sep="\t", header=None)

        # add headers if it is the right size
        if df.shape[1] != 15:
            raise Exception("blastp TSV is not in expected twig format (15 columns)")
        if not header:
            df.columns = NAMES

        # clean query and subject for bad header names (add more here)
        df['query'] = df['query'].apply(lambda x: x.split("\\t")[0])        
        df['subject'] = df['subject'].apply(lambda x: x.split("\\t")[0])
        
        # add labels of name and group
        df['name'] = name
        df['group'] = label
        data.append(df)
    return data


def _get_ogid(query: str, ogs: pd.Series) -> int:
    """Return Orthogroup ID."""
    return np.where(ogs.str.contains(query))[0][0]


@lru_cache
def _get_gene_info(query: str) -> tuple[str, str, str]:
    """Return Ensembl info for a gene ID using NCBI E-utilities."""
    cmd1 = ["esearch", "-db", "protein", "-query", query]
    cmd2 = ["efetch", "-format", "docsum"]
    cmd3 = ["xtract", "-pattern", "DocumentSummary", "-element", "Organism"]
    proc1 = Popen(cmd1, stdout=PIPE)
    proc2 = Popen(cmd2, stdin=proc1.stdout, stdout=PIPE)
    proc3 = Popen(cmd3, stdin=proc2.stdout, stdout=PIPE)
    proc1.stdout.close()
    proc2.stdout.close()
    try:
        out, err = proc3.communicate()
    except CalledProcessError:
        decoded = out.decode().strip()
        if not decoded:
            raise ValueError(f"No results returned for query: {query}")
    return out.decode()


def _get_score(group: pd.DataFrame, additional: list[Path]) -> pd.Series:
    """Return a pd.Series with alien index for this gene and additional info.

    ... (usually NR database). In the
    case that a sample itself is in the database the blast search
    should have kep more than the closest 1 match. So we 
    """
    # must have a match in the outgroup.
    odf = group.loc[group['group'] == 'o'].reset_index(drop=True)
    if not odf.size:
        return pd.Series(index=INDEX, data=np.nan)
    # odf = group.loc[group.groupby('query')['cbitscore'].idxmax()]

    # must have a match in the ingroup
    idf = group[group['group'] == 'i'].reset_index(drop=True)
    if not idf.size:
        return pd.Series(index=INDEX, data=np.nan)

    # store the query and its length
    query, qlen = odf.loc[0, ['query', 'qlen']]

    # store top hit to outgroup and meta info
    oidx = odf['cbitscore'].idxmax()
    otax, ogene, ohit, obit = odf.loc[oidx, ['name', 'subject', 'hit_pct', 'cbitscore']]

    # store top hit to outgroup and meta info
    iidx = idf['cbitscore'].idxmax()
    itax, igene, ihit, ibit = idf.loc[iidx, ['name', 'subject', 'hit_pct', 'cbitscore']]    

    # compute alien score
    alien_score = obit - ibit

    # add data for 'additionals'
    aresults = []
    aindex = []
    for tax in additional:
        if group["name"].str.contains(tax).any():
            agene, ahit, abit = group.loc[group['name'] == tax, ["subject", "hit_pct", "cbitscore"]].iloc[0]
        else:
            agene, ahit, abit = np.nan, np.nan, np.nan
        aresults += [agene, ahit, abit]
        aindex += [f"gene_{tax}", f"hit_{tax}", f"bit_{tax}"]

    # ...
    result = [
        query, qlen, alien_score,
        otax, ogene, ohit, obit,
        itax, igene, ihit, ibit,
    ]
    return pd.Series(result + aresults, index=INDEX + aindex)


######################################################################


def get_concatenated_df(outgroups: list[Path], ingroups: list[Path], additional: list[Path]):
    """...
    """
    # expand path args to dicts of {name: Path, ...}
    idict = _paths_to_path_dicts(_expand_path_group(ingroups))
    odict = _paths_to_path_dicts(_expand_path_group(outgroups))
    adict = _paths_to_path_dicts(_expand_path_group(additional)) if additional else {}

    # convert each Path to a list of DataFrames and concatenate
    idf = _path_dict_to_dataframes(idict, "i")
    odf = _path_dict_to_dataframes(odict, "o")
    adf = _path_dict_to_dataframes(adict, "a")
    data = pd.concat(idf + odf + adf).sort_values(by="query").reset_index(drop=True)
    return data


def run_alien_index(args):
    """..."""
    # get concatenated dataframe of [file, group] labeled blast results
    data = get_concatenated_df(args.outgroup, args.ingroup, args.additional)

    # remove any perfect hits to self. If a sample is in the NR database
    # then it will always hit self best, these samples should have been
    # run to keep >1 hit using -k during blastp search so other hits are
    # present. Then best non-self based on cbitscore will be kept later.
    mask1 = data['query'] != data['subject']
    mask2 = data['hit_pct'] != 100.0
    mask = mask1 & mask2
    data = data.loc[mask]

    # filter additional hits by prefix in subject of 'o' and 'i' hits
    masks = []
    for pre in args.exclude_prefix:
        phits = data['subject'].str.startswith(pre).values
        odb = data["group"] == "o"
        idb = data["group"] == "i"
        masks.append((phits & odb) | (phits & idb))
        logger.debug(f"prefix {pre} occurred {phits.sum()} times in subject")
    if masks:
        mask = np.sum(masks, axis=0, dtype=np.bool_)
        data = data.loc[~mask]
        logger.debug(f"{mask.sum()}/{data.shape[0]} rows removed")

    # extract list of 'additional' taxa not used for alien metric
    adds = sorted(data.loc[data["group"] == "a", "name"].unique())

    # calculate alien metric using ingroup and outgroup groups
    adata = [_get_score(j, adds) for i, j in data.groupby("query")]
    
    # filter for negative aliens
    adf = pd.DataFrame(adata)
    if not args.positive:
        adf = adf.loc[adf["alien"] > 0.]
    
    # sort by alien score
    adf = adf.sort_values(by="alien", ascending=False).reset_index(drop=True)

    # add NR gene info
    # if args.efetch_nr_info:
    #     ...
    return adf


def main():
    parser = get_parser_alien_index()
    args = parser.parse_args()
    adf = run_alien_index(args)
    adf.to_csv(sys.stdout, sep="\t", float_format="%.2f", index=False)



if __name__ == "__main__":
    main()

    # import sys
    # # out, err = Popen(["which", "xtract"], stdout=PIPE).communicate()
    # res = _get_gene_info("KAH6780014")
    # print(res)
    # print("HI")
    # sys.argv.extend(["-i", "/home/deren/Documents/projects/orthoparasites/data/blastp_Pcra_S*.tsv.gz"])
    # sys.argv.extend(["-o", "/home/deren/Documents/projects/orthoparasites/data/blastp_Pcra_NR_only_33090.tsv.gz"])
    # sys.argv.extend(["-a", "/home/deren/Documents/projects/orthoparasites/data/blastp_Pcra_P[j,a]*.tsv.gz"])
    # main()
