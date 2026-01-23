# Welcome to twig

## Tree-based Workflows for Integrative Genomics.

Twig is a command-line tool for phylogenomic and evolutionary analyses of orthogroups.
It is designed in the fashion of unix utilities, with many small tools each intended
to complete a simple task.
It includes tools to trim and filter sequences and alignments to improve orthology;
to filter, root, and modify gene trees in gene tree distributions;
and wrappers for downstream analyses that take both CDS alignments and
gene trees as inputs.

## Subcommands
Use the `-h` flag to view the subcommands available in `twig`.

```bash
twig -h
```

```literal
-----------------------------------------------------
|  twig: tree-based workflows for integrative genomics |
-----------------------------------------------------

options:
  -h, --help         show this help message and exit
  -v, --version      show program's version number and exit

subcommands:
  --------------     -----------------------------------------------------
    csubst           run csubst in wrapper
    diamond-bl       search blastp hits of one protein fasta to another
    diamond-pw       write .tsv blast hits for all pairs of fasta sequence files
    format-fasta     format fasta to sort/filter/subset/relabel/modify sequences
    format-gff       format gff to filter and relabel genes and scaffolds
    genome-table     write table with chrom/scaff [names lengths ngenes]
    filter-concat    filter loci and genes, concatenate down or right, convert .fa/.phy/.nex
    macse-prep       prepare CDS for macse alignment (trim,filter,iso-collapse)
    macse-align      run macse alignment (CDS/AA)
    macse-refine     trim/filter/convert alignments
    partition-cds    write a partition file for phylogenetic inference
    tree-filter      filter and relabel a set of trees for downstream analyses
    tree-rooter      Re-root gene trees by outgroup, species tree, and/or MAD
    tree-skeleton    Force gene tree to match species tree topology
```

## Calling subcommands

Each subcommand has its own help menu describing the utility, its parameters, and with examples.
See the documentation for further examples.

```bash
twig tree-filter -h
```

```
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

options:
  -h, --help                                     show this help message and exit
  -i path, --input path                          newick or multi-newick trees file
  -o path, --out path                            outfile name else printed to stdout
  -d str, --delim str                            delimiter to split tip labels
  -di int [int ...], --delim-idxs int [int ...]  index of delimited name items to keep
  -dj str, --delim-join str                      join character on delimited name items
  -I path, --imap path                           filepath listing names to keep, or assigning (name population)
  -M path, --minmap path                         filepath listing (population mincov) for filters
  -m int, --min-tips int                         min tips after pruning [4]
  -c int, --max-copies int                       filter trees with >c gene copies from a taxon [None]
  -eo float, --edge-outlier-outgroup float       exclude 'outgroup' population edges if >eo stdev from mean [10]
  -ei float, --edge-outlier-ingroup float        exclude non 'outgroup' population edges if >ei stdev from mean [5]
  -rd, --relabel-delim                           relabel tips by their delim parsed names
  -ri, --relabel-imap                            relabel tips to their imap mapped names
  -S, --subsample                                subsample to include only tips in imap
  -E, --exclude-outliers                         exclude tips with outlier edge lengths (>ei or >eo)
  -R, --require-outgroups                        require at least one 'outgroup' sample
  -C, --collapse-outgroups                       keep only the most distant 'outgroup' (assumes rooted trees)
  -O path, --outgroups path                      path listing outgroups for -R (overrides IMAP group 'outgroup')
  -l level, --log-level level                    stderr logging level (DEBUG, [INFO], WARNING, ERROR)

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
```
