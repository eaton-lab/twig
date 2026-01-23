
# tree-filter

The `twig tree-filter` tool can be used to filter one tree, or a set of trees, based on
a set of parameters. You can collapse low-support branches, subselect samples, rename
samples; filter trees by min number of tips or min number of resolved edges; drop edges
that have outlier edges lengths; collapse outgroup clades; exclude trees with paralogs;
and more.


## The tree-filter command
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
```


## I/O: --input and --output
The input to tree-filter is a text file containing one or more newick trees, and the outout is written either to a specified file path `--out` or to STDOUT. If multiple tree are present in the input file the output will write the trees that pass filtering in the same order. A stats summary of how many trees passed filtering and the reasons that each failed is written to STDERR. A useful workflow to capture the output on STDOUT and the stats log on STDERR is to run:
```bash
twig tree-filter -i NWK 1> OUT 2> LOG
```

## Renaming tips
There are many options for relabeling tips using tree-filter, which can be used not only to write trees with new tip labels, but also to perform operations/filters on trees based on both taxon and gene names that can be parsed from labels. ...
```bash
...
```

## Subselect tips
...

## Filtering edges by length
...

## Filtering edges by support
...

## Filtering trees by ntips
...

## Filtering trees by nedges


## Filtering trees by duplicates (paralogs)


## Collapsing outgroup clades


## Requiring outgroups
The -R and -O flags are used together to filter trees by the requirement that one or more
outgroups is present. Outgroup names are specified in a file selected by -O. If the names
are taxon names (e.g., SppA), then you can use -d, -di, and dj to parse tip names to sample
only the relevant taxon name part, in the case that tip labels are longer gene names (e.g., SppA-geneX-isoformY).

```bash
twig tree-filter -i NWK -d "-" -di 1 -O OUTGROUPS.txt -R 1>NWK2 2>LOG
```

A summary of the filtering is printed to stderr.
```
{tree-filtered: 0, ...}
```
And the filtered tree or trees is printed to stdout. If no trees passed filtering then nothing will be printed to stdout.
```
((...))
```
