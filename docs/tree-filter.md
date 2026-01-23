
# tree-filter

The `twig tree-filter` tool can be used to filter one tree, or a set of trees, based on
a set of parameters. You can collapse low-support branches, subselect samples, rename
samples; filter trees by min number of tips or min number of resolved edges; drop edges
that have outlier edges lengths; collapse outgroup clades; exclude trees with paralogs;
and more.


## The tree-filter command


## Running tree-filter


## Renaming tips


## Subselect tips


## Filtering edges by length


## Filtering edges by support


## Filtering trees by ntips


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
