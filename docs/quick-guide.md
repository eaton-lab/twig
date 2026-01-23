
# quick guide


## Example dataset
...


## The fastest guide
You can perform an ortholab analysis to assemble and analyze orthogroups
using the all-in-one analysis command `orthogroup run`. This takes a large
number of options that can be used to specify the run. It will identify
orthology among sequences to group them into orthogroups; align and trim 
orthogroup sequences; infer gene trees; infer a species tree; root gene 
trees and species trees; and report a number of evolutionary analyses of
orthogroups.

```bash
# run the full analysis in one or few commands
ortholab run data/*.faa 
```

<!-- ## Orthogroup identification -->
## Sequence orthology
A primary goal of ortholab is to identify homology among sequences both
within and between biological samples to enable evolutionary analyses.

### Fetch transcriptome datasets
See the data-guide (link) for tips on assembling transcriptome datasets
for use in ortholab, and/or accessing publicly accessible datasets from
online repositories.

```bash
# download an example dataset
ortholab seqs download --database phytozome --accession 23456 --all
```
```
# filetree of example dataset (collapsible)
...
```

### Preprocess data
Many datasets will vary in ... One major aspect is in the sequence labels
store in the protein or cds fasta files (also referred to as fasta headers).

```bash
# create a map ..
ortholab seqs create-relabel-map data/*/*.proteins.faa.gz > map_proteins.tsv

# optionally also crate one for the CDS
ortholab seqs create-relabel-map data/*/*.cds.fna.gz > map_cds.tsv
```

If you create both it is important that the label columns match between the two.
```
# show map TSV collapsible
```

### Homology scores

There are many options that can be implemented in this tool to modify which 
program is used to measure homology, and what the minimum cutoff will be. See
the tools tutorial page (link). This will use the map file to relabel sequences
as they are passed into the blast database.
```bash
ortholab seqs homology-score \
    data/*/*.proteins.faa.gz \
    --map map.tsv \
    --tool diamond-blast \
    --max-evalue 0.001 \
    --outdir blastdir/ \
    --use-labels
```

The blast table format can be modified with options to the tool above,
but we restrict it such that the first X columns will always contain the 
same information, and additional requested features are appended as the
final columns.

From looking at this table you can see that we have scored hits between
sequences sometimes with *very* low sequence similarity scores. This is
sometimes desirable. In the cases where we do not want to allow such 
deep homology, don't worry, we can filter the hits more strictly in the
next step. However, if you modify run parameters of this step, such as
changing the max-evalue or ... , it can have the effect of making the
homology search run faster.
```
# show blast table collapsible
```

### Homology filtering

```bash
# filter homology hits to write a graph
ortholab seqs homology-graph \
    blastdir/*.tsv \
    --min-identity 0.5 \
    --max-... 5 \
    --min-... 0.1 \
    ... \
    --length-normalize \
    --require-rbbh \
    --outdir blastdir \
    --prefix test
```

The graph is represented by pairs of sequence labels and a weight value.
This graph file can be visualized in a tool like ... 
Following the approach used by OrthoMCL and Orthofinder, we will split
this graph into smaller orthogroups using a Markov Clustering method 
in `mcl`. This takes an inflation parameter `-i` which ...

```
# show blast table collapsible
```

### Orthogroup identification
Next we write each group of orthologous sequences (orthogroup) to a file.
This tool will take an orthogroup set split from the graph in the
previous step and write a orthogroup

```bash
# filter homology hits to write a graph
ortholab seqs homology-split \
    --graph blastdir/homology-graph.tsv \
    --databases blastdir/blast*.tsv \
    --map map.tsv \
    --format [json,tsv] \
    > blastdir/orthogroups-i-1.5.tsv
```

### Orthogroup refinement
In many other software tools you now simply proceed to use your orthogroup
sequences in downstream analyses and assume that homology was accurately assessed.
However, in *ortholab* we provide a number of additional tools that allow you to
incorporate information you can learn from downstream analyses to feedback into the
orthology assessment in order to refine your orthogroup assignments. This can be
an iterative process. 

The two primary sources of additional information include gene tree relationships
and synteny. 
```
TODO
``` 

### Orthogroup sequence trimming & alignment
Before running evolutionary analysis of orthogroup sequences we should filter and align 
them. Filters can be used to skip files that do not have taxa of interest. 

```bash
# filter homology hits to write a graph
ortholab seqs orthogroup-align \
    --orthogroups/*.faa \
    --map map.tsv \
    --cds ... \
    --align-tool macse2 \
    --trim-tool ... \
    # require taxon A,B,C (use if --map then use mapped names)
    --sample-require A B \ 
    # exlude sequence with too many variant sites, indels?
    --sequence-filter ... \  # 
```




## Tree analyses
Now that we have assembled sequences into orthogroups we can proceed to infer the evolutionary relationships among sequences. 

```bash
ortholab trees infer-gene-trees \
    --model GTR+G+F4 \
    --bootstraps 100 \
    --max-copy-number 10 \
    --max-size 1000 \
    --cores 40 \
    --threads 10 \
    --outdir gene-trees/
```

```bash
```

```bash
ortholab trees draw gene-trees/OG1.nwk 



## Evolutionary analyses
...

## Utilities
...