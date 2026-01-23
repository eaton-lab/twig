
# macse-prep

The `macse-prep` command is a wrapper tool intended to apply and extend features of
the method `trimNonHomologousFragments` in `macse`. This applies a kmer approach to
compare substrings between sequences to identify and trim fragments or sequences that
do not pass minimum thresholds for homology. This helps to prevent the inclusion of
very low homology sequences in downstream steps where they can greatly affect alignments.
Low homology sequences can end up in orthogroups when using very low thresholds
for homology identification.

This step measures homology separately across full sequences (`-hf`) and the non-edge/internal
parts of sequences (`-hi`) by counting exact matches in subsequences (`-mm` length).
Internal or external sequence fragments that fall below this threshold, and are longer than
some minimum length setting (`-ti` or `-te`) can be trimmed from sequences. Homology is assessed
pairwise between all samples and each fragment must share homology above this min threshold with
a set minimum number of other samples (`-mc`). After non-homologous fragments are removed, if a
sample does not exceed the homology thresholds with others samples it is discarded.

We implement a number of additional filtering steps in addition to macse. Trimmed sequences
that fall below a minimum length (`-ml`) are discarded. Isoforms can also be excluded at this
step by entering an isoform regular expression that will group sequences of the same sample
and gene, but different isoform, by header names. The default looks for Trinity-type naming
conventions with `_i0`, `_i1`, etc. at the end of names. You can likely build a similar expression
for your samples using help from a chatbot. This will select only the isoform that has highest
homology with other samples to keep, and discard extra isoforms. This step can be toggled off
with the `-xi` option.


## The macse-prep subcommand
```bash
twig macse-prep -h
```

```
-------------------------------------------------------------------
| macse-prep: prepare CDS for alignment (trim, filter, iso-collapse)
-------------------------------------------------------------------
| Macse ...
| This implements `macse -prog trimNonHomologousFragments` and
| additional steps to ...
-------------------------------------------------------------------

options:
  -h, --help                                        show this help message and exit
  -i path, --input path                             input CDS (aligned or unaligned)
  -o path, --outpath path                           path to write nt fasta result
  -e [str ...], --exclude [str ...]                 optional names or glob to exclude one or more sequences
  -s [str ...], --subsample [str ...]               optional names or glob to include only a subset sequences
  -ml int, --min-length int                         min nt sequence length after trimming [0]
  -mc int, --min-count int                          min num sequences that must pass filtering to write output [0]
  -hf float, --min-homology-full float              min homology required w/ >=mc others across full sequence [0.1]
  -hi float, --min-homology-internal float          min homology required w/ >=mc others in the internal sequence [0.5]
  -hc int, --min-homology-coverage int              min samples a seq must share homology with at >= mh and mi [3]
  -ti int, --min-trim-length-homology-internal int  trim fragments w/ len <ti and homology <mi w/ <mc sequences [50]
  -te int, --min-trim-length-homology-external int  trim fragments w/ len <tx and homology <mh w/ <mc sequences [50]
  -mm int, --mem-length int                         homology is the prop of aa Maximum Exact Matches of this length [6]
  -is str, --isoform-regex str                      regex used to group isoform sequences ['^([^|]+)\|.*?__(.+?)_i\d+']
  -xi, --skip-isoform-collapse                      skip isoform collapse step
  -v, --verbose                                     print macse progress info to stderr
  -f, --force                                       overwrite existing result files in outdir
  -k, --keep                                        keep tmp files (for debugging)
  -l level, --log-level level                       stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

### Running twig macse-prep
Some commonly used example parameter settings:
```bash
$ twig macse-prep -i CDS -o OUT.nt.fa
$ twig macse-prep -i CDS -hf 0.1 -hi 0.5 -ti 50 -te 50 -hc 15
$ twig macse-prep -i CDS -hf 0.3 -hi 0.8 -ti 25 -te 25 -hc 15
$ twig macse-prep -i CDS -hf 0.5 -hc 10 -k -ml 200 -e '^sppA.*'
$ twig macse-prep -i CDS -s '^sppA.*'
```

## Options

### --input and --outpath
...

### --exclude and --subsample
A list of sample names or glob pattern can be used to select samples by name that should be excluded
from the analysis using `--exclude`, and/or select the subset that you want to be included using
`--subsample`. This is convenient for testing and/or discarding samples.

### --min-length
...

### --min-count
...

###


## Extras

### Why use twig macse-prep
It provides a convenient interface to `macse trimNonHomologousFragments` paired with a number of additional
sequence filtering options to include/exclude samples, filter by length, filter by number of passing sequences,
and select/collapse isoforms by highest homology/length.


### Parallelization
There are two easy ways to parallelize macse-prep over many samples. The first is to
use gnu parallel (easy to install on a linux or mac) and iterate over a set of fasta
files in a directory, and specifying their out path.

```bash
# run parallel jobs on many cds files
$ parallel -j 10 "twig macse-prep -i {} ..."  ::: CDS/*.fa
```

The other way that might be more commonly employed on an HPC cluster is to write a
batch script ...
```bash
...
```

### full pipeline
This step is part of a 3-step pipeline for jointly aligning CDS and AA sequences:
```bash
$ twig macse-prep -i CDS -o TRIM     # {TRIM}
$ twig macse-align -i TRIM -o MSA    # {MSA}
$ twig macse-refine -i MSA -o MSA    # {MSA}.nt.fa, {MSA}.aa.fa
```
