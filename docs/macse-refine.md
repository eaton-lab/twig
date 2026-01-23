
# macse-refine
After sequences have been aligned the `macse-refine` step can be used to trim edges,
mask frameshift or stop codons, perform quality filtering, and write the final
CDS to proteins alignments to file. Edge trimming, masking, and writing is performed
using a wrapper around `macse ...`. We provide an additional filter to examine overlap
among sequences in an alignment and enforce a minimum overlap between samples. This
is very important for preventing terraces during phylogenetic inference, where some
samples share no phylogenetic information. This can arise in orthogroups often when
they were clustered together based on different shared motifs. Sequences that do not
overlap with at least N other sequences by X amount of nucleotides are iteratively
removed until all meet the minimum requirements, or the alignment fails quality check.


## The twig macse-refine subcommand
```bash
twig macse-refine -h
```

```
-------------------------------------------------------------------
| macse-refine: CDS/AA jointly and filter low homology seqs
-------------------------------------------------------------------
| Macse codon-aware alignment of CDS and AA sequences. This method
| takes an alignment from `twig macse-refine` and can subsample/
| filter to remove selected samples, or those with too little
| overlap with others (-mo), then refine the alignment, and write
| final alignments with optional conversion of frameshift and stop
| codons.
-------------------------------------------------------------------

options:
  -h, --help                                     show this help message and exit
  -i path, --input path                          input aligned CDS
  -o path, --outprefix path                      out prefix; creates {prefix}.nt.fa and {prefix}.aa.fa
  -e [str ...], --exclude [str ...]              optional names or glob to exclude one or more sequences
  -s [str ...], --subsample [str ...]            optional names or glob to include only a subset sequences
  -t path, --tree path                           optional newick file to subsample genes present in tree
  -mo int, --min-overlap int                     min sites shared by each pair of samples, else iter prune lowest [0]
  -ms int, --min-samples int                     min samples after filtering else error [0]
  -ac float, --aln-trim-ends-min-coverage float  trim alignment edges to where a min percent of samples have data [0.4]
  -as int, --aln-trim-window-size int            trim alignment edges using a sliding 'half_window_size' [5]
  -if str, --codon-int-fs str                    codon to sub for internal frame shift [NNN]
  -ef str, --codon-ext-fs str                    codon to sub for external frame shift [NNN]
  -fs str, --codon-final-stop str                codon to sub for final stop [NNN]
  -is str, --codon-int-stop str                  codon to sub for internal stop [NNN]
  -r, --refine-alignment                         refine alignment
  -R, --refine-alignment-if                      refine alignment only if >=1 sequences are filtered out
  -ri int, --max-iter-refine-alignment int       max iterations in refine alignment [-1]
  -v, --verbose                                  print macse progress info to stderr
  -f, --force                                    overwrite existing result files in outdir
  -k, --keep                                     keep tmp files (for debugging)
  -l level, --log-level level                    stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

### Running twig macse-refine
```bash
Examples
--------
$ twig macse-refine -i CDS -o PATH/PRE -ac 0.5
$ twig macse-refine -i CDS -o PATH/PRE -t NWK -R
$ twig macse-refine -i CDS -fs XXX -is XXX
```

### Options
...


### Parallelization
```
# run parallel jobs on many cds files
$ parallel -j 10 "twig macse-refine -i {} ..."  ::: CDS/*.msa.nt.fa
```


### Full pipeline
```
# full pipeline
$ twig macse-prep -i CDS -o TRIM        # {TRIM}
$ twig macse-align -i TRIM -o ALN       # {ALN}
$ twig macse-refine -i ALN -o MSA       # {MSA}.nt.fa, {MSA}.aa.fa
```

