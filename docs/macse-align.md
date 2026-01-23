
# macse-align

...

## The macse-align subcommand

```bash
-------------------------------------------------------------------
| macse-align: run macse joint cds/aa alignment
-------------------------------------------------------------------
| Calls `macse -prog AlignSequences` and writes cds alignment
| to {outpath}. This is an intermediate file. You should next
| run `twig macse-refine` to trim/filter/convert to write final
| cds and aa alignments.
-------------------------------------------------------------------

options:
  -h, --help                     show this help message and exit
  -i path, --input path          input CDS (aligned or unaligned)
  -o path, --outpath path        path to write aligned nt fasta
  -m int, --max-refine-iter int  max refinement iterations during optimizing [default -1 = no limit]
  -v, --verbose                  print macse progress info to stderr
  -f, --force                    overwrite existing result files in outdir
  -l level, --log-level level    stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

### Running twig macse-align
```
Examples
--------
$ twig macse-align -i CDS -o /OUT/CDS.msa
```

### Parallelization...
```bash
# run parallel jobs on many cds files
$ parallel -j 10 "twig macse-align -i {} ..." ::: CDS/*.nt.fa
```

### full pipeline
```bash
# full pipeline
$ twig macse-prep -i CDS -o TRIM     # {TRIM}
$ twig macse-align -i TRIM -o MSA    # {MSA}
$ twig macse-refine -i MSA -o MSA    # {MSA}.nt.fa, {MSA}.aa.fa
```

