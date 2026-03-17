# align-cds

## The `align-cds` subcommand

`twig align-cds` wraps `macse -prog alignSequences` to write an intermediate
codon-aware CDS alignment. The output from this step should usually be passed to
`twig align-post` for trimming and filtering.

```bash
twig align-cds -i CDS.nt.fa -o OUT.msa.nt.fa
twig align-cds -i CDS.nt.fa -o OUT.msa.nt.fa -m 100
twig align-cds -i CDS.nt.fa -o OUT.msa.nt.fa --seq-lr '^pseudo_' '^partial_'
twig align-cds -i CDS.nt.fa -o OUT.msa.nt.fa --seq-lr seq_lr_patterns.txt
```

## Less-reliable sequences with `--seq-lr`

Use `--seq-lr` when some input records should be aligned as MACSE
less-reliable sequences via `-seq_lr`.

- Inline mode: pass one or more regex selectors directly.
- File mode: pass one existing file path; each nonblank, non-comment line is
  treated as one selector.
- Matching uses Python regular expressions with `re.search()` against each FASTA
  header.

When `--seq-lr` matches one or more records, `twig` writes temporary FASTAs for
reliable and less-reliable sequences and passes both to MACSE. If nothing
matches, `twig` warns and falls back to a normal `align-cds` run without
`-seq_lr`.

## Options

```text
-i, --input path             input CDS fasta (aligned or unaligned)
-o, --outpath path           output path for aligned nucleotide fasta
-m, --max-refine-iter int    maximum MACSE refinement iterations (-1 uses MACSE default)
--seq-lr str [str ...]       regex selectors or one selector-file path for less-reliable sequences
-v, --verbose                print MACSE progress to stderr
-f, --force                  overwrite existing output files
-l, --log-level level        stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```
