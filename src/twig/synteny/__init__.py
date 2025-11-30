#!/usr/bin/env python

"""Synteny analysis subpackage

1. Generate/format GFF/PEP data files
2. Perform blast
3. dotplot --draw --write
4. refine collinearity using MCScan algorithm
5. Measure Ka/Ks (can we do this w/o calling PAML?) or use HyPhy?
6. BlockInfo?

$ format-genome FA --max-width 80 --sort-len --subset-idx 1 2 3 --relabel-prefix sub > FILE
Write a new genome fasta. Sort scaffolds by current name or length; 
relabel as numeric w/ prefix; wrap to max width; subset scaffolds.

$ format-gff GFF --relabel 
Write a new gff with conversion to gff3; relabel genes as numeric
with prefix; subset to fewer features; write a map file for converting
back to original names.

$ format-fasta FA --max-width 80 --relabel --map ...
Write a new fasta with sequences relabeled using a headers map file; 
set to max-width; ...

$ format-data FA GFF PIP CDS --relabel --relabel-prefix subA --subset-scaffolds 1 2 3 --sort-len --out DATA/subA --max-width 80
Performs the above tasks on all genome dataset files.

$ format-genome-table FA GFF --subset-scaffolds 1 2 3
Writes a genome table with scaffold names, lengths, and ngenes.

$ diamond-blastp ...
$ lastal ...
$ minimap2 ...
...

$ collinear TSV GFF1 GFF2 TABLE1 TABLE2 --write ... --draw ... --draw-args ...
    --max-homologous-hits 1 \
    --min-bitscore 100 \
    --max-evalue 1e-5 \
    --max-repeat-number 10 \
    --position order \
    --reorder1 +1 -2 +3 +4 ... # should we?
    > 

$ collinear ... | dotplot - --label1 A --label2 B --size 5 --color-bins 50 40 20 --width 400 --height 400 --png 

$ collinearity-mcscan TSV \
    --max-homologous-hits 1 --max-evalue None --min-bitscore None \
    --color-bins 50 40 25 \
    --max-gap 40 40 --pvalue 1 --repeat-number 10 \
    --threads 10 \
    > collinear-table.tsv
"""

# from .format_fasta import get_parser_format_fasta
# from .genome_table import get_parser_genome_table

