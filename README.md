
# TWIG

### IN DEVELOPMENT (use with caution)


### Tree-based workflow for integrative genomics
`twig` provides a toolkit of commands for comparative genomic analyses.


### Installation
While in development the best way to install twig is to install dependencies with
conda and the install twig from GitHub using pip.
```bash
# install dependencies into a conda env
conda install toytree diamond macse -c conda-forge -c bioconda

# clone the repo
git clone git@github.com:eaton-lab/twig.git

# install local source code in development mode
cd twig/
pip install -e . --no-deps

# test it
twig -h
```

### Top-level command
```bash
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

### Subcommands
