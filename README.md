
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
twig [subcommand]

subcommands
-----------
    format-data     ...
    format-fasta    ...
    format-genome   ...
    diamond-blastx  ...
    macse-prep      ...
    macse-align     ...
    macse-refine    ...
    tree-filter     ...
    tree-rooting    ...
    alien-index     ...

pipelines?
-----------
    ortholab        ...
    synteny         ...
```

### Subcommands
