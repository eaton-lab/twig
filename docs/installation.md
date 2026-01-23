# installation



## Install from GitHub
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


## Conda installation
The preferred installation method is by using `conda` which will install
the minimal required dependencies along with the `twig` tool.

[THIS IS NOT YET AVAILABLE]
```bash
conda install twig -c conda-forge
```

## Dependencies
The number of dependencies in twig is purposefully kept minimal to make
installation simple and reliable.
The core dependencies are listed below.

```
- toytree
- pandas
- numpy
- scipy
- loguru
- diamond
- ...
```

## Optional dependencies
There are also several optional tools that can be called from
twig but which are not installed by default. This includes alternative
software options for assessing homology among sequences and inferring 
gene trees. If these options are selected during an ortholab analysis
it will ask you to install the required dependency by providing a conda
command. 
