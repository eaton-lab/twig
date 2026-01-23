# installation


## Conda installation
The preferred installation method is by using `conda` which will install
the minimal required dependencies along with the `twig` tool. To install
all

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

```
- csubst
- ...
```
