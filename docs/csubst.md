
# csubst

The `csubst` method by Fukushima and Pollack (2024) is used to detect evidence of adaptive
molecular convergence given a CDS alignment and gene tree. It is implemented in a
command-line tool [csubst GitHub link] that provides many options. The `twig csubst` tool
is a simple wrapper around this tool to make it easier to run in alternative ways.

## The csubst subcommand
```bash
twig csubst -h
```

```
usage: twig csubst -i MSA -t TREE -g FG [options]

-------------------------------------------------------------------
| csubst: A simple wrapper to run and extract results from csubst
-------------------------------------------------------------------
| ...
-------------------------------------------------------------------

options:
  -h, --help                      show this help message and exit
  -a path, --alignment path       input CDS alignment
  -t path, --tree path            input rooted tree file
  -g path, --foreground path      foreground file
  -o path, --outdir path          Out directory will be created at this path [./csubst]
  -m int, --max-arity int         max combinatorial number of branches (K) [2]
  -u int, --exhaustive-until int  perform exhaustive (non-heuristic) search up N branch combs [2]
  -c str, --cutoff-stat str       Cutoff stats for searching higher-order branch combs
                                  [OCNany2spe,2.0|omegaCany2spe,5.0]
  -b, --foreground-table          foreground file is a table (fg_format=2)
  -w, --fg-exclude-wg             sets 'yes' to exclude within-group comparisons
  -s, --fg-stem-only              sets 'yes' to only stem comparisons
  -e path, --env path             conda env name where 'csubst' in installed [csubst]
  -j int, --threads int           number of threads
  -v, --verbose                   print macse progress info to stderr
  -f, --force                     overwrite existing result files in outdir
  -l level, --log-level level     stderr logging level (DEBUG, [INFO], WARNING, ERROR)

Examples
--------
$ twig csubst -a MSA -t NWK -g FG -o OUT/PRE
```

### Running twig csubst

```bash
twig csubst -i MSA -t TREE -g FG ...
```

### Why use twig csubst

I find it to be convenient for the following two reasons:

First, the current installation of csubst enforces some installation constraints, including on the version of numpy, and of the external iqtree tool. For this reason, it is good practice to create a unique conda environment just for installing csubst and iqtree2. You can specify that environment's name using the `-e` flag from `twig csubst` and it will run csubst from that environment without having to activate it. This way you can call and run this tool from your modern updated conda environment.

Second, csubst generates output files in its current working directory. Which requires you to
change directories a lot when running many files. A minor inconvenience. Instead you can set an
output directory path with `-o` in `twig csubst`.

That's it for now.
