
# TWIG


# MICROSYNTENY PHYLOGENETIC ANALYSIS
https://www.nature.com/articles/s41467-021-23665-0

```bash
twig [subcommand]

subcommands
-----------
    format-data     ...
    format-fasta    ...
    format-genome   ...
    diamond-blastx  ...
    macse-pipe      ...
    tree-filter     ...
    tree-rooting    ...
    alien-index     ...
    

pipelines
-----------
    ortholab        ...
    synteny         ...
```


### Modules
- [a] FA/multiFA/GFF manipulation and matching.
- [o] ORTHOLOGY + PHYLO inference
- [s] SYNTENY analysis
- [p] CDS/FA prep [align/trim/frame/mask/SFS/...]
- [e] CDS/FA analysis [iMK,Convergence,]


### DEVELOP
Goal: Develop new tool as extension of OrthoFinder/OrthoMCL.
1. Improve deep homology inference (not always split by first species-MRCA)
2. Improve inference of HGTs by refining clusters.
    a. Use gene tree info (depth, missing, etc.)
    b. Use evol rate info (slow fast)
    c. Use edge trimmed sequences?
3. Readable, documented, and commented code base.
4. Evolutionary analysis of sequences: Identify convergence.


### TODO
[ ] transition cli to Typer
[ ] in get_cognates write database files as chunks.
[ ] explore min-bit-score in get_bitscores called from get_cognates funcs
[ ] add Arabidopsis!?
[ ] download script from phytozome? [no?]
[ ] mpi4py drop-in
[ ] cluster info file(s) for debugging
[ ] alternative trimmers: https://github.com/LKremer/MSA_trimmer?tab=readme-ov-file
[ ] PLM Search/Align https://www.nature.com/articles/s41467-024-46808-5?fromPaywallRec=false
[ ] masking of off-domain regions


### Structure/Synteny analysis
[ ] PC Embedding of distances
[ ] Within sample genes are assigned a distance to all other genes
[ ] Across sample genes are assigned a distance to their cognates, or zero
[ ] ...


### TEMPLATE modular
```bash
$ elab format-fasta ...     [all are shown]                 {messy for growth}
$ elab a format-fasta ...   [modular organized]             {odd}
$ elab a.format-fasta ...   [all are shown, but organized]  {ugly}
```


### TEMPLATE WORKFLOW
```bash
# one shot analysis w/ iterations
$ twig init -n test1 -w ... --min-bit 300 --data data/*.gz 
$ twig seqlab -j test1.json -c cores -t threads
$ twig branch -j test1.json -n test2 --add x --drop y --min-bit 300
$ twig seqlab -j test2.json -c cores -t threads
$ twig treelab -j test1.json -c 10 -t 2 

$ twig run -d ... -o ... -it 2

# step by step
$ twig init -d /data/*.fa -o /tmp/test -i 1.5 --min-bit 30 -n test ... 
$ twig infer-orthogroups -j /tmp/test/test.json 
$ twig infer-trees -j ...
$ twig infer-events -j ...
$ twig refine-orthogroups -j ... -it 2

# continuation w/ parameter or sample changes
$ twig branch -j ... -n ... -a /data/new*.fa -r A B C 

# view stats and filepaths
$ twig summary -j ...
```

# CONSIDER
- blastp search should find even low bitscore hits; we can filter heavily later.
- lastz as alternative mapper?
- ...


### idea: many small tools and some big tools workflow.
```bash
twig run -d ... -w ... --min --max ...               # creates a dir and runs all.

# ...
twig seqlab                                          # creates project dir with all outputs
twig seqlab.score -d data/*.fa --minx .. --maxx ...  # writes blast files to a specified dir
twig seqlab.graph -b blasts/*.match --min --max      # writes graph with filtered edges from blast pairs
twig seqlab.split -g name.graph -i 1.4 -o ...        # writes orthogroups tsv
twig seqlab.build -g name.groups -d data/*.fa --opts # writes orthogroup sequence alignments
twig seqlab.refine ...
twig seqlab.cds --prot data/*.fa --cds data2/*.fa --macs2 ... # align cds with prots.

# ...
twig treelab -a alignments/*.fa --boots --...           # create project dir with all outputs and stats
twig treelab.infer -a alignments/*.fa --raxml --boots       # infers trees, weights, ...
twig treelab.filter -t trees/*.nwk --min-support --relabel  # collapse; trim long edges; re-infer after removal; 
twig treelab.msc -t ftrees/*.nwk --astral-args --gene-tree-rooting-args # infer sptree, gene tree rooting, sptree rooting.
twig treelab.dlc -t ftrees/*.nwk --dlc-args                   # label gene trees w/ dups, losses, and write stats.
twig treelab.metadata -t ftrees/*.nwk --dlc ... --sptree ...  # label gene trees with species tree nodes, 

# ...
twig evolab run -j proj.json --param1 --param2              # run many evol analyses into a project dir.
twig evolab.dlc --
twig evolab.csubst --cds aligned_cds/... --trees ftrees/... --sptree tree.nwk  # ...
twig evolab.diffexpr --cds ... --counts ... --args ...      # 
```


# STEPS
[x] Load protein fasta inputs into a Project JSON schema.  
[ ] Validate protein fasta files: require unique headers.  
[x] Perform all-by-all diamond search  
[x] Parse search TSVs to matrix of [gene x gene] bit-scores * length  
[x] Construct graph from scores  
[x] Dissect graph into homo-clusters and find best inflation parameter  
[x] Write cluster sequences to files.
[x] Align sequences with mafft.
[x] Trim sequences alignments with clipkit.
[x] Infer gene trees by raxml w/ or w/o supports.
[ ] Infer species tree by astral-pro.
[ ] Infer rooting of species tree by stride-like method -> write tree
[ ] Infer rooting of gene trees  
[ ] Count dups/losses/HGTs -> write trees w/ NHX
[ ] Label gene trees with tip names and internal DLC and root prob labels
[ ] Write stats summaries.
[ ] Revise graph and re-run.
[ ] automate removal of long branches / poor alignment

# Continuing a seqlab analysis
[ ] cli 
[x] load json
[ ] add to JSON a simple step-completed section
[ ] cleaning fasta headers log all completed or run remaining
[ ] creating blast database log all completed or run remaining
[ ] pairwise blast searches log all completed for run remaining
[ ] rerun graph to clusters only if any above changed
[ ] if graph changed then run all alignments, else only unaligned

# TO ADD SAMPLES REQUIRES
[ ] cli --branch command to add or drop | use --add List[str] --drop List[~str]
[ ] load existing json
[ ] revise step-completed in JSON
[ ] seqlab using step-completed and detection of files:
    [ ] create relabeled faa
    [ ] create blast db
    [ ] pairwise blast to other relabeled faas
    [ ] write graph, split clusters, align and trim
    - In summary, the only thing reused is blast DBs and relabel faas
      both of which are super fast to recreate. The existing pairwise
      blast hits are the main thing 

# LATER
[ ] Allow adding new samples to an existing dataset.  
[ ] Allow toggling more options.  