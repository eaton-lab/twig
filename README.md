
# TWIG

## IN DEVELOPMENT (use with caution)

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

### Modules
- [a] FA/multiFA/GFF manipulation and matching.
- [o] ORTHOLOGY + PHYLO inference
- [s] SYNTENY analysis
- [p] CDS/FA prep [align/trim/frame/mask/SFS/...]
- [e] CDS/FA analysis [iMK,Convergence,]


# MICROSYNTENY PHYLOGENETIC ANALYSIS
https://www.nature.com/articles/s41467-021-23665-0

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
