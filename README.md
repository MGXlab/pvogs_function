# PVOGs functions interactions


## TL;DR


```
# Clone this repo
$ git clone https://git.science.uu.nl/n.pappas/pvogs_function.git

# Get in there
$ cd pvogs_function

# Optional, if snakemake>=5.14 and conda available
$ conda env create -n my_env --file=environment.yml
$ conda activate my_env


# Dry run to check that it works
(my_env)$ snakemake --use-conda -n
```

## Description
The main purpose of this repository is to host the code necessary for full reproducibility.

* Raw data required are hosted on [zenodo sandbox](https://sandbox.zenodo.org/record/666719#.X1c5qoZS_J8). These are automatically
downloaded when executing the workflow, so no need to get them.

Most of the steps are not easily configurable, unless you take a dive into all the rules and scripts. This is by choice.

## Requirements

* A working [conda](https://docs.conda.io/en/latest/) installation
* `snakemake >= 5.14` (any version with jupyter integration should do)
  * Optional: `mamba == 0.5.1` (speeds up environment dependency resolution and creation)

You can create the same `conda` environment used during development with the provided `environment.yml`.
```
$ conda env create -n my_env --file=./environment.yml
```

Make sure you activate it before you launch snakemake
```
$ conda activate my_env
(my_env)$ snakemake --version
5.23.0
```

## Configuration

The configuration options are included in the `config/config.yml`.
These include:
- `negatives`: Specifies the number of negative datasets to create. 10 is used in the manuscript.
  >Changing this will most likely break the workflow

- the zenodo dois Until the workflow gets published, I am using the zenodo sandbox for testing.
- `threads` per rule
  For the most resource demanding rules included in the config, you can specify the number of cores
  each rule will utilize at runtime. I have set these to reasonable values for my own local
  setup (`Ubuntu 16.04.1 x86_64` with `120Gb` of RAM and `20` processors) for a good
  parallelization/runtime balance. **You should adjust these according to your own local
  setup.**

## Usage
---
Currently, this workflow was built and tested on a local machine with an X server available 
(i.e. you can do stuff in a GUI).

>If you run this on a remote machine, make sure that you (can) ssh with `ssh -X ...`.
>This is required for the `summarize_intact.py` script, that uses the `ete3` package
>to do some plotting.

### **Option 1. This repo**

`cd` into the root directory of this repo.

- Dry run

Always a good idea before launching the whole worfklow
```
$ snakemake --use-conda -j16 -np
```

If the dry run completed with no errors you can execute the worfklow by removing the `-n` flag. 
* Adjust the number of parallel jobs (`-j`) according to your setup
* Remove the `-p` flag if you don't want the commands to be printed.
```
$ snakemake --use-conda -j16 -p
```
- Speed up environment creation with mamba

If `mamba` is available in your snakemake environment, or if you created a new environment with the `environment.yml`
provided here:
```
$ snakemake --use-conda -j16 --conda-frontend mamba
```

- Jupyter integration

A central notebook is used for all visualization and machine learning (model search) purposes.
Its main output is the `results/RF/best_model.pkl` file.

If you want to fiddle around with it yourself
```
$ snakemake --use-conda -j16 --conda-frontend mamba --edit-notebook results/RF/best_model.pkl
```
Once the `results/RF/best_model.pkl` is written you can save the changes, and quit the server
([more info here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration) and
you can always [see this demo](https://snakemake.readthedocs.io/en/stable/_images/snakemake-notebook-demo.gif)).
This will trigger the execution of the rest of the workflow.

The resulting notebook will be saved as `results/logs/processed_notebook.py.ipynb`.

**Note that depending on the changes you make, the results you might get will differ from the default, non-interactive run**.


### Option 2. Archived workflow from zenodo (TO DO).

Something along the [guidelines from snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#sustainable-and-reproducible-archiving).


## Output

The output of the whole workflow is produced and stored within a `results` directory. 
This has the structure shown below.
(several directories and files have been omitted)
The most prominent ones are marked with a short description:
```
# Skipping several thousands of intermediate files with the -I option
$ tree -n -I '*NC*.fasta|*_genes.*|*.gff|*.log' results

results
├── annotations.tsv
├── filtered_scores.tsv -------------------- * Table containing feature values for all interactions passing filtering
├── final_training_set.tsv
├── interaction_datasets
│   ├── 01_filter_intact
│   ├── 02_summarize_intact
│   ├── 03_uniprot
│   ├── 04_process_uniprot
│   ├── 05_genomes
│   ├── 06_map_proteins_to_pvogs
│   ├── N1  --------------------------------  
....                                        | * Features, interactions, proteins, and pvogs are stored per dataset
│   └── positives --------------------------  
│       ├── positives.features.tsv
│       ├── positives.interactions.tsv
│       ├── positives.proteins.faa
│       └── positives.pvogs_interactions.tsv
├── logs
├── predictions.tsv ------------------------- * Final predictions made
├── pre_process
│   ├── all_genomes
│   ├── comparem  --------------------------- * Directory with the final AAI matrix used
...
│   ├── fastani  ---------------------------- * Directory with the final ANI matrix used
│   ├── hmmsearch  -------------------------- * HMMER search results for all pvogs profiles agains the translated genomes
│   ├── reflist.txt
│   └── transeq
│       └── transeq.genomes.fasta
├── RF
│   ├── best_model_id.txt ------------------- * Contains the id of the negative dataset
│   ├── best_model.pkl ---------------------- * The best model obtained.
│   ├── features_stats.tsv ------------------ * Mean, max, min. std for feature importances
│   ├── features.tsv ------------------------ * Exact values of features importances for each combination of training/validation
│   ├── figures ----------------------------- * Figures used in the manuscript.       
│   │   ├── Figure_1a.svg
        ....
....
│   ├── metrics.pkl
│   ├── metrics.stats.tsv ------------------- * Mean. max, min, std across all models
│   ├── metrics.tsv ------------------------- * Exact values of metrics for each combination of training/validation
│   └── models
│       ├── N10.RF.pkl ---------------------- * Best model obtained when optimizing with each negative set
        .....
.....		
└── scores.tsv  ----------------------------- * Master table with feature values for all possible pVOGs combinations

```

