# PVOGs functional associations

A proof-of-concept, automated and reproducible pipeline for predicting 
functional associations between
[pVOGs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5210652/).

Pre-print: TBD

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.14.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

## TL;DR

```
# Clone this repo
$ git clone https://github.com/MGXlab/pvogs_function.git

# Get in here
$ cd pvogs_function

# Create the environment with conda
$ conda env create -n pvogs --file=environment.yml
$ conda activate pvogs

# Dry run 
(pvogs)$ snakemake --use-conda -j 16 -np

# If dry run finished successfully
(pvogs)$ snakemake --use-conda -j 16 --conda-frontend mamba
```

## Description

The main purpose of this repository is to host the code necessary for full 
reproducibility.

* Raw data required are [hosted on zenodo](https://zenodo.org/record/4576599). 
These are automatically downloaded when executing the workflow, 
so no need to get them.
Alterantively, you can download the archive from zenodo and unpack it in a 
directory called `data` in here.

Most of the steps are not easily configurable, unless you take a dive into all
the rules and scripts. This is by choice.

## Requirements

* A working [conda](https://docs.conda.io/en/latest/) installation
* `snakemake >= 5.14` (developed with `5.24.1`)
  * Optional: `mamba == 0.5.1` (speeds up environment dependency resolution 
  and creation)

You can create the same `conda` environment used during development with the 
provided `environment.yml`.

```
$ conda env create -n pvogs --file=./environment.yml
```

Make sure you activate it before you launch snakemake
```
$ conda activate pvogs
(pvogs)$ snakemake --version
5.24.1
```

## Configuration

The configuration options are included in the `config/config.yaml`. This file 
must be present in the `config` directory and named `config.yaml`.

Options in that file:

- `negatives`: Specifies the number of negative datasets to create. 
   10 is used in the manuscript.

 > Changing this will most likely break the workflow

- the Zenodo `doi` 
  
  All necessary input is provided as a gzipped tape archive (`.tar.gz`) 
  [available on Zenodo](https://zenodo.org/record/4576599) with record id 
  `4576599`. This record id must be used if you wish to get the same data. 
  
- `threads` per rule

  For the most resource demanding rules included in the config, you can 
  specify the number of cores each rule will utilize at runtime. These are set 
  to reasonable values for the setup used during development (`Ubuntu 
  20.04.2 LTS` with `120G` of RAM and `20` processors) for a reasonable 
  parallelization/runtime balance. 
  **You should adjust these according to your own local setup.**


## Usage

This workflow was built and tested locally. No cluster execution was tested,
although it should be feasible to scale it with
[the appropriate profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). 

All commands assume you have `cd`ed in the root directory of this repo and you
are using a conda environment created with the provided `environment.yml`. Here
it is named `pvogs`, and the prefix `(pvogs)` in front of the prompt (`$`) 
shows it is activated.

1. **Dry run** - _Always a good idea before launching the whole worfklow_

```
(pvogs)$ snakemake --use-conda -j 16 -np
```

2. If the dry run completed with no errors,

* You can execute the workflow by removing the `-n` flag. 
* Adjust the number of parallel jobs (`-j`) according to your setup
* Remove the `-p` flag if you don't want the commands to be printed.

```
(pvogs)$ snakemake --use-conda -j16 -p
```
2a. Speed up environment creation with mamba

If `mamba` is available in your snakemake environment, or if you created a new
environment with the `environment.yml` provided here:

```
(pvogs)$ snakemake --use-conda -j16 --conda-frontend mamba
```

## Jupyter integration

A central [notebook](workflow/notebooks/analysis.py.ipynb) is used for all 
visualization and machine learning (model search) purposes. Its main output is
the `results/RF/best_model.pkl` file. This file stores the final model, which 
is then picked up by the rest of the workflow to make predictions on the 
full target set and produce final output.

In the workflow, this comes after all pre-processing has been done, 
i.e. all feature values have been calculated, all negative datasets have been 
generated. That means you can either:

- Execute the whole workflow once in its default, non-interactive mode
and (hopefully) get the same output (yay for reproducibility). Then re-run
the workflow with the option for interactive data exploration on.

- Enable the interaction with the notebook already from the start. You need
to produce the target output of the rule (`results/RF/best_model_id.txt`).
Once the `results/RF/best_model.pkl` is written you can save the changes, 
and quit the server
([more info here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration) 
and you can always 
[watch this demo](https://snakemake.readthedocs.io/en/stable/_images/snakemake-notebook-demo.gif)).
This will trigger the execution of the rest of the workflow.

In both cases the command you are looking for is

```
(pvogs)$ snakemake --use-conda -j16 --edit-notebook random_forest
```

Changes are stored in the original notebook.

Another jupyter notebook is available and used for all annotation statistics
reported, as well as the generation of some supplementary info. The main output
of this is `results/figures/fig_6.pdf`. For successful execution you need to 
create this file.

Similarly, you can run this one with: 

```
(pvogs)$ snakemake --use-conda -j16 --edit-notebook annotation_stats
```

**Note: Depending on the changes you make, the results you get will 
differ from the default, non-interactive run and what is reported on the 
paper**.


## Output

The output of the whole workflow is produced and stored within a `results` 
directory, with the structure shown below. (several directories and files have 
been omitted).

The most prominent ones are marked with a short description:

```
# Skipping a lot of intermediate files with the -I option
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
        ├── positives.search_stats.tsv        * Supplementary table 1.
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
│   │   ├── fig_1.[pdf|svg|eps|png]
        ....
....
│   ├── metrics.pkl
│   ├── metrics.stats.tsv ------------------- * Mean. max, min, std across all models
│   ├── metrics.tsv ------------------------- * Exact values of metrics for each combination of training/validation
│   └── models
│       ├── N10.RF.pkl ---------------------- * Best model obtained when optimizing with each negative set
        .....
.....		
│── scores.tsv  ----------------------------- * Table with feature values for all possible pVOGs combinations
│── predictions_annotations_features.tsv ---- * Master table that contains all results. (Supplementary table 2)
```

