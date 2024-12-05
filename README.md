# Installation
First, install conda and nextflow on your system
```
git clone --recurse-submodules https://github.com/dmckeow/crg-bcaortho.git
```
# Usage
```
# Prepare your input fasta file directory by making symbolic soft links (ln -s) to your protein fastas (one file per genome) and place them all within a single directory
# Running locally with example data provided
nextflow run main.nf -profile local

# Running using SLURM with example data provided
sbatch submit_nf.sh main.nf -profile slurm
```
## Changing parameters, resources, etc
The following files control the paramters for the pipeline:
    * nextflow.config - parameters unlikely to be changed. If you change anything in this it will be the memory and cpu paramaters
        Recommendation - create new profiles in this file (see local, slurm for examples) and load them from the command line to change memory and cpus (see below)
    * params.config - parameters that are likely to be changed frequently, such as arguments for software, etc. It is loaded by the nextflow.config by default
        Recommendation - copy this file with a memorable name and change the parameters within, then provide to command line with -c
You can run the pipeline with a different config.file e.g.:
```
sbatch submit_nf.sh main.nf -profile slurm -c params-fullhmms.config

```

# Notes
### How Broccoli was incorporated in the pipeline

```
# Broccoli added as a git submodule
git submodule add https://github.com/rderelle/Broccoli.git broccoli

# Broccoli submodule committed
git add .gitmodules broccoli
git commit -m "Add Broccoli as a submodule"

# To update the Broccoli submodule
cd broccoli
git pull origin master
cd ..
git add broccoli
git commit -m "Update Broccoli submodule"


# If the broccoli module mysteriously disappears
git config --file=.gitmodules --get-regexp path # check if the submodule still setup
git submodule init
git submodule update
```
### How much resources?
## OrthoFinder

## Broccoli
8 threads & 0.1-0.5 Gb per eukaryotic genome (from Broccoli documentation)


# Running with test dataset

* Target gene families:
    * Homeodomains
    * Ets
    * TFs
    * zf-met
    * zf-c2h2

* Target genomes:
    * **List of genomes to include:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/species_list_annotation.txt
    * **Genomes:** /users/asebe/xgraubove/genomes/data/*long.pep.fasta

* HMM profies:
    * From Grisha (has an extra field for domain name):
        * **gene_family_info:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/results_phylogenies/gene_families_searchinfo.csv
        * **hmm_dir:** /users/asebe/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/results_phylogenies/hmms
    * From Xavi (larger and newer dataset):
        * **gene_family_info:** /users/asebe/xgraubove/climate-adaptation-corals-ops/results_annotation/data/gene_families_searchinfo.csv
        * **hmm_dir:** /users/asebe/xgraubove/climate-adaptation-corals-ops/results_annotation/data/hmms


# ORTHOFINDER PROBLEM

## sbatch submit_nf.sh main.nf -profile slurm -params-file params-MetazoansTest_Homeodomains-part1_3.json  Failed
* parts 1 and 3 success separately, but together FAIL
* Problem genomes:
    * 41: Coeast_long.pep.None.Homeodomains.domains.fasta
    * 90: Sros_long.pep.None.Homeodomains.domains.fasta
* Errors:
    * WARNING: Too few hits between species 90 and species 41 to normalise the scores, these hits will be ignored
    * WARNING: program called by OrthoFinder produced output to stderr

        Command: mcl /users/asebe/dmckeown/projects/crg-bcaortho/work/95/84bcd09d04906d600a2b2417c5b0b7/input/OrthoFinder/Results_orthofinder/WorkingDirectory/OrthoFinder_graph.txt -I 1.5 -o /users/asebe/dmckeown/projects/crg-bcaortho/work/95/84bcd09d04906d600a2b2417c5b0b7/input/OrthoFinder/Results_orthofinder/WorkingDirectory/clusters_OrthoFinder_I1.5.txt -te 2 -V all

        stdout
        ------
        b''
        stderr
        ------
        b'[mcl] added <2> garbage entries\n'

    * Analysing Orthogroups
        =====================

        Calculating gene distances
        --------------------------
        2024-12-05 08:37:49 : Done
        /usr/local/lib/python3.12/site-packages/numpy/core/fromnumeric.py:3504: RuntimeWarning: Mean of empty slice.
        return _methods._mean(a, axis=axis, dtype=dtype,
        /usr/local/lib/python3.12/site-packages/numpy/core/_methods.py:129: RuntimeWarning: invalid value encountered in scalar divide
        ret = ret.dtype.type(ret / rcount)

    * ERROR: Species tree inference failed
        cannot access local variable 'speciesTree' where it is not associated with a value

## sbatch submit_nf.sh main.nf -profile slurm -params-file params-MetazoansTest_Homeodomains-part2-2.json  Failed
## sbatch submit_nf.sh main.nf -profile slurm -params-file params-MetazoansTest_Homeodomains-part2-3.json  Failed
* Part 2 failed so split it into three parts - parts 2 and 3 failed
* Errors:
    * Analysing Orthogroups
        =====================

        Calculating gene distances
        --------------------------
        2024-12-05 08:16:30 : Done
        /usr/local/lib/python3.12/site-packages/numpy/core/fromnumeric.py:3504: RuntimeWarning: Mean of empty slice.
        return _methods._mean(a, axis=axis, dtype=dtype,
        /usr/local/lib/python3.12/site-packages/numpy/core/_methods.py:129: RuntimeWarning: invalid value encountered in scalar divide
        ret = ret.dtype.type(ret / rcount)
        2024-12-05 08:16:30 : Done 0 of 62
        2024-12-05 08:16:30 : Done 10 of 62
        2024-12-05 08:16:30 : Done 20 of 62
        2024-12-05 08:16:30 : Done 30 of 62
        2024-12-05 08:16:30 : Done 40 of 62
        2024-12-05 08:16:31 : Done 50 of 62
        Using fallback species tree inference method

        Inferring gene and species trees
        --------------------------------

        Best outgroup(s) for species tree
        ---------------------------------
        2024-12-05 08:16:36 : Starting STRIDE
        File "/usr/local/bin/scripts_of/__main__.py", line 1790, in main
            GetOrthologues(speciesInfoObj, options, prog_caller)
        File "/usr/local/bin/scripts_of/__main__.py", line 1550, in GetOrthologues
            orthologues.OrthologuesWorkflow(speciesInfoObj.speciesToUse,
        File "/usr/local/bin/scripts_of/orthologues.py", line 1047, in OrthologuesWorkflow
            roots, clusters_counter, rootedSpeciesTreeFN, nSupport, _, _, stride_dups = stride.GetRoot(spTreeFN_ids, files.FileHandler.GetOGsTreeDir(), stride.GeneToSpecies_dash, nHighParallel, qWriteRootedTree=True)
                                                                                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        File "/usr/local/bin/scripts_of/stride.py", line 513, in GetRoot
            species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
                                                                ^^^^^^^^^^^
        ERROR: Species tree inference failed
        cannot access local variable 'speciesTree' where it is not associated with a value
        Traceback (most recent call last):
        File "/usr/local/bin/orthofinder", line 7, in <module>
            main(args)
        File "/usr/local/bin/scripts_of/__main__.py", line 1790, in main
            GetOrthologues(speciesInfoObj, options, prog_caller)
        File "/usr/local/bin/scripts_of/__main__.py", line 1550, in GetOrthologues
            orthologues.OrthologuesWorkflow(speciesInfoObj.speciesToUse,
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        File "/usr/local/bin/scripts_of/orthologues.py", line 1047, in OrthologuesWorkflow
            roots, clusters_counter, rootedSpeciesTreeFN, nSupport, _, _, stride_dups = stride.GetRoot(spTreeFN_ids, files.FileHandler.GetOGsTreeDir(), stride.GeneToSpecies_dash, nHighParallel, qWriteRootedTree=True)
                                                                                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        File "/usr/local/bin/scripts_of/stride.py", line 513, in GetRoot
            species, dict_clades, clade_names = AnalyseSpeciesTree(speciesTree)
                                                                ^^^^^^^^^^^
        UnboundLocalError: cannot access local variable 'speciesTree' where it is not associated with a value

# Remove no overlap genome combinations from orthofinder
## If a genome with no overlaps with anohter genome is identified, we want to start over with Orthofinder with those genomes removed from the input channel
input/OrthoFinder/Results_orthofinder/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv
