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


# Propagating sample names across init_ortho

## SEARCH

python bin/orthogroup_summary.py --input results/example_data_fullhmms_SEARCH/deflines/deflines_combined.txt --orthofinder results/example_data_fullhmms_SEARCH/orthofinder/Orthogroups/Orthogroups.tsv --broccoli results/example_data_fullhmms_SEARCH/broccoli/dir_step3/table_OGs_protein_names.txt --diamond_mcl results/example_data_fullhmms_SEARCH/dmnd_mcl/dmnd_mcl.orthogroups.txt --mmseqs results/example_data_fullhmms_SEARCH/mmseqs/mmseqs.orthogroups.txt --search

(tabs between OG and genome column, space delimited within genome columns)

results/example_data_fullhmms_SEARCH/deflines/deflines_combined.txt
sample  fasta   defline
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P000887:32-802
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P000943:1-720
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P001526:1-695
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P002149:1-147
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P002687:12-782
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P003048:1-508
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P003207:1-242
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P003434:1-631
Pygbif  Pygbif_long.pep.None.Myosin.domains.fasta       Pygbif_EP00033_Pygsuia_biforma_P003707:46-188



## without SEARCH

(tabs between OG and genome column, space delimited within genome columns)

results/example_data_fullhmms_ALL/deflines/deflines_combined.txt
sample  fasta   defline
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000001
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000002
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000003
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000004
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000005
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000006
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000007
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000008
Pygbif  Pygbif_long.pep.fasta   Pygbif_EP00033_Pygsuia_biforma_P000009


## OrthoFinder
(tabs between OG and genome column, comma + space (, ) delimited within genome columns)

results/example_data_fullhmms_SEARCH/orthofinder/Orthogroups/Orthogroups.tsv
Orthogroup      Otau_long.pep.None.Myosin.domains       Pygbif_long.pep.None.Myosin.domains     Salkve_long.pep.None.Myosin.domains     Xesp_long.pep.None.Myosin.domains
OG0000000               Pygbif_EP00033_Pygsuia_biforma_P002149_1-147, Pygbif_EP00033_Pygsuia_biforma_P003207_1-242, Pygbif_EP00033_Pygsuia_biforma_P003707_46-188, Pygbif_EP00033_Pygsuia_biforma_P005362_36-359, Pygbif_EP00033_Pygsuia_biforma_P005776_16-178, Pygbif_EP00033_Pygsuia_biforma_P007503_1-99, Pygbif_EP00033_Pygsuia_biforma_P008049_1-146    Salkve_EP00050_Salpingoeca_kvevrii_P001348_5-151, Salkve_EP00050_Salpingoeca_kvevrii_P001349_1-140, Salkve_EP00050_Salpingoeca_kvevrii_P003763_36-275, Salkve_EP00050_Salpingoeca_kvevrii_P014413_1-402, Salkve_EP00050_Salpingoeca_kvevrii_P020062_11-158, Salkve_EP00050_Salpingoeca_kvevrii_P021405_8-154        Xesp_009839-T1_1-467, Xesp_009998-T1_1-115, Xesp_010000-T1_1-135, Xesp_017515-T1_31-215, Xesp_017536-T1_31-211
OG0000001       Otau_XM_022984064.1_33-805      Pygbif_EP00033_Pygsuia_biforma_P002687_12-782, Pygbif_EP00033_Pygsuia_biforma_P003048_1-508, Pygbif_EP00033_Pygsuia_biforma_P004222_1-263, Pygbif_EP00033_Pygsuia_biforma_P004321_16-779, Pygbif_EP00033_Pygsuia_biforma_P007950_1-185        Salkve_EP00050_Salpingoeca_kvevrii_P001267_18-778, Salkve_EP00050_Salpingoeca_kvevrii_P003767_1-515, Salkve_EP00050_Salpingoeca_kvevrii_P010140_29-782, Salkve_EP00050_Salpingoeca_kvevrii_P020058_1-196, Salkve_EP00050_Salpingoeca_kvevrii_P021406_1-298    Xesp_001364-T1_16-774, Xesp_014913-T1_30-804, Xesp_020905-T1_19-774
OG0000002               Pygbif_EP00033_Pygsuia_biforma_P000887_32-802, Pygbif_EP00033_Pygsuia_biforma_P005501_37-643, Pygbif_EP00033_Pygsuia_biforma_P008030_1-196      Salkve_EP00050_Salpingoeca_kvevrii_P001266_40-819, Salkve_EP00050_Salpingoeca_kvevrii_P017667_1-325, Salkve_EP00050_Salpingoeca_kvevrii_P017671_1-417 Xesp_004923-T1_432-1252, Xesp_016506-T1_33-806, Xesp_017509-T1_31-807, Xesp_017520-T1_120-471
OG0000003                       Salkve_EP00050_Salpingoeca_kvevrii_P001260_40-836, Salkve_EP00050_Salpingoeca_kvevrii_P001264_1-759, Salkve_EP00050_Salpingoeca_kvevrii_P001265_1-745, Salkve_EP00050_Salpingoeca_kvevrii_P001270_1-759, Salkve_EP00050_Salpingoeca_kvevrii_P001273_1-762, Salkve_EP00050_Salpingoeca_kvevrii_P001275_1-735, Salkve_EP00050_Salpingoeca_kvevrii_P001276_1-752, Salkve_EP00050_Salpingoeca_kvevrii_P001277_18-799, Salkve_EP00050_Salpingoeca_kvevrii_P001279_1-761, Salkve_EP00050_Salpingoeca_kvevrii_P001281_1-762
OG0000004               Pygbif_EP00033_Pygsuia_biforma_P000943_1-720, Pygbif_EP00033_Pygsuia_biforma_P003434_1-631, Pygbif_EP00033_Pygsuia_biforma_P004619_1-334, Pygbif_EP00033_Pygsuia_biforma_P006409_1-125Salkve_EP00050_Salpingoeca_kvevrii_P001344_1-331        Xesp_004732-T1_18-759, Xesp_009997-T1_1-666
OG0000005               Pygbif_EP00033_Pygsuia_biforma_P001526_1-695, Pygbif_EP00033_Pygsuia_biforma_P004568_1-556      Salkve_EP00050_Salpingoeca_kvevrii_P006115_1-744        Xesp_003099-T1_1-667, Xesp_003393-T1_1-769
OG0000006               Pygbif_EP00033_Pygsuia_biforma_P005558_15-459, Pygbif_EP00033_Pygsuia_biforma_P008300_1-297     Salkve_EP00050_Salpingoeca_kvevrii_P020061_1-287        Xesp_016418-T1_53-994
OG0000007                       Salkve_EP00050_Salpingoeca_kvevrii_P001274_1-343, Salkve_EP00050_Salpingoeca_kvevrii_P001280_1-341, Salkve_EP00050_Salpingoeca_kvevrii_P001343_1-419, Salkve_EP00050_Salpingoeca_kvevrii_P016899_1-318
OG0000008                       Salkve_EP00050_Salpingoeca_kvevrii_P001258_1-782, Salkve_EP00050_Salpingoeca_kvevrii_P001261_1-769, Salkve_EP00050_Salpingoeca_kvevrii_P001278_1-767


## Broccoli:

(tabs between OG and genome column, space delimited within genome columns)

results/example_data_fullhmms_SEARCH/broccoli/dir_step3/table_OGs_protein_names.txt
#OG_name        Salkve_long.pep.None.Myosin.domains.fasta       Pygbif_long.pep.None.Myosin.domains.fasta       Xesp_long.pep.None.Myosin.domains.fasta Otau_long.pep.None.Myosin.domains.fasta
OG_1    Salkve_EP00050_Salpingoeca_kvevrii_P017667:1-325 Salkve_EP00050_Salpingoeca_kvevrii_P017671:1-417 Salkve_EP00050_Salpingoeca_kvevrii_P000069:49-801     Pygbif_EP00033_Pygsuia_biforma_P005501:37-643 Pygbif_EP00033_Pygsuia_biforma_P000887:32-802   Xesp_017520-T1:120-471 Xesp_017515-T1:31-215 Xesp_017536-T1:31-211 Xesp_017509-T1:31-807 Xesp_016506-T1:33-806
OG_2    Salkve_EP00050_Salpingoeca_kvevrii_P001256:1-357 Salkve_EP00050_Salpingoeca_kvevrii_P020533:1-620 Salkve_EP00050_Salpingoeca_kvevrii_P001343:1-419 Salkve_EP00050_Salpingoeca_kvevrii_P021347:1-388 Salkve_EP00050_Salpingoeca_kvevrii_P001353:354-874 Salkve_EP00050_Salpingoeca_kvevrii_P001352:354-796 Salkve_EP00050_Salpingoeca_kvevrii_P017772:1-144 Salkve_EP00050_Salpingoeca_kvevrii_P018593:1-143 Salkve_EP00050_Salpingoeca_kvevrii_P010140:29-782 Salkve_EP00050_Salpingoeca_kvevrii_P016793:1-177 Salkve_EP00050_Salpingoeca_kvevrii_P021348:1-387 Salkve_EP00050_Salpingoeca_kvevrii_P001272:1-154 Salkve_EP00050_Salpingoeca_kvevrii_P001862:1-687 Salkve_EP00050_Salpingoeca_kvevrii_P018595:37-690 Salkve_EP00050_Salpingoeca_kvevrii_P001260:40-836 Salkve_EP00050_Salpingoeca_kvevrii_P001264:1-759 Salkve_EP00050_Salpingoeca_kvevrii_P001265:1-745 Salkve_EP00050_Salpingoeca_kvevrii_P001266:40-819 Salkve_EP00050_Salpingoeca_kvevrii_P001267:18-778 Salkve_EP00050_Salpingoeca_kvevrii_P001270:1-759 Salkve_EP00050_Salpingoeca_kvevrii_P001273:1-762 Salkve_EP00050_Salpingoeca_kvevrii_P001275:1-735 Salkve_EP00050_Salpingoeca_kvevrii_P001276:1-752 Salkve_EP00050_Salpingoeca_kvevrii_P001277:18-799 Salkve_EP00050_Salpingoeca_kvevrii_P001279:1-761 Salkve_EP00050_Salpingoeca_kvevrii_P001281:1-762 Salkve_EP00050_Salpingoeca_kvevrii_P016794:1-590 Salkve_EP00050_Salpingoeca_kvevrii_P001258:1-782 Salkve_EP00050_Salpingoeca_kvevrii_P001261:1-769 Salkve_EP00050_Salpingoeca_kvevrii_P001274:1-343 Salkve_EP00050_Salpingoeca_kvevrii_P001278:1-767 Salkve_EP00050_Salpingoeca_kvevrii_P001280:1-341 Salkve_EP00050_Salpingoeca_kvevrii_P017773:278-464 Salkve_EP00050_Salpingoeca_kvevrii_P008789:351-1182 Salkve_EP00050_Salpingoeca_kvevrii_P001268:1-754       Pygbif_EP00033_Pygsuia_biforma_P005776:16-178 Pygbif_EP00033_Pygsuia_biforma_P005362:36-359 Pygbif_EP00033_Pygsuia_biforma_P004321:16-779 Pygbif_EP00033_Pygsuia_biforma_P002687:12-782       Xesp_010949-T1:113-988 Xesp_014913-T1:30-804 Xesp_001364-T1:16-774
OG_3    Salkve_EP00050_Salpingoeca_kvevrii_P001348:5-151 Salkve_EP00050_Salpingoeca_kvevrii_P001349:1-140 Salkve_EP00050_Salpingoeca_kvevrii_P014275:1-341 Salkve_EP00050_Salpingoeca_kvevrii_P001344:1-331 Salkve_EP00050_Salpingoeca_kvevrii_P001346:1-327 Salkve_EP00050_Salpingoeca_kvevrii_P014414:1-347 Salkve_EP00050_Salpingoeca_kvevrii_P014413:1-402 Salkve_EP00050_Salpingoeca_kvevrii_P021224:1-328 Salkve_EP00050_Salpingoeca_kvevrii_P021407:1-325 Salkve_EP00050_Salpingoeca_kvevrii_P021405:8-154 Salkve_EP00050_Salpingoeca_kvevrii_P001347:1-298 Salkve_EP00050_Salpingoeca_kvevrii_P021406:1-298 Salkve_EP00050_Salpingoeca_kvevrii_P014277:1-410      Pygbif_EP00033_Pygsuia_biforma_P004222:1-263 Pygbif_EP00033_Pygsuia_biforma_P004619:1-334 Pygbif_EP00033_Pygsuia_biforma_P002149:1-147 Pygbif_EP00033_Pygsuia_biforma_P003048:1-508 Pygbif_EP00033_Pygsuia_biforma_P000943:1-720 Pygbif_EP00033_Pygsuia_biforma_P003434:1-631 Pygbif_EP00033_Pygsuia_biforma_P008457:1-150    Xesp_009839-T1:1-467 Xesp_004732-T1:18-759 Xesp_009997-T1:1-666 Xesp_009998-T1:1-115 Xesp_010000-T1:1-135
OG_4    Salkve_EP00050_Salpingoeca_kvevrii_P021223:1-692 Salkve_EP00050_Salpingoeca_kvevrii_P021221:12-779 Salkve_EP00050_Salpingoeca_kvevrii_P012733:14-784 Salkve_EP00050_Salpingoeca_kvevrii_P000845:1-421 Salkve_EP00050_Salpingoeca_kvevrii_P006115:1-744 Salkve_EP00050_Salpingoeca_kvevrii_P003767:1-515       Pygbif_EP00033_Pygsuia_biforma_P007447:1-548 Pygbif_EP00033_Pygsuia_biforma_P003825:34-904 Pygbif_EP00033_Pygsuia_biforma_P004568:1-556 Pygbif_EP00033_Pygsuia_biforma_P001526:1-695 Pygbif_EP00033_Pygsuia_biforma_P003207:1-242 Pygbif_EP00033_Pygsuia_biforma_P003707:46-188 Pygbif_EP00033_Pygsuia_biforma_P004908:1-444 Pygbif_EP00033_Pygsuia_biforma_P005558:15-459  Xesp_003393-T1:1-769 Xesp_003099-T1:1-667 Xesp_011648-T1:10-756 Xesp_020905-T1:19-774   Otau_XM_022984064.1:33-805
OG_5    Salkve_EP00050_Salpingoeca_kvevrii_P022611:16-775               Xesp_004923-T1:432-1252


## DIAMOND MCL

(tabs between OG and genome column, space delimited within genome columns)

results/example_data_fullhmms_SEARCH/dmnd_mcl/combined_input.orthogroups.txt
#placeholder    header
combined_input.dmnd.csv.1.mcl   Pygbif_EP00033_Pygsuia_biforma_P000887:32-802 Xesp_017509-T1:31-807 Xesp_016506-T1:33-806 Salkve_EP00050_Salpingoeca_kvevrii_P001266:40-819 Pygbif_EP00033_Pygsuia_biforma_P005501:37-643 Otau_XM_022984064.1:33-805 Salkve_EP00050_Salpingoeca_kvevrii_P001267:18-778 Xesp_020905-T1:19-774 Pygbif_EP00033_Pygsuia_biforma_P002687:12-782 Pygbif_EP00033_Pygsuia_biforma_P004321:16-779 Pygbif_EP00033_Pygsuia_biforma_P000943:1-720 Pygbif_EP00033_Pygsuia_biforma_P003434:1-631 Xesp_004732-T1:18-759 Pygbif_EP00033_Pygsuia_biforma_P001526:1-695 Salkve_EP00050_Salpingoeca_kvevrii_P006115:1-744 Xesp_009997-T1:1-666 Xesp_003099-T1:1-667 Xesp_003393-T1:1-769 Pygbif_EP00033_Pygsuia_biforma_P003048:1-508 Pygbif_EP00033_Pygsuia_biforma_P004568:1-556 Pygbif_EP00033_Pygsuia_biforma_P002149:1-147 Salkve_EP00050_Salpingoeca_kvevrii_P014413:1-402 Salkve_EP00050_Salpingoeca_kvevrii_P014277:1-410 Xesp_009839-T1:1-467 Xesp_014913-T1:30-804 Salkve_EP00050_Salpingoeca_kvevrii_P010140:29-782 Xesp_001364-T1:16-774 Salkve_EP00050_Salpingoeca_kvevrii_P002988:1-682 Salkve_EP00050_Salpingoeca_kvevrii_P001277:18-799 Salkve_EP00050_Salpingoeca_kvevrii_P021223:1-692 Pygbif_EP00033_Pygsuia_biforma_P003207:1-242 Salkve_EP00050_Salpingoeca_kvevrii_P000845:1-421 Pygbif_EP00033_Pygsuia_biforma_P003707:46-188 Salkve_EP00050_Salpingoeca_kvevrii_P001349:1-140 Salkve_EP00050_Salpingoeca_kvevrii_P001348:5-151 Salkve_EP00050_Salpingoeca_kvevrii_P001350:43-189 Xesp_010000-T1:1-135 Xesp_009998-T1:1-115 Pygbif_EP00033_Pygsuia_biforma_P003825:34-904 Salkve_EP00050_Salpingoeca_kvevrii_P012733:14-784 Salkve_EP00050_Salpingoeca_kvevrii_P021221:12-779 Xesp_011648-T1:10-756 Pygbif_EP00033_Pygsuia_biforma_P004222:1-263 Salkve_EP00050_Salpingoeca_kvevrii_P001352:354-796 Salkve_EP00050_Salpingoeca_kvevrii_P021406:1-298 Pygbif_EP00033_Pygsuia_biforma_P004619:1-334 Salkve_EP00050_Salpingoeca_kvevrii_P001344:1-331 Salkve_EP00050_Salpingoeca_kvevrii_P021407:1-325 Salkve_EP00050_Salpingoeca_kvevrii_P021224:1-328 Salkve_EP00050_Salpingoeca_kvevrii_P014414:1-347 Salkve_EP00050_Salpingoeca_kvevrii_P001346:1-327 Pygbif_EP00033_Pygsuia_biforma_P004908:1-444 Pygbif_EP00033_Pygsuia_biforma_P004961:1-148 Salkve_EP00050_Salpingoeca_kvevrii_P000846:1-340 Pygbif_EP00033_Pygsuia_biforma_P005362:36-359 Pygbif_EP00033_Pygsuia_biforma_P005558:15-459 Pygbif_EP00033_Pygsuia_biforma_P005776:16-178 Pygbif_EP00033_Pygsuia_biforma_P006409:1-125 Salkve_EP00050_Salpingoeca_kvevrii_P000847:1-242 Pygbif_EP00033_Pygsuia_biforma_P007447:1-548 Salkve_EP00050_Salpingoeca_kvevrii_P003767:1-515 Pygbif_EP00033_Pygsuia_biforma_P007503:1-99 Salkve_EP00050_Salpingoeca_kvevrii_P018595:37-690 Salkve_EP00050_Salpingoeca_kvevrii_P008789:351-1182 Pygbif_EP00033_Pygsuia_biforma_P007745:1-125 Pygbif_EP00033_Pygsuia_biforma_P009341:1-159 Salkve_EP00050_Salpingoeca_kvevrii_P001278:1-767 Salkve_EP00050_Salpingoeca_kvevrii_P001261:1-769 Salkve_EP00050_Salpingoeca_kvevrii_P022611:16-775 Pygbif_EP00033_Pygsuia_biforma_P007950:1-185 Pygbif_EP00033_Pygsuia_biforma_P007985:1-232 Salkve_EP00050_Salpingoeca_kvevrii_P001265:1-745 Salkve_EP00050_Salpingoeca_kvevrii_P001276:1-752 Salkve_EP00050_Salpingoeca_kvevrii_P001264:1-759 Salkve_EP00050_Salpingoeca_kvevrii_P001281:1-762 Salkve_EP00050_Salpingoeca_kvevrii_P001260:40-836 Pygbif_EP00033_Pygsuia_biforma_P008030:1-196 Pygbif_EP00033_Pygsuia_biforma_P008049:1-146 Salkve_EP00050_Salpingoeca_kvevrii_P017671:1-417 Pygbif_EP00033_Pygsuia_biforma_P008300:1-297 Pygbif_EP00033_Pygsuia_biforma_P008383:1-155 Salkve_EP00050_Salpingoeca_kvevrii_P001862:1-687 Pygbif_EP00033_Pygsuia_biforma_P008457:1-150 Salkve_EP00050_Salpingoeca_kvevrii_P017667:1-325 Salkve_EP00050_Salpingoeca_kvevrii_P020058:1-196 Xesp_010949-T1:113-988 Xesp_002868-T1:1-735 Salkve_EP00050_Salpingoeca_kvevrii_P001270:1-759 Xesp_004923-T1:432-1252 Salkve_EP00050_Salpingoeca_kvevrii_P003763:36-275 Salkve_EP00050_Salpingoeca_kvevrii_P001275:1-735 Xesp_016418-T1:53-994 Salkve_EP00050_Salpingoeca_kvevrii_P000069:49-801 Xesp_017520-T1:120-471 Xesp_017515-T1:31-215 Xesp_017536-T1:31-211 Salkve_EP00050_Salpingoeca_kvevrii_P001256:1-357 Salkve_EP00050_Salpingoeca_kvevrii_P001280:1-341 Salkve_EP00050_Salpingoeca_kvevrii_P001279:1-761 Salkve_EP00050_Salpingoeca_kvevrii_P020533:1-620 Salkve_EP00050_Salpingoeca_kvevrii_P001343:1-419 Salkve_EP00050_Salpingoeca_kvevrii_P001258:1-782 Salkve_EP00050_Salpingoeca_kvevrii_P001273:1-762 Salkve_EP00050_Salpingoeca_kvevrii_P001268:1-754 Salkve_EP00050_Salpingoeca_kvevrii_P001272:1-154 Salkve_EP00050_Salpingoeca_kvevrii_P021348:1-387 Salkve_EP00050_Salpingoeca_kvevrii_P001274:1-343 Salkve_EP00050_Salpingoeca_kvevrii_P016794:1-590 Salkve_EP00050_Salpingoeca_kvevrii_P016899:1-318 Salkve_EP00050_Salpingoeca_kvevrii_P021347:1-388 Salkve_EP00050_Salpingoeca_kvevrii_P014275:1-341 Salkve_EP00050_Salpingoeca_kvevrii_P001347:1-298 Salkve_EP00050_Salpingoeca_kvevrii_P020062:11-158 Salkve_EP00050_Salpingoeca_kvevrii_P021405:8-154 Salkve_EP00050_Salpingoeca_kvevrii_P001353:354-874 Salkve_EP00050_Salpingoeca_kvevrii_P001861:10-157 Salkve_EP00050_Salpingoeca_kvevrii_P017773:278-464 Salkve_EP00050_Salpingoeca_kvevrii_P016793:1-177 Salkve_EP00050_Salpingoeca_kvevrii_P005345:1-679 Salkve_EP00050_Salpingoeca_kvevrii_P009558:1-296 Salkve_EP00050_Salpingoeca_kvevrii_P016414:311-717 Salkve_EP00050_Salpingoeca_kvevrii_P009898:1-54 Salkve_EP00050_Salpingoeca_kvevrii_P016320:1-112 Salkve_EP00050_Salpingoeca_kvevrii_P016900:1-211 Salkve_EP00050_Salpingoeca_kvevrii_P017772:1-144 Salkve_EP00050_Salpingoeca_kvevrii_P020059:1-146 Salkve_EP00050_Salpingoeca_kvevrii_P017774:1-471 Salkve_EP00050_Salpingoeca_kvevrii_P018593:1-143 Salkve_EP00050_Salpingoeca_kvevrii_P020061:1-287 Salkve_EP00050_Salpingoeca_kvevrii_P020531:1-80 Salkve_EP00050_Salpingoeca_kvevrii_P020532:87-853 Salkve_EP00050_Salpingoeca_kvevrii_P021218:1-742
combined_input.dmnd.csv.2.mcl   Salkve_EP00050_Salpingoeca_kvevrii_P013164:1-103


## MMSEQS2

(tabs between OG and genome column, space delimited within genome columns)

results/example_data_fullhmms_SEARCH/mmseqs/clusters.orthogroups.txt
#placeholder    header
clusters.tsv.1.mmseqs   Salkve_EP00050_Salpingoeca_kvevrii_P021405:8-154
clusters.tsv.2.mmseqs   Salkve_EP00050_Salpingoeca_kvevrii_P013164:1-103
clusters.tsv.3.mmseqs   Pygbif_EP00033_Pygsuia_biforma_P002149:1-147 Salkve_EP00050_Salpingoeca_kvevrii_P001272:1-154
clusters.tsv.4.mmseqs   Salkve_EP00050_Salpingoeca_kvevrii_P014413:1-402 Salkve_EP00050_Salpingoeca_kvevrii_P001352:354-796 Salkve_EP00050_Salpingoeca_kvevrii_P014277:1-410 Pygbif_EP00033_Pygsuia_biforma_P005558:15-459 Salkve_EP00050_Salpingoeca_kvevrii_P021348:1-387 Salkve_EP00050_Salpingoeca_kvevrii_P001353:354-874 Salkve_EP00050_Salpingoeca_kvevrii_P017671:1-417 Salkve_EP00050_Salpingoeca_kvevrii_P000845:1-421
clusters.tsv.5.mmseqs   Pygbif_EP00033_Pygsuia_biforma_P003207:1-242
clusters.tsv.6.mmseqs   Xesp_009839-T1:1-467
clusters.tsv.7.mmseqs   Salkve_EP00050_Salpingoeca_kvevrii_P016320:1-112
clusters.tsv.8.mmseqs   Xesp_009998-T1:1-115
clusters.tsv.9.mmseqs   Salkve_EP00050_Salpingoeca_kvevrii_P016414:311-717

