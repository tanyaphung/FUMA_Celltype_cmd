# FUMA_Celltype_cmd
- This repo hosts codes for running FUMA Cell Type on the command line. 
**Disclaimer:** use as own risk. I do not guarantee that the results will match 100% with the outputs from FUMA. 
- The script `code/fuma_celltype_step2and3.R` is based off of the script https://github.com/vufuma/FUMA-webapp/blob/master/scripts/magma_celltype/magma_celltype.R (credit: Kyoko Watanabe). 
    - The original script was written to be run as part of FUMA Celltype.
    - I adapted the original script to make it into a standalone script that can be run locally on the command line, instead of having to submit to FUMA

- I set up a snakemake pipeline to run the Rscript: `code/fuma_celltype_step2_3.smk`. Examples to run: 
```
snakemake -s code/fuma_celltype_step2_3.smk --configfile code/config_fuma_ct_step2_3.json -j 
```

## Instructions to run FUMA Cell type original workflow
### Step 1: 
- snakemake file for running step 1 is: `code/fuma_celltype_step1.smk`
- the configfile for step 1 is: `code/config_fuma_ct_step1.json`. Note the following: 
    - Line 2: give the full path to magma executable
    - Line 3 and 4 can be left as is if you processed the scRNAseq datasets with `protein_coding` and `ensemble`. If not, update as appropriate. See https://github.com/tanyaphung/scrnaseq_viewer/tree/main/code/postprocessing for how to process the scRNAseq data 
    - Line 5: Content of the `gwas_dir` directory `/Processed_GWAS/magma/out/` contains a folder called `protein_coding_ensemble`. In it there are outputs from MAGMA. 
    ```
    ── protein_coding_ensemble
        ├── ald_protein_coding_ensemble_magma.genes.raw
    ```
    - Content of the `scrna_dir` directory `/Processed_scRNA/data/magma/` contains a folder for each scRNAseq. For example: 
    ```
    ├── 70_Siletti_Hypothalamus.HTHma.MN_Human_2022
    │   ├── cell_type_level_1
    │   │   ├── means_cell_log_counts_pM.tsv
    │   │   ├── place_holder.txt
    │   │   └── spec_cell_log_counts_pM.tsv
    │   ├── cell_type_level_2
    │   │   ├── means_cell_log_counts_pM.tsv
    │   │   ├── place_holder.txt
    │   │   └── spec_cell_log_counts_pM.tsv
    │   └── cell_type_level_3
    │       └── place_holder.txt
    ```
- run command: `snakemake -s fuma_celltype_step1.smk --rerun-incomplete -j5 --configfile config_fuma_ct_step1.json`
    - modify the value after the `-j` flag as appropriate for your data

- The outputs are stored in the folder specified on line 6 of the config file. In this example, it is in `out/fuma_ct`: 
```
└── ald
    ├── 70_Siletti_Hypothalamus.HTHma.MN_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 71_Siletti_Hypothalamus.HTHpo.HTHso_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 72_Siletti_Hypothalamus.HTHma.HTHtub_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 73_Siletti_Hypothalamus.HTHso.HTHtub_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 74_Siletti_Hypothalamus.HTHpo_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 75_Siletti_Hypothalamus.HTHtub_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    ├── 76_Siletti_Hypothalamus.HTHso_Human_2022
    │   └── cell_type_level_2
    │       ├── magma.gsa.out
    │       ├── magma.log
    │       └── place_holder.txt
    └── 77_Siletti_Hypothalamus.HTHma_Human_2022
        └── cell_type_level_2
            ├── magma.gsa.out
            ├── magma.log
            └── place_holder.txt
```

### Step 2 and 3: 
- The snakemake file is currently set up so that it will run step 2 and 3 for all of the scRNAseq specified in a file, which is `{region}_scrna_ids.txt`
- Create files `{region}_scrna_ids.txt` where region corresponds to the region of interest. In this example, the `region` is specified in the config file `code/config_fuma_ct_step2_3.json` is Hypothalamus. So create a file called `Hypothalamus_scrna_ids.txt`. Contents of the file should include the IDs and the cell type level 
```
cat Hypothalamus_scrna_ids.txt
ids
70_Siletti_Hypothalamus.HTHma.MN_Human_2022_level_2
71_Siletti_Hypothalamus.HTHpo.HTHso_Human_2022_level_2
72_Siletti_Hypothalamus.HTHma.HTHtub_Human_2022_level_2
73_Siletti_Hypothalamus.HTHso.HTHtub_Human_2022_level_2
74_Siletti_Hypothalamus.HTHpo_Human_2022_level_2
75_Siletti_Hypothalamus.HTHtub_Human_2022_level_2
76_Siletti_Hypothalamus.HTHso_Human_2022_level_2
77_Siletti_Hypothalamus.HTHma_Human_2022_level_2
```

- Create a directory `/out/fuma_ct/ald/Hypothalamus` 
    - Move the outputs from step 1 to appropriate folder
    - Make sure that in this directory, the files are named `magma_celltype_{ds}.gsa.out` where `ds` is the dataset in the file `{region}_scrna_ids.txt`. Example: 
    ```
    pwd
    /out/fuma_ct/ald/Hypothalamus

    for i in 70_Siletti_Hypothalamus.HTHma.MN_Human_2022 71_Siletti_Hypothalamus.HTHpo.HTHso_Human_2022 72_Siletti_Hypothalamus.HTHma.HTHtub_Human_2022 73_Siletti_Hypothalamus.HTHso.HTHtub_Human_2022 74_Siletti_Hypothalamus.HTHpo_Human_2022 75_Siletti_Hypothalamus.HTHtub_Human_2022 76_Siletti_Hypothalamus.HTHso_Human_2022 77_Siletti_Hypothalamus.HTHma_Human_2022; do
    cp ../${i}/cell_type_level_2/magma.gsa.out magma_celltype_${i}_level_2.gsa.out
    done;
    ```

- Create a directory `/magma_covs/celltype`. The files in this directory are for example: `14_Siletti_CerebralCortex.SMG_Human_2022_level_2.txt`. Example:

```
pwd
/magma_covs/celltype

for i in 70_Siletti_Hypothalamus.HTHma.MN_Human_2022 71_Siletti_Hypothalamus.HTHpo.HTHso_Human_2022 72_Siletti_Hypothalamus.HTHma.HTHtub_Human_2022 73_Siletti_Hypothalamus.HTHso.HTHtub_Human_2022 74_Siletti_Hypothalamus.HTHpo_Human_2022 75_Siletti_Hypothalamus.HTHtub_Human_2022 76_Siletti_Hypothalamus.HTHso_Human_2022 77_Siletti_Hypothalamus.HTHma_Human_2022; do
cp ../../Processed_scRNA/data/magma/${i}/cell_type_level_2/means_cell_log_counts_pM.tsv ${i}_level_2.txt
done;
```
    - Note that the reason why this was set up this way was because the file structures for the processed scRNAseq data was set up first before implementing step2 and 3. TODO in the future is to make the file structure easier/less redundant. 

- Make sure that the file `Processed_GWAS/magma/out/protein_coding_ensemble/ald_protein_coding_ensemble_magma.genes.raw` is copied over to the folder `/out/fuma_ct/ald/Hypothalamus` as `magma.genes.raw`

- Make sure to put in the correct path for magma on line 13 of `code/original_pipeline/fuma_celltype_step2_3_example.smk`

- Run step 2 and 3: 
    - Copy and modify the snakemake file: 
    ```
    cp fuma_celltype_step2_3_example.smk fuma_celltype_step2_3.smk
    ```
    
    ```
    snakemake -s fuma_celltype_step2_3.smk -j1 --configfile config_fuma_ct_step2_3.json --rerun-incomplete
    ```

# Modified workflow
- Modified FUMA cell type that was used in Phung et al. 202x
![alt text](FUMA_Celltype_mod.png)

- Note that in order to adapt these following scripts to new analyses, one need to modify variables and paths, etc... Please contact Tanya if you need help with this. 
- Step 1: 
    - snakemake: `code/fuma_celltype_step1.smk`
    - configfile: `code/config_fuma_ct_step1.json`
    - run command: `snakemake -s fuma_celltype_step1.smk --rerun-incomplete -j -k --configfile config_fuma_ct_step1.json`
- Step 2 and 3: 
    - snakemake: `code/fuma_celltype_step2and3_mod_all_regions.R`
    - run command: `snakemake -s modified_fuma_celltype_step2_3.smk --rerun-incomplete -j -k`