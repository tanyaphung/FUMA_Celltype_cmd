# FUMA_Celltype_cmd
- This repo hosts codes for running FUMA Cell Type on the command line
- The script `code/fuma_celltype_step2and3.R` is based off of the script https://github.com/vufuma/FUMA-webapp/blob/master/scripts/magma_celltype/magma_celltype.R (credit: Kyoko Watanabe). 
    - The original script was written to be run as part of FUMA Celltype.
    - I adapted the original script to make it into a standalone script that can be run locally on the command line, instead of having to submit to FUMA

- I set up a snakemake pipeline to run the Rscript: `code/fuma_celltype_step2_3.smk`. Examples to run: 
```
snakemake -s code/fuma_celltype_step2_3.smk --configfile code/config_fuma_ct_step2_3.json -j 
```

# Modified workflow
- Modified FUMA cell type that was used in Phung et al. 202x
![alt text](FUMA_Celltype_mod.png)

- Note that in order to adapt these following scripts to new analyses, one need to modify variables and paths, etc... Please contact Tanya if you need help with this. 
- Step 1: 
    - snakemake: `code/fuma_celltype_step1.smk`
    - configfile: `/home/tnphung/FUMA-dev/FUMA_Celltype_cmd/code/config_fuma_ct_step1.json`
    - run command: `snakemake -s fuma_celltype_step1.smk --rerun-incomplete -j -k --configfile config_fuma_ct_step1.json`
- Step 2 and 3: 
    - snakemake: `code/fuma_celltype_step2and3_mod_all_regions.R`
    - run command: `snakemake -s modified_fuma_celltype_step2_3.smk --rerun-incomplete -j -k`