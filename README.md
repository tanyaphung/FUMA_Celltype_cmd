# FUMA_Celltype_cmd
- This repo hosts codes for running FUMA Cell Type on the command line
- The script `code/fuma_celltype_step2and3.R` is based off of the script https://github.com/vufuma/FUMA-webapp/blob/master/scripts/magma_celltype/magma_celltype.R (credit: Kyoko Watanabe). 
    - The original script was written to be run as part of FUMA Celltype.
    - I adapted the original script to make it into a standalone script that can be run locally on the command line, instead of having to submit to FUMA

- I set up a snakemake pipeline to run the Rscript: `code/fuma_celltype_step2_3.smk`. Examples to run: 
```
snakemake -s code/fuma_celltype_step2_3.smk --configfile code/config_fuma_ct_step2_3.json -j 
```