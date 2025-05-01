import os

rule all:
    input:
        expand("celltype/{trait}/{region}/place_holder.txt", trait=config["traits"], region=config["regions"])
        
rule magma_step2_3:
    output:
        "celltype/{trait}/{region}/place_holder.txt"
    params:
        filedir = "celltype/{trait}/{region}/",
        ids = "{region}_scrna_ids.txt",
        magmadir = "MAGMA",
        magmafiles = "magma_covs"
    run: 
        shell("touch {output}");
        shell("Rscript fuma_celltype_step2and3.R --filedir {params.filedir} --ids_file {params.ids} --magmadir {params.magmadir} --magmafiles {params.magmafiles}");