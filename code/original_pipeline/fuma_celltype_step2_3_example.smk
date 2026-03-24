import os

rule all:
    input:
        expand("{base_dir}/out/fuma_ct/{trait}/{region}/place_holder.txt", trait=config["traits"], region=config["regions"], base_dir=config["base_dir"])
        
rule magma_step2_3:
    output:
        "{base_dir}/out/fuma_ct/{trait}/{region}/place_holder.txt"
    params:
        filedir = "{base_dir}/out/fuma_ct/{trait}/{region}/",
        ids = "{region}_scrna_ids.txt",
        magmadir = "magma",
        magmafiles = "{base_dir}/magma_covs"
    run: 
        shell("touch {output}");
        shell("Rscript fuma_celltype_step2and3.R --filedir {params.filedir} --ids_file {params.ids} --magmadir {params.magmadir} --magmafiles {params.magmafiles}");