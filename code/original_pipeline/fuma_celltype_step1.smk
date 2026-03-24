import os

rule all:
    input:
        expand(os.path.join(config["out_dir"], "{trait}/{scrna}/{level}/place_holder.txt"), trait=config["traits"], scrna=config["scrnas"], level=config["levels"])
        
rule magma_step1: #magma_mean_linear
    input:
        gene_results = os.path.join(config["gwas_dir"], config["gene_region"] + "_" + config["gene_name"],  "{trait}_" + config["gene_region"] + "_" + config["gene_name"] + "_magma.genes.raw")
    output:
        os.path.join(config["out_dir"], "{trait}/{scrna}/{level}/place_holder.txt")
    params:
        magma_executable = config["magma_executable"],
        out_prefix = config["out_dir"] + "/{trait}/{scrna}/{level}/magma",
        gene_covar = os.path.join(config["scrna_dir"], "{scrna}/{level}/means_cell_log_counts_pM.tsv")
    run: 
        shell("touch {output}");
        shell("{params.magma_executable} --gene-results {input.gene_results} --gene-covar {params.gene_covar} --model condition-hide=Average direction=greater --out {params.out_prefix}");
        shell("mv {params.out_prefix}.gsa.out {params.out_prefix}.gsa.out.tmp"); 
        shell("grep -v '#' {params.out_prefix}.gsa.out.tmp > {params.out_prefix}.gsa.out"); 
        shell("rm -r {params.out_prefix}.gsa.out.tmp"); 

