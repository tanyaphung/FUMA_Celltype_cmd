import os

# this is a snakemake file to run modified FUMA cell type step 2 and 3 specifiaclly for Phung et al. 202x. 
# for other applications, please make sure to modified the snakemake file accordingly

traits = ["gwas_1", "gwas_2", "gwas_3", "gwas_4", "gwas_5", "gwas_6", "gwas_7", "gwas_8", "gwas_9", "gwas_10", "gwas_11", "gwas_12", "gwas_13", "gwas_14"]

# for prenatal
base_dir="/projects/gwas_celltype_atlas/out/prenatal"
regions = ["brain", "forebrain", "midbrain", "hindbrain", "transientStructuresOfForebrain", "telencephalon", "diencephalon", "cerebralGyriAndLobules", "cerebellum", "pons", "medulla", "neocortex", "allocortex", "periallocortex", "frontalNeocortex", "parietalNeocortex", "temporalNeocortex", "cingulateNeocortex", "prefrontalCortex", "dorsolateralPrefrontalCortex", "primaryMotorCortex", "orbitalFrontalCortex", "primarySomatosensoryCortex", "primaryVisualCortex", "meninges"] #prenatal

# # for postnatal
# base_dir="/gpfs/work5/0/vusr0480/projects/gwas_celltype_atlas/out/postnatal"
# regions = ["midbrain", "cerebralGyriAndLobules", "middleTemporalGyrus", "whiteMatter", "myelencephalon", "cerebralNuclei", "hypothalamus", "thalamus", "cerebellum", "pons", "neocortex", "allocortex", "periallocortex", "hippocampalGyrusFormation", "frontalNeocortex", "parietalNeocortex", "temporalNeocortex", "occipitalNeocortex", "cingulateNeocortex", "insularNeocortex", "vagalNucleus", "prefrontalCortex", "dorsolateralPrefrontalCortex", "ventrolateralPrefrontalCortex", "primaryMotorCortex", "orbitalFrontalCortex", "primarySomatosensoryCortex", "primaryAuditoryCortex", "primaryVisualCortex", "spinalCord"] #postnatal

rule all:
    input:
        expand("{base_dir}/{trait}/{region}/magma_celltype_step1.txt", base_dir=base_dir, trait=traits, region=regions),
        expand("{base_dir}/{trait}/all_regions/scrna_ds.txt", base_dir=base_dir, trait=traits),
        expand("{base_dir}/{trait}/all_regions/celltype/place_holder.txt", base_dir=base_dir, trait=traits),
        expand("{base_dir}/{trait}/all_regions/place_holder_2.txt", base_dir=base_dir, trait=traits),
        expand("{base_dir}/{trait}/all_regions/magma_celltype_step3_filtered.txt", base_dir=base_dir, trait=traits),
        expand("{base_dir}/{trait}/all_regions/magma_celltype_step3_filtered_sig.txt", base_dir=base_dir, trait=traits)

rule post_step1_sig:
    input:
        ds_list = "{base_dir}/{trait}/{region}/scrna_ds.txt"
    output:
        "{base_dir}/{trait}/{region}/magma_celltype_step1.txt"
    params:
        base_dir = "{base_dir}",
        filedir = "{base_dir}/{trait}/{region}/",
        magmadir = "/gpfs/home6/tphung/software/magma",
        magmafiles = "{base_dir}/{trait}/{region}/",
        trait = "{trait}"
    run:
        shell("Rscript fuma_celltype_postStep1.R --filedir {params.filedir} --adjPmeth bonferroni --magmadir {params.magmadir} --magmafiles {params.magmafiles} --ds_list {input.ds_list}");

        shell("cp {params.filedir}/celltype/* {params.base_dir}/{params.trait}/all_regions/celltype/");

rule combine_regions:
    output:
        "{base_dir}/{trait}/all_regions/scrna_ds.txt"
    params:
        trait = "{trait}",
        base_dir = "{base_dir}"
    run:
        shell("touch {output}");
        shell("python combine_step1_all_regions.py --trait {params.trait} --base_dir {params.base_dir}");

rule prepare_parent_ds:
    output:
        "{base_dir}/{trait}/all_regions/celltype/place_holder.txt"
    params:
        celltype_dir = "{base_dir}/{trait}/all_regions/celltype",
        ds_list = "{base_dir}/{trait}/all_regions/scrna_ds.txt"
    run:
        shell("touch {output}");
        shell("python create_parent_ds_input.py --ds_magma_dir {params.celltype_dir} --parent_magma_dir {params.celltype_dir} --ds_list {params.ds_list}");

# IMPORTANT NOTES: after the rule prepare_parent_ds, run these following commands on the command line: 

# run both for prenatal and postnatal
# for i in {1..14}; do
# cp gwas_${i}/all_regions/magma_celltype_step1.txt gwas_${i}/all_regions/magma_celltype_step1.txt.bk
# sed 's/[(|)]//g' gwas_${i}/all_regions/magma_celltype_step1.txt.bk > gwas_${i}/all_regions/magma_celltype_step1.txt
# sed -i 's/[(|)]//g' gwas_${i}/all_regions/celltype/Sepp2023.txt
# sed -i 's/[(|)]//g' gwas_${i}/all_regions/celltype/OteroGarcia.txt
# sed -i 's/[(|)]//g' gwas_${i}/all_regions/celltype/TorresFlores.txt
# done;

# for i in {1..14}; do
# cp /magma/out/protein_coding_ensemble/gwas_${i}_protein_coding_ensemble_magma.genes.raw gwas_${i}/all_regions/magma.genes.raw
# done;

rule modified_step23:
    input:
    output:
        "{base_dir}/{trait}/all_regions/place_holder_2.txt"
    params:
        filedir = "{base_dir}/{trait}/all_regions/",
        magmadir = "/gpfs/home6/tphung/software/magma",
        ds_list = "{base_dir}/{trait}/all_regions/scrna_ds.txt"
    run:
        shell("touch {output}");
        shell("Rscript fuma_celltype_step2and3_mod_all_regions.R --filedir {params.filedir} --adjPmeth bonferroni --magmadir {params.magmadir} --magmafiles {params.filedir} --ds_list {params.ds_list}");

rule filter_step23:
    input:
    output:
        "{base_dir}/{trait}/all_regions/magma_celltype_step3_filtered.txt"
    params:
        step1_2 = "{base_dir}/{trait}/all_regions/step1_2_summary.txt",
        step3 = "{base_dir}/{trait}/all_regions/magma_celltype_step3.txt"
    run:
        shell("Rscript filter_step3_all_regions.R --step1_2_sum_fp {params.step1_2} --step3_fp {params.step3} --out_fp {output}");

rule remove_drop:
    input:
    output:
        "{base_dir}/{trait}/all_regions/magma_celltype_step3_filtered_sig.txt"
    params:
        "{base_dir}/{trait}/all_regions/magma_celltype_step3_filtered.txt"
    run:
        shell("grep -v drop {params} > {output}");