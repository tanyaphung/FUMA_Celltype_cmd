# Date: 2025-07-15
# The goal of this script is to combine the file magma_celltype_step1.txt for all regions

import argparse
import os
import shutil


parser = argparse.ArgumentParser()
parser.add_argument('--trait', required=True)
parser.add_argument('--base_dir', required=True) #/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow
args = parser.parse_args()

regions = ["forebrain", "midbrain", "hindbrain", "telencephalon", "diencephalon", "myelencephalon", "neocortex", "allocortex", "periallocortex", "pons", "cerebellum", "medulla", "thalamus", "hypothalamus", "frontalNeocortex", "parietalNeocortex", "occipitalNeocortex", "cingulateNeocortex", "temporalNeocortex", "insularNeocortex", "prefrontalCortex", "dorsolateralPrefrontalCortex", "orbitalFrontalCortex", "ventrolateralPrefrontalCortex", "middleTemporalGyrus", "primaryMotorCortex", "primarySomatosensoryCortex", "primaryAuditoryCortex", "primaryVisualCortex", "transientStructuresOfForebrain", "cerebralGyriandLobules", "whiteMatter", "cerebralNuclei", "vagalNucleus", "hippocampalGyrusFormation", "spinalCord"]
# regions = ["forebrain", "primaryVisualCortex"]

trait = args.trait
base_dir = args.base_dir

outfile = open(os.path.join(base_dir, trait, "all_regions/magma_celltype_step1.txt"), "w")

header = ["Dataset", "Cell_type", "NGENES", "BETA", "BETA_STD", "SE", "P", "P.adj.pds", "P.adj"]
print("\t".join(header), file=outfile)

scrna_ds_set = set()

for i in regions: 
    if os.path.exists(os.path.join(base_dir, trait, i, "magma_celltype_step1_sig.txt")):
        with open(os.path.join(base_dir, trait, i, "magma_celltype_step1_sig.txt"), "r") as f:
            for line in f:
                if line.startswith("Dataset"):
                    continue
                print(line.rstrip("\n"), file=outfile)
                ds = line.rstrip("\n").split("\t")[0]
                # tmp = i + "/" + line.rstrip("\n").split("\t")[0]
                scrna_ds_set.add(ds)
                
                ds_path = os.path.join(base_dir, trait, "all_regions", "celltype", ds + ".txt")
                if not os.path.exists(ds_path):
                    orig_path = os.path.join(base_dir, trait, i, "celltype", ds + ".txt")
                    dest = shutil.copyfile(orig_path, ds_path)

outfile.close()

scrna_ds = open(os.path.join(base_dir, trait, "all_regions/scrna_ds.txt"), "w")
for i in scrna_ds_set:
    print(i, file=scrna_ds)
scrna_ds.close()