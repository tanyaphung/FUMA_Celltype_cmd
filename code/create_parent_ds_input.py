import pandas as pd
import os
from collections import defaultdict
import argparse

# This script combines the columns from the same parent study to use in modified step2 where we are doing within parent study conditional analysis.

parser = argparse.ArgumentParser()
parser.add_argument('--ds_magma_dir', required=True) #example: /mnt/d/reference_data/MAGMA/celltype
parser.add_argument('--parent_magma_dir', required=True) #example: /home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/fuma_632461/celltype/
parser.add_argument('--ds_list', required=True)
args = parser.parse_args()

ds_magma_dir = args.ds_magma_dir
parent_magma_dir = args.parent_magma_dir

datasets = []
with open(args.ds_list, "r") as f:
    for line in f:
        datasets.append(line.rstrip("\n"))

parent_ds = defaultdict(list)

for dataset in datasets:
    parent_ds[dataset.split("_")[1]].append(dataset)
    

for parent in parent_ds:
    parent_df = pd.DataFrame()
    for dataset in parent_ds[parent]:
        id = dataset.split("_")[0]
        data_df = pd.read_csv(os.path.join(ds_magma_dir, dataset + ".txt"), sep="\t")
        data_df = data_df.drop('Average', axis=1)
        colnames = data_df.columns.values.tolist()
        for colname in colnames:
            if colname=="GENE":
                continue
            new_colname = id + "_" + colname
            data_df.rename({colname:new_colname}, axis=1, inplace=True)
            
        if len(parent_df) == 0:
            parent_df = data_df
        else:
            parent_df = pd.merge(parent_df, data_df, on="GENE")
            
    parent_df["Average"] = parent_df.iloc[:,1:].mean(axis=1)
    parent_df.to_csv(os.path.join(parent_magma_dir, parent + ".txt"), sep="\t", index=False)
            
