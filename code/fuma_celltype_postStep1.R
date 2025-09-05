library(data.table)
library(kimisc)
library(dplyr)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = '')

# Add arguments
parser$add_argument("--filedir", type = "character", required = TRUE,
                    help = "Path")
parser$add_argument("--adjPmeth", type = "character", required = TRUE,
                    help = "adjPmeth")
parser$add_argument("--magmadir", type = "character", required = TRUE,
                    help = "magmadir")
parser$add_argument("--magmafiles", type = "character", required = TRUE,
                    help = "magmafiles")
parser$add_argument("--ds_list", type = "character", required = TRUE,
                    help = "ds_list")

# Parse arguments
args <- parser$parse_args()


filedir = args$filedir #/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/transientStructuresOfForebrain/
adjPmeth = args$adjPmeth #bonferroni
magmadir = args$magmadir #/home/tnphung/FUMA-dev/FUMA_Celltype_cmd/code
magmafiles = args$magmafiles #/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/transientStructuresOfForebrain
ds_list = args$ds_list #/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/transientStructuresOfForebrain/scrna_ds.txt

datasets_df = read.table(ds_list)
colnames(datasets_df) = "datasets"
datasets = datasets_df$datasets

unique_ds_count = c()

step1 <- data.frame()
for(ds in datasets){
    f <- paste0(filedir, "magma_celltype_", ds, ".gsa.out")
    if (file.exists(f)) {
	tmp <- fread(cmd=paste0("grep -v '^#' ", filedir, "magma_celltype_", ds, ".gsa.out"), data.table=F)
	if("FULL_NAME" %in% colnames(tmp)){
		tmp$VARIABLE <- tmp$FULL_NAME #convert VARIABLE to FULL_NAME
		tmp <- tmp[,-ncol(tmp)] #remove the FULL_NAME column
	}
    unique_ds_count = c(unique_ds_count, tmp$VARIABLE)
    tmp$VARIABLE <- paste0(strsplit(ds, "_")[[1]][1], "_", tmp$VARIABLE) #update variable to append the data set id at the front
	tmp <- tmp[order(tmp$P),] #order by p values
	tmp$ds <- ds #add ds column
	tmp$P.adj.pds <- p.adjust(tmp$P, method=adjPmeth) #add P.adj.pds column 
	if(nrow(step1)==0){step1 <- tmp}
	else{step1 <- rbind(step1, tmp)}
}
}

step1$P.adj <- p.adjust(step1$P, method=adjPmeth) #add P.adj column 
tmp_out <- step1[,c("ds", "VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "P.adj.pds", "P.adj")]
colnames(tmp_out)[1:2] <- c("Dataset", "Cell_type")
write.table(tmp_out, paste0(filedir, "magma_celltype_step1.txt"), quote=F, row.names=F, sep="\t")
# rm(tmp_out)
print(length(unique(unique_ds_count)))
unique_ds_count_len = length(unique(unique_ds_count)) #TODO: fix this because the VARIABLE is now unique because I appended the dataset id
tmp_out = tmp_out[which(tmp_out$P<(0.05/unique_ds_count_len)),] #modified to do bonferroni correction for the number of unique cell types
print(nrow(tmp_out))
if (nrow(tmp_out) == 0) {print("There is no significant cell type after step 1")}

write.table(tmp_out, paste0(filedir, "magma_celltype_step1_sig.txt"), quote=F, row.names=F, sep="\t")