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


filedir = args$filedir
adjPmeth = args$adjPmeth
magmadir = args$magmadir
magmafiles = args$magmafiles
ds_list = args$ds_list

# filedir = "/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/all_regions/"
# adjPmeth = "bonferroni"
# magmadir = "/home/tnphung/FUMA-dev/FUMA_Celltype_cmd/code"
# magmafiles = "/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/all_regions"
# ds_list = "/home/tnphung/onedrive_documents/ctg/projects/gwas_celltype_atlas/analyses/modified_workflow/scz/all_regions/scrna_ds.txt"

datasets_df = read.table(ds_list)
colnames(datasets_df) = "datasets"
datasets = datasets_df$datasets

# unique_ds_count = c()

# step1 <- data.frame()
# for(ds in datasets){
#     f <- paste0(filedir, "magma_celltype_", ds, ".gsa.out")
#     if (file.exists(f)) {
# 	tmp <- fread(cmd=paste0("grep -v '^#' ", filedir, "magma_celltype_", ds, ".gsa.out"), data.table=F)
# 	if("FULL_NAME" %in% colnames(tmp)){
# 		tmp$VARIABLE <- tmp$FULL_NAME #convert VARIABLE to FULL_NAME
# 		tmp <- tmp[,-ncol(tmp)] #remove the FULL_NAME column
# 	}
#     unique_ds_count = c(unique_ds_count, tmp$VARIABLE)
#     tmp$VARIABLE <- paste0(strsplit(ds, "_")[[1]][1], "_", tmp$VARIABLE) #update variable to append the data set id at the front
# 	tmp <- tmp[order(tmp$P),] #order by p values
# 	tmp$ds <- ds #add ds column
# 	tmp$P.adj.pds <- p.adjust(tmp$P, method=adjPmeth) #add P.adj.pds column 
# 	if(nrow(step1)==0){step1 <- tmp}
# 	else{step1 <- rbind(step1, tmp)}
# }
# }

# step1$P.adj <- p.adjust(step1$P, method=adjPmeth) #add P.adj column 
# tmp_out <- step1[,c("ds", "VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P", "P.adj.pds", "P.adj")]
# colnames(tmp_out)[1:2] <- c("Dataset", "Cell_type")
# write.table(tmp_out, paste0(filedir, "magma_celltype_step1.txt"), quote=F, row.names=F, sep="\t")
# # rm(tmp_out)
# print(length(unique(unique_ds_count)))
# unique_ds_count_len = length(unique(unique_ds_count)) #TODO: fix this because the VARIABLE is now unique because I appended the dataset id
# tmp_out = tmp_out[which(tmp_out$P<(0.05/unique_ds_count_len)),] #modified to do bonferroni correction for the number of unique cell types
# print(nrow(tmp_out))
# if (nrow(tmp_out) == 0) {print("There is no significant cell type after step 1")}

# write.table(tmp_out, paste0(filedir, "magma_celltype_step1_sig.txt"), quote=F, row.names=F, sep="\t")

step1 = fread(paste0(filedir, "magma_celltype_step1.txt"))

step1$cond_state <- "single"
step1$cond_cell_type <- NA

step1 <- step1 %>%
  mutate(parent_ds = sapply(strsplit(Dataset, "_"), `[`, 2)) #add a collumn parent_ds for example Bhaduri 2021

step1$parent_ds <- as.character(step1$parent_ds)

if(nrow(step1)>1){
    step2_parent <- table(step1$parent_ds)
    filtered_names <- names(step2_parent[step2_parent > 1]) 
    if(length(which(step2_parent>1))==0){ #if there is only 1 count per parent study, then all are single
        step1$cond_state <- "single"
        step1$cond_cell_type <- NA
        step2_out <- data.frame()
    }else{
        step2_command <- c()
        for (parent in filtered_names) {
            tmp_df = step1 %>% filter(parent_ds== !!parent)
            all_cts = c()
            for (ct in tmp_df$Cell_type) {
                all_cts = c(all_cts, ct)
            }
            step2_command <- c(step2_command, paste0(magmadir, "/magma --gene-results ", filedir, "magma.genes.raw",
                        " --gene-covar ", magmafiles, "/celltype/", parent, ".txt --model condition-hide=Average direction=greater",
                        " analyse=list,", paste(all_cts, collapse=","), " joint-pairs",
                        " --out ", filedir, "magma_celltype_step2_", parent))

        }
        write.table(step2_command, paste0(filedir, "step2.sh"), quote=F, row.names=F, col.names=F)
        system(paste0("bash ", filedir, "step2.sh"))
        rm(step2_command)
        step1$cond_state <- NA
        step1$cond_state[step1$parent_ds %in% names(step2_parent)[step2_parent==1]] <- "single"
        step1$cond_cell_type <- NA
        step2_out <- data.frame()
        for(ds in names(step2_parent)[step2_parent>1]){
            if (ds != "Phan2024") {
            tmp.sig <- step1[step1$parent_ds==ds,] #tmp.sig is a df that contain significant rows per ds 
            #!!! if the file doesn't exist
            tmp <- fread(cmd=paste0("grep -v '^#' ", filedir, "magma_celltype_step2_", ds, ".gsa.out"), data.table=F)
            if("FULL_NAME" %in% colnames(tmp)){
            tmp$VARIABLE <- tmp$FULL_NAME
            tmp <- tmp[,-ncol(tmp)]
            }
            tmp$Marginal.P <- tmp.sig$P[match(tmp$VARIABLE, tmp.sig$Cell_type)] #get the marginal P from step 1
            tmp <- tmp[with(tmp, order(MODEL, Marginal.P)),] #for each pair, order by lowest marginal P first
            tmp$PS <- -log10(tmp$P)/-log10(tmp$Marginal.P) #calc PS 
            if(nrow(step2_out)==0){step2_out <- data.frame(tmp, ds=ds)}
            else{step2_out <- rbind(step2_out, data.frame(tmp, ds=ds))}
            checked <- c()
            while(length(checked)<nrow(tmp.sig)){
            top <- tmp.sig$Cell_type[!tmp.sig$Cell_type %in% checked][1]
            ### check if there is main drover for this cell type
            top.check <- c()
            for(m in unique(tmp$MODEL[tmp$VARIABLE==top])){
                t <- tmp[tmp$MODEL==m,]
                if(all(is.na(t$P))){next}
                if((t$PS[2]>=0.2 & t$PS[1]<0.2)){
                top.check <- c(top.check, t$VARIABLE[2])
                }
            }
            if(length(top.check)>0){
                top <- tmp.sig$Cell_type[tmp.sig$Cell_type %in% top.check][1]
            }
            checked <- c(checked, top)
            for(m in unique(tmp$MODEL[tmp$VARIABLE==top])){
                t <- tmp[tmp$MODEL==m,]
                if(t$VARIABLE[1]!=top){t <- t[2:1,]}
                if(all(is.na(t$P))){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "colinear"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "colinear", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "colinear-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(all(t$PS>=0.8)){
                #indep
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(all(t$P>=0.05)){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "joint-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(all(t$PS<0.2)){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "joint-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(t$PS[1]>=0.5 & t$P[2]>=0.05){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "main"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "main", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(t$PS[1]>=0.8 & t$PS[2]<0.2){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "main"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "main", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "partial-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(all(t$PS>=0.5)){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "partial-joint"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                #checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(all(t$PS>=0.2)){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "partial-joint-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else if(t$PS[1]>=0.2 & t$P[2]>=0.05){
                    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                    }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                    }
                    tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "joint-drop"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                    checked <- c(checked, t$VARIABLE[2])
                    tmp <- tmp[tmp$MODEL!=m,]
                }else if(t$PS[1]>=0.2 & t$PS[2]<0.2){
                if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$VARIABLE[2]
                }else{
                    tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
                    tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$VARIABLE[2], sep=";")
                }
                tmp.sig$cond_state[tmp.sig$Cell_type==t$VARIABLE[2]] <- "partial-joint-drop"
                tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$VARIABLE[2]] <- top
                checked <- c(checked, t$VARIABLE[2])
                tmp <- tmp[tmp$MODEL!=m,]
                }else{
                print(t)
                }
            }
            tmp <- tmp[!tmp$MODEL %in% unique(tmp$MODEL[tmp$VARIABLE %in% checked]),]
            }
            tmp.sig$cond_state[is.na(tmp.sig$cond_state)] <- "indep"
            step1$cond_state[step1$parent_ds==ds] <- tmp.sig$cond_state
            step1$cond_cell_type[step1$parent_ds==ds] <- tmp.sig$cond_cell_type
        }
        }
    }
    rm(tmp, tmp.sig, t)
    step1_out <- step1[,-12]
    # step1_out <- step1
    # step1_out <- step1_out[,c(7,1:6,8:11)]
    # colnames(step1_out)[1:2] <- c("Dataset", "Cell_type")
    step1_out$step3 <- with(step1_out, ifelse(grepl("drop", cond_state), 0, 1))
    write.table(step1_out, paste0(filedir, "step1_2_summary.txt"), quote=F, row.names=F, sep="\t")
    if(nrow(step2_out)>0){
        step2_out <- step2_out[,-2]
        step2_out <- step2_out[,c(10,1:9)]
        colnames(step2_out)[1:2] <- c("Dataset", "Cell_type")
        step2_out$MODEL <- rep(1:(nrow(step2_out)/2), each=2)
        write.table(step2_out, paste0(filedir, "magma_celltype_step2.txt"), quote=F, row.names=F, sep="\t")
    }
    rm(step1_out)
}

## Step 3
if(length(unique(step1$parent_ds))>1){
step1 <- step1[!grepl("drop", step1$cond_state),]
step3_ds <- unique(step1$parent_ds)
### condition average of the other dataset and pair-wise
step3_avg <- data.frame()
step3_cond <- data.frame()
for(i in 1:(length(step3_ds)-1)){
    ds1 <- step3_ds[i]
    exp1 <- fread(paste0(magmafiles, "/celltype/", ds1, ".txt"), data.table=F)
    exp1 <- exp1[,c("GENE", step1$Cell_type[step1$parent_ds==ds1], "Average")]
    colnames(exp1)[2:ncol(exp1)] <- paste(ds1, colnames(exp1)[2:ncol(exp1)], sep=":")
    colnames(exp1)[ncol(exp1)] <- "Average1"
    for(j in (i+1):length(step3_ds)){
    ds2 <- step3_ds[j]
    exp2 <- fread(paste0(magmafiles, "/celltype/", ds2, ".txt"), data.table=F)
    exp2 <- exp2[,c("GENE", step1$Cell_type[step1$parent_ds==ds2], "Average")]
    colnames(exp2)[2:ncol(exp2)] <- paste(ds2, colnames(exp2)[2:ncol(exp2)], sep=":")
    colnames(exp2)[ncol(exp2)] <- "Average2"
    exp <- cbind(exp1, exp2[match(exp1$GENE, exp2$GENE), -1])
    exp <- exp[!is.na(exp$Average2),]
    write.table(exp, paste0(filedir, "step3_exp.txt"), quote=F, row.names=F, sep="\t")
    step3_command <- paste0(magmadir, "/magma --gene-results ", filedir, "magma.genes.raw",
                            " --gene-covar ", filedir, "step3_exp.txt max-miss=0.1 --model condition-hide=Average1,Average2 direction=greater",
                            " --out ", filedir, "magma_celltype_step3_avg")
    res <- system(step3_command, ignore.stdout = T)
    if(res>0){
        #!!! implement for error
        print(paste("error: Average ", ds1, ds2))
    }else{
        tmp <- fread(cmd=paste0("grep -v '^#' ", filedir, "magma_celltype_step3_avg.gsa.out"), data.table=F)
        if("FULL_NAME" %in% colnames(tmp)){
        tmp$VARIABLE <- tmp$FULL_NAME
        tmp <- tmp[,-ncol(tmp)]
        }
        tmp$ds <- sub("(.+):.+", "\\1", tmp$VARIABLE)
        tmp$cond_ds <- ifelse(tmp$ds==ds1, ds2, ds1)
        tmp$VARIABLE <- sub(".+:(.+)", "\\1", tmp$VARIABLE)
        step3_avg <- rbind(step3_avg, tmp)
    }
    step3_command <- paste0(magmadir, "/magma --gene-results ", filedir, "magma.genes.raw",
                            " --gene-covar ", filedir, "step3_exp.txt max-miss=0.1 --model condition-hide=Average1,Average2 direction=greater joint-pairs",
                            " --out ", filedir, "magma_celltype_step3")
    res <- system(step3_command, ignore.stdout = T)
    if(res>0){
        tmp_ts <- colnames(exp)[-1]
        tmp_ts <- tmp_ts[!grepl("Average",tmp_ts)]
        tmp <- data.frame()
        for(tmp_i in 1:(length(tmp_ts)-1)){
        for(tmp_j in (tmp_i+1):length(tmp_ts)){
            tmp <- rbind(tmp, data.frame(VARIABLE=tmp_ts[c(tmp_i,tmp_j)], TYPE="COVAR", MODEL=1, NGENES=NA, BETA=NA, BETA_STD=NA, SE=NA, P=NA))
        }
        }
        tmp$MODEL <- rep(1:(nrow(tmp)/2), each=2)
    }else{
        tmp <- fread(cmd=paste0("grep -v '^#' ", filedir, "magma_celltype_step3.gsa.out"), data.table=F)
        if("FULL_NAME" %in% colnames(tmp)){
        tmp$VARIABLE <- tmp$FULL_NAME
        tmp <- tmp[,-ncol(tmp)]
        }
    }
    tmp$ds <- sub("(.+):.+", "\\1", tmp$VARIABLE)
    tmp$VARIABLE <- sub(".+:(.+)", "\\1", tmp$VARIABLE)
    check.model <- with(tmp, aggregate(ds, list(MODEL), function(x){length(unique(x))}))
    tmp <- tmp[tmp$MODEL %in% check.model$Group.1[check.model$x==2],]
    step3_cond <- rbind(step3_cond, tmp)
    }
}

rm(exp, exp1, exp2)
system(paste0("rm ", filedir, "step3_exp.txt"))

### add within dataset conditional analyses
step3_cond <- step3_cond[,-2]
step3_cond <- step3_cond[,c(8,1:7)]
colnames(step3_cond)[1:2] <- c("Dataset", "Cell_type")
for(ds in step3_ds){
    if(length(which(step1$parent_ds==ds))>1){
    tmp <- step2_out[step2_out$Dataset==ds & step2_out$Cell_type %in% step1$Cell_type[step1$parent_ds==ds],]
    tmp.model <- table(tmp$MODEL)
    tmp <- tmp[tmp$MODEL %in% names(tmp.model)[tmp.model==2],]
    if(nrow(tmp)>0){
        step3_cond <- rbind(step3_cond, tmp[,1:8])
    }
    }
}
step3_cond$MODEL <- rep(1:(nrow(step3_cond)/2), each=2)

step3_cond$label <- paste(step3_cond$Dataset, step3_cond$Cell_type, sep=":")
step3_cond$label[seq(1,nrow(step3_cond),2)] <- paste(step3_cond$Dataset[seq(2,nrow(step3_cond),2)], step3_cond$label[seq(1,nrow(step3_cond),2)], sep=":")
step3_cond$label[seq(2,nrow(step3_cond),2)] <- paste(step3_cond$Dataset[seq(1,nrow(step3_cond),2)], step3_cond$label[seq(2,nrow(step3_cond),2)], sep=":")
step3_avg$label <- paste(step3_avg$cond_ds, step3_avg$ds, step3_avg$VARIABLE, sep=":")
tmp <- step3_avg[match(step3_cond$label, step3_avg$label),c(4:7,9)]
colnames(tmp)[ncol(tmp)] <- "ds"
colnames(tmp) <- paste0("CDM.", colnames(tmp))
step3_cond <- cbind(step3_cond[,-ncol(step3_cond)], tmp)
step3_cond$Marginal.P <- step1$P[match(paste(step3_cond$Dataset, step3_cond$Cell_type, sep=":"), paste(step1$parent_ds, step1$Cell_type, sep=":"))]
step3_cond$PS <- -log10(step3_cond$P)/-log10(with(step3_cond, ifelse(is.na(CDM.P), Marginal.P, CDM.P)))
step3_cond$PS.avg <- -log10(step3_cond$CDM.P)/-log10(step3_cond$Marginal.P)
write.table(step3_cond, paste0(filedir, "magma_celltype_step3.txt"), quote=F, row.names=F, sep="\t")
}
