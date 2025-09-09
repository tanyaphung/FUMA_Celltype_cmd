library(dplyr)
library(data.table)
library(argparse)

# run this script after the file `magma_celltype_step3.txt` is generated

# Create argument parser
parser <- ArgumentParser(description = '')

# Add arguments
parser$add_argument("--step1_2_sum_fp", type = "character", required = TRUE,
                    help = "Path to the file step1_2_summary.txt")
parser$add_argument("--step3_fp", type = "character", required = TRUE,
                    help = "Path to the file magma_celltype_step3.txt")
parser$add_argument("--out_fp", type = "character", required = TRUE,
                    help = "Path to the output file magma_celltype_step3_filtered.txt")

# Parse arguments
args <- parser$parse_args()

# step1_2_sum_fp = "scz/all_regions/step1_2_summary.txt"
# step3_fp = "scz/all_regions/magma_celltype_step3.txt"
# out_fp = "scz/all_regions/magma_celltype_step3_filtered.txt"

step1_2_sum_fp = args$step1_2_sum_fp
step3_fp = args$step3_fp
out_fp = args$out_fp

infile = fread(step1_2_sum_fp)
infile_filtered = infile %>% filter(step3==1)
tmp.sig = infile_filtered %>% select("Dataset", "Cell_type", "P", "P.adj.pds", "P.adj")

tmp.sig$cond_state <- "single"
tmp.sig$cond_cell_type <- NA

step3_out <- data.frame()
tmp <- fread(step3_fp)
tmp <- tmp[with(tmp, order(MODEL, Marginal.P)),] #for each pair, order by lowest marginal P first
step3_out <- data.frame(tmp)

checked <- c()
while(length(checked)<nrow(tmp.sig)){
top <- tmp.sig$Cell_type[!tmp.sig$Cell_type %in% checked][1]
### check if there is main driver for this cell type
top.check <- c()
for(m in unique(tmp$MODEL[tmp$Cell_type==top])){
    t <- tmp[tmp$MODEL==m,]
    if(all(is.na(t$P))){next}
    if((t$PS[2]>=0.2 & t$PS[1]<0.2)){
    top.check <- c(top.check, t$Cell_type[2])
    }
}
if(length(top.check)>0){
    top <- tmp.sig$Cell_type[tmp.sig$Cell_type %in% top.check][1]
}
checked <- c(checked, top)
for(m in unique(tmp$MODEL[tmp$Cell_type==top])){
    t <- tmp[tmp$MODEL==m,]
    if(t$Cell_type[1]!=top){t <- t[2:1,]}
    if(all(is.na(t$P))){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "colinear"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "colinear", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "colinear-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(all(t$PS>=0.8)){
    #indep
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(all(t$P>=0.05)){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "joint-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(all(t$PS<0.2)){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "joint-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(t$PS[1]>=0.5 & t$P[2]>=0.05){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "main"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "main", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(t$PS[1]>=0.8 & t$PS[2]<0.2){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "main"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "main", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "partial-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(all(t$PS>=0.5)){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "partial-joint"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    #checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(all(t$PS>=0.2)){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "partial-joint-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else if(t$PS[1]>=0.2 & t$P[2]>=0.05){
        if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
        }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
        }
        tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "joint-drop"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
        checked <- c(checked, t$Cell_type[2])
        tmp <- tmp[tmp$MODEL!=m,]
    }else if(t$PS[1]>=0.2 & t$PS[2]<0.2){
    if(is.na(tmp.sig$cond_state[tmp.sig$Cell_type==top])){
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- "partial-joint"
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- t$Cell_type[2]
    }else{
        tmp.sig$cond_state[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_state[tmp.sig$Cell_type==top], "partial-joint", sep=";")
        tmp.sig$cond_cell_type[tmp.sig$Cell_type==top] <- paste(tmp.sig$cond_cell_type[tmp.sig$Cell_type==top], t$Cell_type[2], sep=";")
    }
    tmp.sig$cond_state[tmp.sig$Cell_type==t$Cell_type[2]] <- "partial-joint-drop"
    tmp.sig$cond_cell_type[tmp.sig$Cell_type==t$Cell_type[2]] <- top
    checked <- c(checked, t$Cell_type[2])
    tmp <- tmp[tmp$MODEL!=m,]
    }else{
    print(t)
    }
}
tmp <- tmp[!tmp$MODEL %in% unique(tmp$MODEL[tmp$Cell_type %in% checked]),]
}
tmp.sig$cond_state[is.na(tmp.sig$cond_state)] <- "indep"

write.table(tmp.sig, out_fp, quote = FALSE, row.names = FALSE, sep="\t")

print("Finish")