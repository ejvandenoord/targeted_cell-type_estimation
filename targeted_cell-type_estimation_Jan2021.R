### If you use this code please cite the following publication:
### A targeted solution for estimating the cell-type composition of bulk samples 
### Edwin van den Oord, Lin Y. Xie, Charles J. Tran, Min Zhao, Karolina A. Aberg
### @@@ add full citation when available @@@
  

rm(list = ls())

### Set work directory
work_dir <- "your_work_directory"
setwd(work_dir)

### Read in required info
music             = T # if T music else standard OLS estimate
data              = read.csv('example_data.csv',row.names = 1) # Table S3 with exclusion of the 2nd row in the table
celltype_samples  = readLines("example_celltype_samples.txt") # samples from sorted cells to be used as reference panel
bulk_samples      = readLines("example_bulk_samples.txt") # samples for which to estimate celltype proportion
pool_exp          = read.csv("example_expected_proportions.csv",row.names = 1,head=T) # expected celltype proportions, Table S5
amplicons         = read.csv("example_amplicons.csv",stringsAsFactors = F) # all amplicons in the data file,, first row for each amplicon in Table S1 where  Chromosome is labeled "chr" and Start is labeled "pos". 
used_sites        = readLines("example_used_sites_best_model.txt") # specific sites in the amplicons to be used
label             = "best_model" # a string describing the used model

### Load required libraries and code
source("estimate_celltype_proportions_July2020_functions.R")
library(MuSiC)
library(Biobase) # need to create "eset" data set used by Music
library(xbioc)   # used by Music
library(mltools) # to calulate rmse

##### Prepare information for estimation of celltype proportions 
n_short_pos      = 4
amplicon_type    = parse_amplicon_names(amplicons,colnames(data),n_short_pos )
ind              = match( colnames(data),paste0('chr',amplicon_type$chr,'_',amplicon_type$pos) )
amplicon_type    = amplicon_type[ind,]
colnames(data)   = paste0(colnames(data),'_',amplicon_type$amplicon_type)
data<-data[,used_sites]
ref_data         = data[ match(celltype_samples,rownames(data)),order(colnames(data))]
bulk_data        = data[ match(bulk_samples,rownames(data)),order(colnames(data))]
bulk_data.eset   = ExpressionSet(assayData = as.matrix(t(bulk_data)))

### Create ref panel
temp = create_reference( ref_data )
ref_panel      = temp[[1]]
ref_panel_eset = temp[[2]]
write.csv(ref_panel,'ref_panel.csv',row.names=T) # Reference panel 

### Estimate cell type proportions
if (music) {
  cell_types = colnames(ref_panel)
  temp       = music_prop(bulk.eset=bulk_data.eset, sc.eset=ref_panel_eset,clusters='clusters',samples='cell_id',select.ct=cell_types,verbose=T)
  estimates  = temp$Est.prop.weighted } else estimates = reg_celltype(bulk_data,ref_panel,F)$estimates 

#### Calculate rmse
# Across the entire set of proportions
tmp_all<-rmse(c(pool_exp[,1],pool_exp[,2],pool_exp[,3],pool_exp[,4]),
              as.numeric(c(estimates[1:8,1],estimates[1:8,2],estimates[1:8,3],estimates[1:8,4])))
# For each cell type separately
tmp_cd3<-rmse(pool_exp$CD3,estimates[1:8,"CD3"])
tmp_cd15<-rmse(pool_exp$CD15,estimates[1:8,"CD15"])
tmp_cd14<-rmse(pool_exp$CD14,estimates[1:8,"CD14"])
tmp_cd19<-rmse(pool_exp$CD19,estimates[1:8,"CD19"])
# Mean across the cell types
tmp_all_mean_rmse<-mean(c(tmp_cd3,tmp_cd15,tmp_cd14,tmp_cd19))
# Summary
RMSE<-c(tmp_all,tmp_cd3,tmp_cd14,tmp_cd15,tmp_cd19,tmp_all_mean_rmse)

#### Calculate correlation
# Across the entire set of proportions
tmp_all<-cor(c(pool_exp[,1],pool_exp[,2],pool_exp[,3],pool_exp[,4]),
             as.numeric(c(estimates[1:8,1],estimates[1:8,2],estimates[1:8,3],estimates[1:8,4])),
             use="complete.obs",method="pearson")
# For each cell type separately - (the diagonal o fthe matrix)
tmp<-cor(pool_exp,estimates[1:8,],use="complete.obs",method="pearson")
# Mean across the cell types
tmp_all_mean_cor<-mean(c(tmp[1],tmp[11],tmp[6],tmp[16])) #cd3, cd14, cd15, cd19
# Summary
cor<-c(tmp_all,tmp[c(1,6,11,16)],tmp_all_mean_cor)

### Summary of the evaluation
eval_res<-c("label","rmse_all_data","rmse_cd3", "rmse_cd14", "rmse_cd15", "rmse_cd19", "rmse_mean_across_cds",
                   "cor_all_data","cor_cd3", "cor_cd14", "cor_cd15", "cor_cd19", "cor_mean_across_cds")
tmp<-c(label,RMSE,cor)
eval_res<-rbind(eval_res,tmp)
write.table(eval_res,'eval_res.csv',row.names=F,col.names=F,quote=F,sep=",") # summary file with evaluation measures



