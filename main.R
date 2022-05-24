library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(reshape)
library(Seurat)
library(HGNChelper)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load auto-detection function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

# load database
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

ctx <- tercenCtx()
scRNAseqData<-as.matrix(ctx)

#.ri <- seq(from = 0, to = length(as.matrix(rselect(ctx))) - 1)
#.ci <- seq(from = 0, to = length(as.matrix(cselect(ctx))) - 1)
#rownames(scRNAseqData)<-.ri 
#colnames(scRNAseqData)<-.ci
rownames(scRNAseqData)<-as.matrix(rselect(ctx))
colnames(scRNAseqData)<-as.matrix(cselect(ctx))

#scRNAseqData <- read.csv("./pbmc3k.csv",row.names = 1)
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
# assign cell types
es.max <- sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

es.max.long <- melt(es.max)
df_test<-data.frame(.ci = seq(from=0,to=length(es.max.long)-1), f=es.max.long) 
df_test %>%
  ctx$addNamespace() %>%
  ctx$save()