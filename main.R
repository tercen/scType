library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(reshape)
library(Seurat)
library(HGNChelper)
library(openxlsx)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load auto-detection function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

# load database
db_ = "./ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

ctx <- tercenCtx()
scRNAseqData<-as.matrix(ctx)

ctx %>% dplyr::select(.ci) 
clusters<-as.matrix(unique(ctx$select(unlist(list(ctx$colors, '.ci')))))
colnames(clusters)<-c(".cluster",".ci")

rownames(scRNAseqData)<-as.matrix(rselect(ctx))
colnames(scRNAseqData)<-as.matrix(ctx%>% select(.ci)%>%unique(.))

#scRNAseqData <- read.csv("./pbmc3k.csv",row.names = 1)
# prepare gene sets
gs_list <- suppressWarnings(suppressMessages(gene_sets_prepare(db_, tissue)))
# assign cell types
es.max <- sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(clusters[,".cluster"]), function(cl){
  cl.tmp<-clusters[clusters[,1]==cl,]
  es.max.cl = sort(rowSums(es.max[,colnames(es.max)==as.integer(cl.tmp[,2]), drop=FALSE]), decreasing = !0)
  head(data.frame(.cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(clusters==cl)), 10)
}))

merge.cl_results<-merge(cL_resutls,clusters,all=TRUE) 
sctype_scores = merge.cl_results %>% group_by(.cluster) %>% top_n(n = 1, wt = scores) 

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

sctype_scores<-sctype_scores%>%
  mutate(.ci = as.integer(.ci))%>%
  ctx$addNamespace()

names(sctype_scores)[1]<-ctx$colors[[1]]

sctype_scores%>%
  ctx$save()

