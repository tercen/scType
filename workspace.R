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

options("tercen.workflowId" = "bb1b63cf1a12a9aabb0ce7b33d018aee")
options("tercen.stepId"     = "bdf24f26-8dd1-4cca-ab9c-80a09c2a6a64")
#options("tercen.stepId"     = "eb27b255-8edc-4597-ad32-8154aa690062")

getOption("tercen.workflowId")
getOption("tercen.stepId")

# load database
db_ = "./ScTypeDB_full.xlsx"
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

ctx <- tercenCtx()
scRNAseqData<-as.matrix(ctx)

ctx %>% dplyr::select(.ci) 
clusters<-as.matrix(unique(ctx$select(unlist(list(ctx$colors, '.ci')))))
colnames(clusters)<-c(".cluster",".ci")
#rownames(scRNAseqData)<-.ri 
#colnames(scRNAseqData)<-.ci
rownames(scRNAseqData)<-as.matrix(rselect(ctx))
#colnames(scRNAseqData)<-as.matrix(cselect(ctx))
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
#sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = merge.cl_results %>% group_by(.cluster) %>% top_n(n = 1, wt = scores) 

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])

#es.max.long <- melt(es.max)
#colnames(es.max.long)<-c("population","ID","sctype_score")
#df_test<-data.frame(.ci = seq(from=0,to=length(rownames(es.max.long))-1), f=es.max.long) 

#df_test %>%

#sctype_scores<-sctype_scores %>% 
#  select(everything()) %>% 
#  mutate(across(.ci, as.integer))
names(sctype_scores)[1]<-ctx$colors[[1]]

sctype_scores<-sctype_scores%>%
  mutate(.ci = as.integer(.ci))%>%
  ctx$addNamespace()

names(sctype_scores)[1]<-ctx$colors[[1]]


sctype_scores%>%
  ctx$save()

