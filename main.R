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
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

ctx <- tercenCtx()
scRNAseqData<-as.matrix(ctx)

ctx %>% dplyr::select(.ci) 
clusters<-as.data.frame(unique(ctx$select(unlist(list(ctx$colors, '.ci')))))
colnames(clusters)<-c(".cluster",".ci")
clusters[,".ci"]<-as.integer(clusters[,".ci"])

rownames(scRNAseqData)<-as.matrix(rselect(ctx))
colnames(scRNAseqData)<-as.matrix(ctx%>% select(.ci)%>%unique(.))

# prepare gene sets
gs_list <- suppressWarnings(suppressMessages(gene_sets_prepare(db_, tissue)))
# assign cell types
es.max <- sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(clusters[,".cluster"]), function(cl){
  cl.tmp<-clusters[clusters[,1]==cl,]
  es.max.cl = sort(rowSums(es.max[,colnames(es.max)==as.integer(cl.tmp[,2]), drop=FALSE]), decreasing = !0)
  head(data.frame(.cluster = cl, population = names(es.max.cl), scores = es.max.cl, ncells = sum(clusters==cl)), 10)
}))

merge.cl_results<-merge(cL_resutls,clusters,all=TRUE)
colnames(es.max)<-seq(from=0,to=length(colnames(es.max))-1)
es.max.long <- melt(es.max)
colnames(es.max.long)<-c("population",".ci","solo_score")
merge.tmp<-merge(es.max.long,clusters,all=TRUE)
merge.solo.cl_results<-merge(merge.tmp,merge.cl_results, by = c(".ci",".cluster","population"),all=TRUE)

sctype_scores = merge.solo.cl_results %>% group_by(.cluster) %>% top_n(n = 1, wt = scores) 

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$population <- as.character(sctype_scores$population)
sctype_scores$population[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4]<- "Unknown"
sctype_scores$population <- as.factor(sctype_scores$population)

sctype_scores<-sctype_scores%>%
  mutate(.ci = as.integer(.ci))%>%
  ctx$addNamespace()

sctype_scores%>%
  ctx$save()

