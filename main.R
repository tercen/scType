library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(reshape)
library(Seurat)
library(HGNChelper)
library(openxlsx)
source("gene_sets_prepare.R")
source("sctype_score_.R")

### FUNCTION
gene_sets_prepare_custom <- function(table_in, cell_type){
  
  cell_markers = table_in
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

#######

ctx <- tercenCtx()

### load database
doc.id.tmp<-as_tibble(ctx$select())
doc.id<-doc.id.tmp[[grep("documentId" , colnames(doc.id.tmp))]][1]

# prepare gene sets
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

if(is.null(doc.id)){
  db_ = "./ScTypeDB_full.xlsx"
  gs_list <- suppressWarnings(suppressMessages(gene_sets_prepare(db_, tissue)))
}else{
  doc.id<-doc.id.tmp[[grep("documentId" , colnames(doc.id.tmp))]][1]
  table.pop<-ctx$client$tableSchemaService$select(doc.id)
  tbl_pop<-as_tibble(table.pop)
  
  gs_list<-suppressWarnings(suppressMessages(gene_sets_prepare_custom(tbl_pop, tissue)))
}

scRNAseqData<-as.matrix(ctx)

ctx %>% dplyr::select(.ci) 
clusters<-as.data.frame(unique(ctx$select(unlist(list(ctx$colors, '.ci')))))
colnames(clusters)<-c(".cluster",".ci")
clusters[,".ci"]<-as.integer(clusters[,".ci"])

rownames(scRNAseqData)<-as.matrix(rselect(ctx))
colnames(scRNAseqData)<-as.matrix(ctx%>% select(.ci)%>%unique(.))

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

