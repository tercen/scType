library(tercen)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(reshape)
library(Seurat)
library(HGNChelper)
library(openxlsx)
#source("gene_sets_prepare.R")
#source("sctype_score_.R")


#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load auto-detection function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

options("tercen.workflowId" = "bb1b63cf1a12a9aabb0ce7b33d018aee")
#options("tercen.stepId"     = "bdf24f26-8dd1-4cca-ab9c-80a09c2a6a64")
options("tercen.stepId"     = "9224c792-7098-4b97-9ff7-0edadfaa1364")

options("tercen.workflowId" = "9af47f50c3160717a5378977da00e254")
options("tercen.stepId"     = "cc131013-caed-4fb2-a3c0-8392b16ee6ad")




getOption("tercen.workflowId")
getOption("tercen.stepId")

### FUNCTION

table_in<-tbl_pop
cell_type<-tissue
gene_sets_prepare_custom <- function(table_in, cell_type){
  
  cell_markers = table_in
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    # if(length(markers_all) > 0){
    #   markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
    #   paste0(markers_all, collapse=",")
    # } else {
    #   ""
    # }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    # if(length(markers_all) > 0){
    #   markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
    #   paste0(markers_all, collapse=",")
    # } else {
    #   ""
    # }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

sctype_score_custom <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
    rownames(scRNAseqData) = gsub(" ","",rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
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
  #db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  gs_list <- suppressWarnings(suppressMessages(gene_sets_prepare(db_, tissue)))
}else{
  #doc = ctx$client$fileService$get(doc.id)
  doc.id<-doc.id.tmp[[grep("documentId" , colnames(doc.id.tmp))]][1]
  table.pop<-ctx$client$tableSchemaService$select(doc.id)
  tbl_pop<-as_tibble(table.pop)
  #tbl_pop<-as.list(tbl_pop)

  gs_list<-suppressWarnings(suppressMessages(gene_sets_prepare_custom(tbl_pop, tissue)))
}


#split(x = testse, f = testse$columns[[1]])

scRNAseqData<-as.matrix(ctx)

ctx %>% dplyr::select(.ci) 
clusters<-as.data.frame(unique(ctx$select(unlist(list(ctx$colors, '.ci')))))
colnames(clusters)<-c(".cluster",".ci")
clusters[,".ci"]<-as.integer(clusters[,".ci"])

#rownames(scRNAseqData)<-.ri 
#colnames(scRNAseqData)<-.ci
rownames(scRNAseqData)<-as.matrix(rselect(ctx))
#colnames(scRNAseqData)<-as.matrix(cselect(ctx))
colnames(scRNAseqData)<-as.matrix(ctx%>% select(.ci)%>%unique(.))

#scRNAseqData <- read.csv("./pbmc3k.csv",row.names = 1)

# assign cell types
es.max <- sctype_score_custom(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(clusters[,".cluster"]), function(cl){
  cl.tmp<-clusters[clusters[,1]==cl,]
  es.max.cl = sort(rowSums(es.max[,colnames(es.max)==as.integer(cl.tmp[,2]), drop=FALSE]), decreasing = !0)
  head(data.frame(.cluster = cl, population = names(es.max.cl), scores = es.max.cl, ncells = sum(clusters==cl)), 10)
}))

merge.cl_results<-merge(cL_resutls,clusters,all=TRUE)
colnames(es.max)<-seq(from=0,to=length(colnames(es.max))-1)
es.max.long <- melt(es.max)
#df_test<-data.frame(.ci = seq(from=0,to=length(rownames(es.max.long))-1), f=es.max.long) 
colnames(es.max.long)<-c("population",".ci","solo_score")
merge.tmp<-merge(es.max.long,clusters,all=TRUE)
merge.solo.cl_results<-merge(merge.tmp,merge.cl_results, by = c(".ci",".cluster","population"),all=TRUE)

#merge.cl_results[(merge.cl_results[,".ci"]==0),]
#merge.solo.cl_results[(merge.solo.cl_results[,".ci"]==0),]
#sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = merge.solo.cl_results %>% group_by(.cluster) %>% top_n(n = 1, wt = scores) 

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$population <- as.character(sctype_scores$population)
sctype_scores$population[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4]<- "Unknown"
sctype_scores$population <- as.factor(sctype_scores$population)
#print(sctype_scores[,1:3])


#df_test<-data.frame(.ci = seq(from=0,to=length(rownames(es.max.long))-1), f=es.max.long) 

#df_test %>%

#sctype_scores<-sctype_scores %>% 
#  select(everything()) %>% 
#  mutate(across(.ci, as.integer))
#names(sctype_scores)[1]<-ctx$colors[[1]]

sctype_scores<-sctype_scores%>%
  mutate(.ci = as.integer(.ci))%>%
  ctx$addNamespace()

#names(sctype_scores)[1]<-ctx$colors[[1]]


sctype_scores%>%
  ctx$save()

