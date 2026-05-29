library(ggplot2)

filter_degs = function(degs, p_val_limit = 0.05, fold_limit = 0, top = 50, pct.1_limit_low = 0, pct.2_limit_low = 0, pct.1_limit_high = 1.0, pct.2_limit_high = 1.0) {
  degs_1 = subset(degs, p_val < p_val_limit)
  degs_2 = subset(degs_1, avg_log2FC > fold_limit)
  degs_3 = subset(degs_2, pct.1 > pct.1_limit_low)
  degs_4 = subset(degs_3, pct.2 > pct.2_limit_low)
  degs_5 = subset(degs_4, pct.1 < pct.1_limit_high)
  degs_6 = subset(degs_5, pct.2 < pct.2_limit_high)
  degs_filtered = degs_6
  return (degs_filtered)
}

assign_GO_cluster <- function(tmp, markers, file) {
  all_genes <- rownames(tmp)
  
  per.cluster <- split(markers, markers$cluster)
  sig.up <- lapply(per.cluster, function(x) as.vector(x[x$avg_log2FC > 0 && x$p_val_adj<0.01,]$gene))
  
  combined_res <- vector(mode = "list", length = length(sig.up))
  
  for (i in 2:length(sig.up)) {
    geneList <- factor(as.integer (all_genes %in% sig.up[[i]]))
    names(geneList) <- all_genes
    
    GO       <- new("topGOdata",
                    description = "Cluster", ontology = "BP",
                    allGenes = geneList,
                    nodeSize = 10,
                    annot=annFUN.org, mapping="org.Dr.eg.db", ID = "alias")
    
    fisher <- runTest(GO, algorithm = "classic", statistic = "fisher")
    res <- GenTable(GO, classicFisher = fisher, topNodes = 10)
    combined_res[[i]] <- res
    #    kable(res, row.names = F)
    showSigOfNodes(GO, score(fisher), firstSigNodes = 5, useInfo = "def")
  } 
  for(i in 1:length(combined_res)) {
    write.table(combined_res[[i]], file = paste0(file,i))
  }
}

assign_GO_cluster_pair <- function(tmp, markers, file) {
  all_genes <- rownames(tmp)
  sig.up <- as.vector(markers[markers$p_val_adj<0.01,]$gene)
  geneList <- factor(as.integer (all_genes %in% sig.up))
  names(geneList) <- all_genes
  
  GO       <- new("topGOdata",
                  description = "Cluster", ontology = "BP",
                  allGenes = geneList,
                  nodeSize = 10,
                  annot=annFUN.org, mapping="org.Dr.eg.db", ID = "alias")
  
  fisher <- runTest(GO, algorithm = "classic", statistic = "fisher")
  res <- GenTable(GO, classicFisher = fisher, topNodes = 10)
  kable(res, row.names = F)
  showSigOfNodes(GO, score(fisher), firstSigNodes = 5, useInfo = "def")
  
  write.table(res, file = paste0(file))
}

plotClusters_marker = function(tmp, marker){
  clusters = as.factor(levels(tmp))
  counts = matrix(nrow=1, ncol=0)
  for(i in 1:length(clusters)) {
    cluster = clusters[i]
    counts_cluster = matrix(nrow=1, ncol=ncol(tmp@assays$RNA@data[,tmp@active.ident == cluster]))
    counts_cluster[1,0:ncol(counts_cluster)] = tmp@assays$RNA@data[marker,tmp@active.ident == cluster]
    colnames(counts_cluster) = rep(cluster, ncol(counts_cluster))
    counts = cbind(counts, counts_cluster)
  }
  
  
  counts_df <- data.frame(cluster = colnames(counts),colSums_counts=colSums(counts))
  counts_df$idu <- as.numeric(row.names(counts_df))
  
  p<-ggplot(data=counts_df, aes(x=idu, y=colSums_counts, fill=cluster)) + 
    geom_bar(stat="identity")+theme_minimal() + xlab("Cluster") + ylab("Counts")
 
  return(p)
}

overlap_degs = function(degs_1, degs_2) {
  degs_overlap = degs_1[(degs_1$gene %in% degs_2$gene),]
  return (degs_overlap)
}

overlap_only_first = function(degs_1, degs_2) {
  degs_only = degs_1[!(degs_1$gene %in% degs_2$gene),]
  return (degs_only)
}