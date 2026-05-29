#1. MULTI or two way (below) DEG cluster comparison

DEG_Multi36 <- FindAllMarkers(adata_36h, assay = "RNA")

#2. Cut offs choose

#Besser, laut Anna:
DEG_Multi36_adjpval01 <- subset(DEG_Multi36, p_val_adj < 0.1)
#or use also p_val_adj < 0.01 or 1
DEG_Multi36_adjpval01_FC <- subset(DEG_Multi36_adjpval01, avg_log2FC > 0)

write.csv(DEG_Multi36_adjpval01_FC, "/Users/margl916/Documents/DEG_Multi36_adjpval01_FC.csv", row.names=TRUE)

#Muss per cluster sein:
DEG_Multi36_adjpval01_FC_cluster0 <- DEG_Multi36_adjpval01_FC %>% filter(cluster == "0")

tibble(DEG_Multi36_adjpval01_FC_cluster0)

DEG_Multi36_adjpval01_FC_cluster0 = DEG_Multi36_adjpval01_FC_cluster0 %>% dplyr::arrange(desc(avg_log2FC))

write.csv(DEG_Multi36_adjpval01_FC_cluster0, "/Users/margl916/Documents/DEG_Multi36_adjpval01_FC_cluster0.csv", row.names=TRUE)


#3. Which method to use (Agata or other/Anna-similar)

#USE OTHER METHOD!!!!! - gives more results !!! (But maybe more unspecific?)
#Seems that if this Other method used, similar results to 2way/individual cluster comparisons or multi way cluster comparison!!!!
#also much quicker
  
#Other method:
  
#convert zebrafish gene symbols to Entrez gene ID's:
library(org.Dr.eg.db)
keytypes(org.Dr.eg.db)
library(clusterProfiler)


DEG_Multi36_adjpval01_FC_cluster0_list = scan("/Users/margl916/Documents/DEG_Multi36_adjpval01_FC_cluster0_list.txt", character(), quote = "")
t <- c(DEG_Multi36_adjpval01_FC_cluster0_list)

et <- bitr(t, fromType="SYMBOL", toType=(c("ENTREZID","PATH","GO","ALIAS","GENENAME")), OrgDb="org.Dr.eg.db")
head(et)
write.csv(et, "/Users/margl916/Documents/et.csv", row.names=TRUE)

entrezID_et <- select(et, ENTREZID)

#remove duplicate rows across entire data frame 
entrezID_et_single <- entrezID_et%>%distinct(.keep_all = TRUE)
#save
write.csv(entrezID_et_single, "/Users/margl916/Documents/entrezID_et_single.csv", row.names=TRUE)
#scan as text file only one column - entrezID_et_single
entrezID_et_single_list = scan("/Users/margl916/Documents/entrezID_et_single_list.txt", character(), quote = "")


#4. Choose pvalueCutoff 

ans.go <- enrichGO(gene = entrezID_et_single_list, ont = "BP",
                   OrgDb ="org.Dr.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.5)

#check if needed
ans.go

#or other gseaPvalCO 0.01 or 0.1 or 1

tab.go <- as.data.frame(ans.go)
View(tab.go)
write.csv(tab.go, "/Users/margl916/Documents/tab.go.csv", row.names=TRUE)
#RENAME FILE ACCORDINGLY


#visualization
library(enrichplot)
dotplot(ans.go, showCategory=20) + ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1)
#file save as 600x600


dotplot(ans.go, showCategory=10, orderBy = "qvalue") + 
  ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1) + 
  scale_y_discrete(limits = rev)


dotplot(ans.go, showCategory= c("blood vessel development", "angiogenesis", "endothelial cell differentiation", "venous blood vessel development", "lymph vessel development", "lymphangiogenesis", "sprouting angiogenesis", "endothelial cell migration", "vascular endothelial growth factor receptor signaling pathway", "branching involved in blood vessel morphogenesis", "blood vessel lumenization", "endothelial cell proliferation", "DNA replication", "mitotic cell cycle process", "cell population proliferation", "otolith development"), orderBy = "qvalue") + 
  ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1) + 
  scale_y_discrete(limits = rev)
#save as 800x354 #...all

dotplot(ans.go, showCategory= c("blood vessel development", "angiogenesis", "endothelial cell differentiation", "venous blood vessel development", "lymph vessel development", "lymphangiogenesis", "sprouting angiogenesis", "endothelial cell migration", "vascular endothelial growth factor receptor signaling pathway", "branching involved in blood vessel morphogenesis", "blood vessel lumenization", "endothelial cell proliferation", "DNA replication", "mitotic cell cycle process", "cell population proliferation", "otolith development")) + 
  ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1) + 
  scale_y_discrete(limits = rev)
#save as 800x354 #...all

#"regulation of cell population proliferation"

dotplot(ans.go, showCategory= c("blood vessel development", "angiogenesis", "sprouting angiogenesis", "response to growth factor", "regulation of catalytic activity", "regulation of kinase activity", "branching involved in blood vessel morphogenesis", "blood vessel lumenization", "endothelial cell proliferation", "mitotic cell cycle process")) + 
 ggtitle("GO BP") +
 theme(text = element_text(size = 10)) + 
 theme(axis.text.y = element_text(size=10)) +  
 viridis::scale_fill_viridis(direction = -1)


dotplot(ans.go, showCategory= c("DNA replication", "chromosome organization", "mitotic cell cycle process", "nuclear division", "DNA recombination", "chromatin organization", "chromatin remodeling", "regulation of cell cycle process")) + 
  ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1)
#save as 500x350 #...all


cnetplot(ans.go, color.params = list(foldChange = DEG_Multi36_adjpval01_FC_cluster0_list, list(edge = TRUE)), showCategory = c("blood vessel development", "lymphangiogenesis"), circular = TRUE) +
 theme(text = element_text(size = 10)) + 
 theme(axis.text.y = element_text(size=10)) +  
 viridis::scale_fill_viridis(direction = -1)


upsetplot(ans.go, showCategory=20) + ggtitle("GO BP")


#1. Multi (above) or TWO WAY DEG cluster comparison

Comp36LEC0vs2 <- subset(adata_36h, idents = c("0", "2"))

DEG_Comp36LEC0vs2 <- FindAllMarkers(Comp36LEC0vs2, assay = "RNA")

#2. Cut offs choose

#Besser, laut Anna:
DEG_Comp36LEC0vs2_adjpval01 <- subset(DEG_Comp36LEC0vs2, p_val_adj < 0.1)
#or use also p_val_adj < 0.01 or 1
DEG_Comp36LEC0vs2_adjpval01_FC <- subset(DEG_Comp36LEC0vs2_adjpval01, avg_log2FC > 0)

write.csv(DEG_Comp36LEC0vs2_adjpval01_FC, "/Users/margl916/Documents/DEG_Comp36LEC0vs2_adjpval01_FC.csv", row.names=TRUE)

#Muss per cluster sein:
DEG_Comp36LEC0vs2_adjpval01_FC_cluster2 <- DEG_Comp36LEC0vs2_adjpval01_FC %>% filter(cluster == "2")

tibble(DEG_Comp36LEC0vs2_adjpval01_FC_cluster2)

DEG_Comp36LEC0vs2_adjpval01_FC_cluster2 = DEG_Comp36LEC0vs2_adjpval01_FC_cluster2 %>% dplyr::arrange(desc(avg_log2FC))

write.csv(DEG_Comp36LEC0vs2_adjpval01_FC_cluster2, "/Users/margl916/Documents/DEG_Comp36LEC0vs2_adjpval01_FC_cluster2.csv", row.names=TRUE)


#3. Which method to use (Agata or other/Anna-similar)

#USE OTHER METHOD!!!!! - gives more results !!! (But maybe more unspecific?)
#Seems that if other method used, similar results to 2way/individual cluster comparisons or multi way cluster comparison!!!!

#Other method:

#convert zebrafish gene symbols to Entrez gene ID's:
library(org.Dr.eg.db)
keytypes(org.Dr.eg.db)
library(clusterProfiler)


DEG_Comp36LEC0vs2_adjpval01_FC_cluster2_list = scan("/Users/margl916/Documents/DEG_Comp36LEC0vs2_adjpval01_FC_cluster2_list.txt", character(), quote = "")
t <- c(DEG_Comp36LEC0vs2_adjpval01_FC_cluster2_list)

et <- bitr(t, fromType="SYMBOL", toType=(c("ENTREZID","PATH","GO","ALIAS","GENENAME")), OrgDb="org.Dr.eg.db")
head(et)
write.csv(et, "/Users/margl916/Documents/et.csv", row.names=TRUE)

entrezID_et <- select(et, ENTREZID)

#remove duplicate rows across entire data frame 
entrezID_et_single <- entrezID_et%>%distinct(.keep_all = TRUE)
#save
write.csv(entrezID_et_single, "/Users/margl916/Documents/entrezID_et_single.csv", row.names=TRUE)
#scan as text file only one column - entrezID_et_single
entrezID_et_single_list = scan("/Users/margl916/Documents/entrezID_et_single_list.txt", character(), quote = "")


#4. Choose pvalueCutoff 

ans.go <- enrichGO(gene = entrezID_et_single_list, ont = "BP",
                   OrgDb ="org.Dr.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.5)

#or other gseaPvalCO 0.01 or 0.1 or 1

tab.go <- as.data.frame(ans.go)
View(tab.go)
write.csv(tab.go, "/Users/margl916/Documents/tab.go.csv", row.names=TRUE)
#RENAME FILE ACCORDINGLY

#visualization
library(enrichplot)
dotplot(ans.go, showCategory=20) + ggtitle("GO BP") +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1)
#file save as 600x600

upsetplot(ans.go, showCategory=20) + ggtitle("GO BP")


#Agata method: 
  
avg_log2FC_vect=DEG_Multi36_adjpval01_FC_cluster0$avg_log2FC
names(avg_log2FC_vect)=DEG_Multi36_adjpval01_FC_cluster0$gene 
head(avg_log2FC_vect)

#4. Choose pvalueCutoff 

gseaPvalCO=0.5
gsea.GO.res = gseGO(geneList=avg_log2FC_vect, 
                    ont ="BP", 
                    keyType = "SYMBOL",
                    nPermSimple = 500000,
                    minGSSize = 3, 
                    maxGSSize = 800,
                    eps=0,
                    pvalueCutoff = gseaPvalCO, 
                    verbose = TRUE, 
                    OrgDb = org.Dr.eg.db, 
                    pAdjustMethod = "BH",
                    scoreType = "pos")

#or other gseaPvalCO 0.01 or 0.1 or 1
#or use scoreType = “std” (positive and negative enriched)


gsea.GO.summary=gsea.GO.res@result%>%dplyr::select(Description,ID,NES,p.adjust)
head(gsea.GO.summary)
library(kableExtra)
kable(gsea.GO.summary %>% dplyr::slice_head(n=15), booktabs = TRUE, row.names = FALSE, caption = "Results of GO term GSEA.")%>% 
  kableExtra::kable_minimal(full_width = FALSE)


library(enrichplot, quietly = TRUE)
library(viridis, quietly = TRUE)

gsea_for_plot=gsea.GO.res@result %>% dplyr::select(Description,ID,NES,p.adjust) %>%
  dplyr::arrange(desc(abs(NES))) %>%
  dplyr::slice_head(n=15)

categories_dot=gsea_for_plot$Description

dotplot(gsea.GO.res, showCategory=categories_dot, orderBy = "NES", label_format=50) +
  theme(text = element_text(size = 10)) + 
  theme(axis.text.y = element_text(size=10)) +  
  viridis::scale_fill_viridis(direction = -1)


#using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).


#NOT adj_p_val ordered (geht so nicht - error, order muss absteigend sein, geht nicht bei adj_p_val)
  