
#RAID v2.0 all rna-rna interactions in human
setwd("/pool1/zhangyan/projects/CAD/data")
raid<-read.table("raid.v2_rna-rna.txt", sep="\t", header=T)
raid.human<-raid[intersect(which(raid[,"Species1"] %in% "Homo sapiens"), which(raid[,"Species2"] %in% "Homo sapiens")),-c(1,12)]
raid.human.score<-raid.human[which(raid.human[,"Score"]>0.3),]
raid.human.score.miRNA_circRNA<-raid.human.score[intersect(which(raid.human.score[,"Category1"]=="miRNA"), which(raid.human.score[,"Category2"]=="circRNA")),]
raid.human.score.miRNA_mRNA<-raid.human.score[intersect(which(raid.human.score[,"Category1"]=="miRNA"), which(raid.human.score[,"Category2"]=="mRNA")),]
raid.human.score.miRNA_circRNA_mRNA<-merge(raid.human.score.miRNA_circRNA, raid.human.score.miRNA_mRNA, by.x="Interactor1", by.y="Interactor1")
dim(raid.human.score.miRNA_circRNA_mRNA) #[1] 4297346      21

miRNA_circRNA_mRNA_uniq<-unique(raid.human.score.miRNA_circRNA_mRNA[,c("Interactor1","Interactor2.x","Interactor2.y")])
circRNA_mRNA_sep<-apply(miRNA_circRNA_mRNA_uniq, 1, function(x){paste0(x[2],"_",x[3])})
miRNA_circRNA_mRNA_sep<-apply(miRNA_circRNA_mRNA_uniq, 1, function(x){paste0(x[2],"_",x[3],"_",x[1])})
miRNA_circRNA_mRNA_sep2<-cbind(circRNA_mRNA_sep,miRNA_circRNA_mRNA_sep)
miRNA_circRNA_mRNA_sep2_5<-miRNA_circRNA_mRNA_sep2[which(as.character(miRNA_circRNA_mRNA_sep2[,1]) %in% names(table(miRNA_circRNA_mRNA_sep2[,1]))[which(as.numeric(table(miRNA_circRNA_mRNA_sep2[,1]))>2)]),]
head(miRNA_circRNA_mRNA_sep2_5)
length(unique(miRNA_circRNA_mRNA_sep2_5[,2])) #2089725
length(unique(miRNA_circRNA_mRNA_sep2_5[,1])) #486136

write.table(miRNA_circRNA_mRNA_sep2_5,"/pool1/zhangyan/projects/CAD/2/miRNA_circRNA_mRNA_sep2_3.txt", sep="\t", quote=F, row.names=F, col.names=T)

#gwas variant genes and circs
TableS1<-read.table("/pool1/zhangyan/projects/CAD/2/TableS1_nature2024.txt", sep="\t", header=T)
length(unique(TableS1[,"Lead_SNP_rsID"])) #307
gwas.genes<-unique(TableS1[,"gene"]) #2720
gwas.circs_withcirc<-paste0("circRNA_",gwas.genes) #2720
gwas.circs_withcirc<-intersect(gwas.circs_withcirc, paste0("circRNA_",raid.human[which(raid.human[,"Category2"]=="circRNA"),"Interactor2"])) #90


#known genes/circs/mirnas
setwd("/pool1/zhangyan/projects/CAD/data/known")
#known
known_circs<-read.table("known_circs.txt", sep="\t", header=T)
known_circs_withcirc<-paste0("circRNA_",known_circs[,1]) #29
known_genes<-read.table("known_genes.txt", sep="\t", header=T)
known_genes<-unique(known_genes[,1]) #54
known_miRNAs<-read.table("known_miRNAs.txt", sep="\t", header=T)
known_miRNAs<-unique(known_miRNAs[,1]) #179

#DE genes
setwd("/pool1/zhangyan/projects/CAD/data/DE")
DE_genes<-read.csv("DE_genes.tsv", sep="\t", header=T)
DE_genes_pval_FC<-unique(DE_genes[intersect(which(DE_genes[,"adj.P.Val"]<0.05), which(abs(DE_genes[,"logFC"])>log(1.5,2))),"ORF"])
write.table(DE_genes_pval_FC, "/pool1/zhangyan/projects/CAD/2/DE_genes_pval_FC.txt", quote=F, sep="\t", row.names=F, col.names=F)

DE_genes_pval_FC_up<-unique(DE_genes[intersect(which(DE_genes[,"adj.P.Val"]<0.05), which(DE_genes[,"logFC"]>log(1.5,2))),"ORF"])
DE_genes_pval_FC_dw<-unique(DE_genes[intersect(which(DE_genes[,"adj.P.Val"]<0.05), which(DE_genes[,"logFC"]<(-log(1.5,2)))),"ORF"])
DE_genes_pval_FC_up<-setdiff(DE_genes_pval_FC_up, DE_genes_pval_FC_dw)
# write.table(DE_genes_pval_FC_up,"DE_genes_pval_FC_up.txt",quote=F, sep="\t",row.names=F, col.names=F)
DE_genes_pval_FC_dw<-setdiff(DE_genes_pval_FC_dw, DE_genes_pval_FC_up)
# write.table(DE_genes_pval_FC_dw,"DE_genes_pval_FC_dw.txt",quote=F, sep="\t",row.names=F, col.names=F)

#DE circs
DE_circRNAs<-read.csv("DE_circRNAs.tsv", sep="\t", header=T)
DE_circRNAs_pval_FC<-unique(DE_circRNAs[intersect(which(DE_circRNAs[,"adj.P.Val"]<0.05), which(abs(DE_circRNAs[,"logFC"])>log(1.5,2))),"ID"])
circ_GPL22121<-read.table("circ_GPL22121.txt",sep="\t",header=T)
DE_circRNAs_pval_FC_genesymbol<-unique(circ_GPL22121[which(circ_GPL22121[,"ID"] %in% DE_circRNAs_pval_FC),"GENE_SYMBOL"])
# write.table(DE_circRNAs_pval_FC_genesymbol,"DE_circRNAs_pval_FC_genesymbol.txt",quote=F, sep="\t",row.names=F, col.names=F)
DE_circRNAs_pval_FC_genesymbol_withcirc<-paste0("circRNA_",DE_circRNAs_pval_FC_genesymbol)

DE_circRNAs_pval_FC_up<-unique(DE_circRNAs[intersect(which(DE_circRNAs[,"adj.P.Val"]<0.05), which(DE_circRNAs[,"logFC"]>log(1.5,2))),"ID"])
DE_circRNAs_pval_FC_up_genesymbol<-unique(circ_GPL22121[which(circ_GPL22121[,"ID"] %in% DE_circRNAs_pval_FC_up),"GENE_SYMBOL"])
# write.table(DE_circRNAs_pval_FC_up_genesymbol,"DE_circRNAs_pval_FC_up_genesymbol.txt",quote=F, sep="\t",row.names=F, col.names=F)

DE_circRNAs_pval_FC_dw<-unique(DE_circRNAs[intersect(which(DE_circRNAs[,"adj.P.Val"]<0.05), which(DE_circRNAs[,"logFC"]<(-log(1.5,2)))),"ID"])
DE_circRNAs_pval_FC_dw_genesymbol<-unique(circ_GPL22121[which(circ_GPL22121[,"ID"] %in% DE_circRNAs_pval_FC_dw),"GENE_SYMBOL"])
# write.table(DE_circRNAs_pval_FC_dw_genesymbol,"DE_circRNAs_pval_FC_dw_genesymbol.txt",quote=F, sep="\t",row.names=F, col.names=F)


##manhaton plot
TableS1_nature2024_PoPSScore<-read.table("TableS1_nature2024_PoPSScore.txt", sep="\t", header=T)

TableS1_nature2024_PoPSScore_raid.circs<-TableS1_nature2024_PoPSScore[which(TableS1_nature2024_PoPSScore[,"gene"] %in% raid.human[which(raid.human[,"Category2"]=="circRNA"),"Interactor2"]),]
TableS1_nature2024_PoPSScore_raid.circs[,"gene"]<-paste0("circRNA_", TableS1_nature2024_PoPSScore_raid.circs[,"gene"])

TableS1_nature2024_PoPSScore_genesandcircs<-rbind(TableS1_nature2024_PoPSScore, TableS1_nature2024_PoPSScore_raid.circs)
dim(TableS1_nature2024_PoPSScore_genesandcircs) #[1] 3939    5

write.table(TableS1_nature2024_PoPSScore_genesandcircs,"TableS1_nature2024_PoPSScore_genesandcircs.txt",quote=F, sep="\t",row.names=F, col.names=T)


##CAD risk circ-mRNA network
#kept the pairs with gwas.genes or known or DE genes or cirrnas

# gwas.circs_withcirc
# gwas.genes

# known_circs_withcirc
# known_genes
# known_miRNAs

# DE_circRNAs_pval_FC_genesymbol_withcirc
# DE_genes_pval_FC


##VNN plot
setwd("/pool1/zhangyan/projects/CAD/2")

# install.packages("VennDiagram")
library(VennDiagram)

## Define your gene lists
set1 <- gwas.circs_withcirc
set2 <- known_circs_withcirc
set3 <- DE_circRNAs_pval_FC_genesymbol_withcirc

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2, Set3 = set3),
  category.names = c("GWAS circRNAs", "known circRNAs", "DE circRNAs"),
  filename = NULL,  # Draw on screen
  col = c("red", "blue", "green"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue", "green")
)

# Display the plot
pdf("venn_circRNAs.pdf")
grid.draw(venn.plot)
dev.off()

intersect(known_circs_withcirc,gwas.circs_withcirc)
intersect(known_circs_withcirc,DE_circRNAs_pval_FC_genesymbol_withcirc)

## Define your gene lists
set1 <- gwas.genes
set2 <- known_genes
set3 <- DE_genes_pval_FC

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(Set1 = set1, Set2 = set2, Set3 = set3),
  category.names = c("GWAS genes", "known genes", "DE genes"),
  filename = NULL,  # Draw on screen
  col = c("red", "blue", "green"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue", "green")
)

# Display the plot
pdf("venn_genes.pdf")
grid.draw(venn.plot)
dev.off()


setwd("/pool1/zhangyan/projects/CAD/2")

triplets<-miRNA_circRNA_mRNA_sep2_5[,c(2,2,2)]
triplets[,1]<-sapply(strsplit(miRNA_circRNA_mRNA_sep2_5[,2],"_"), function(x){x[1]})
triplets[,2]<-sapply(strsplit(miRNA_circRNA_mRNA_sep2_5[,2],"_"), function(x){x[2]})
triplets[,3]<-sapply(strsplit(miRNA_circRNA_mRNA_sep2_5[,2],"_"), function(x){x[3]})
triplets[,1]<-paste0("circRNA_",triplets[,1])
write.table(triplets, "triplets.txt", sep="\t",quote=F,row.names=F,col.names=F)

node_attri<-cbind(c(unique(triplets[,1]), unique(triplets[,2]), unique(triplets[,3])), c(rep("circ",length(unique(triplets[,1]))), rep("mrna",length(unique(triplets[,2]))), rep("mirna",length(unique(triplets[,3])))))
write.table(node_attri, "node_attri.txt", sep="\t",quote=F,row.names=F,col.names=F)

triplets_interactions<-rbind(unique(triplets[,c(1,3)]),unique(triplets[,c(3,2)]))
write.table(triplets_interactions, "triplets_interactions.txt", sep="\t", quote=F, row.names=F, col.names=F)
triplets_node_catolog<-rbind(rbind(cbind(unique(triplets[,1]),"circ"), cbind(unique(triplets[,2]),"gene")), cbind(unique(triplets[,3]),"mirna"))
write.table(triplets_node_catolog, "triplets_node_catolog.txt", sep="\t", quote=F, row.names=F, col.names=F)


#circ-mrna-mirna network
circs_CAD_withcirc<-unique(c(gwas.circs_withcirc, known_circs_withcirc, DE_circRNAs_pval_FC_genesymbol_withcirc))
genes_CAD<-unique(c(gwas.genes, known_genes, DE_genes_pval_FC))

triplets_CAD<-triplets[union(which(triplets[,1] %in% circs_CAD_withcirc), which(triplets[,1] %in% genes_CAD)),]
dim(triplets_CAD) #[1] 306559     3
write.table(triplets_CAD,"./triplets_CAD.txt",sep="\t",quote=F,row.names=F,col.names=F)
length(unique(triplets_CAD[,1])) #69
length(unique(triplets_CAD[,2])) #8540
length(unique(triplets_CAD[,3])) #230

#circ-mrna network
circ_mrna<-unique(triplets_CAD[,1:2])
length(unique(circ_mrna[,1])) #69
length(unique(circ_mrna[,2])) #8540
dim(circ_mrna) #[1] 74307   2
write.table(circ_mrna,"./circ_mrna_CAD.txt",sep="\t",quote=F,row.names=F,col.names=F)
degree_circ_pair<-as.data.frame(table(circ_mrna[,1]))
degree_mrna_pair<-as.data.frame(table(circ_mrna[,2]))
degree_circ_mrna<-rbind(degree_circ_pair, degree_mrna_pair)
write.table(degree_circ_mrna,"./circ_mrna_CAD_degree_circ_mrna.txt",sep="\t",quote=F,row.names=F,col.names=F)

intersect(degree_circ_pair[order(degree_circ_pair[,"Freq"],decreasing=T),1], known_circs_withcirc)

#GO analysis
library(clusterProfiler)
library("org.Hs.eg.db")

circ_mrna.geneid <- bitr(unique(circ_mrna[,2]), fromType="SYMBOL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)[,2]
circ_mrna.go<-enrichGO(circ_mrna.geneid, OrgDb="org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
head(as.data.frame(circ_mrna.go))
write.table(as.data.frame(circ_mrna.go), "circ_mrna.go.txt", sep="\t", quote=F, row.names=F, col.names=T)

circ_mrna.go_sub <- circ_mrna.go[circ_mrna.go$ID %in% c("GO:0042060","GO:0016055","GO:0008380","GO:2001020","GO:0006282","GO:1901987","GO:0010506","GO:2001233","GO:0002520","GO:0048041","GO:0003158","GO:0006338","GO:0198738"),]
pdf("circ_mrna.go_sub.pdf")
dotplot(circ_mrna.go_sub, title="circ_mrna.GO", font.size=10)
dev.off()


##13 EC-specific programs from nature 2024
EC.geneset <-read.table("/pool1/zhangyan/projects/CAD/2/TableS14.txt", sep="\t", header=T)

# EC.geneset_8<-EC.geneset[which(EC.geneset[,1] %in% "K60_8"),2]
# circ_mrna_8<-circ_mrna[which(as.character(circ_mrna[,2]) %in% EC.geneset_8),]
# length(unique(circ_mrna_8[,1]))#[1] 61
# length(unique(circ_mrna_8[,2]))#[1] 164

EC.geneset.all_2<-EC.geneset[which(EC.geneset[,1] %in% c("K60_8","K60_50","K60_48","K60_47","K60_41","K60_39","K60_35","K60_31","K60_29","K60_28","K60_25","K60_2","K60_15")),]
write.table(EC.geneset.all_2,"./TableS14_EC.txt",sep="\t",quote=F,row.names=F,col.names=F)

#heatmap
EC.genelist<-unique(EC.geneset.all_2[,"Gene"])
EC_matrix<-matrix(0,nrow=length(EC.genelist),ncol=13)
rownames(EC_matrix)<-EC.genelist
colnames(EC_matrix)<-unique(EC.geneset.all_2[,1])
for (i in 1:13){
    glist<-EC.geneset.all_2[which(EC.geneset.all_2[,1] %in% colnames(EC_matrix)[i]),2]
    EC_matrix[which(rownames(EC_matrix) %in% glist),i]<-1
}

library(pheatmap)
pdf("pheatmap_EC_matrix.pdf")
pheatmap(EC_matrix, cluster_cols=F, cluster_rows=F, color = colorRampPalette(c("white", "navy"))(50), show_rownames = F)
dev.off()


EC.geneset.all<-EC.geneset[which(EC.geneset[,1] %in% c("K60_8","K60_50","K60_48","K60_47","K60_41","K60_39","K60_35","K60_31","K60_29","K60_28","K60_25","K60_2","K60_15")),2]
circ_mrna_EC<-circ_mrna[which(as.character(circ_mrna[,2]) %in% EC.geneset.all),]
write.table(circ_mrna_EC,"./circ_mrna_EC.txt",sep="\t",quote=F,row.names=F,col.names=F)
dim(circ_mrna_EC) #[1] 13778     2
length(unique(circ_mrna_EC[,1]))#[1] 68
length(unique(circ_mrna_EC[,2]))#[1] 1333

degree_circ_mrna_EC<-rbind(as.data.frame(table(circ_mrna_EC[,1])),as.data.frame(table(circ_mrna_EC[,2])))
write.table(degree_circ_mrna_EC,"./degree_circ_mrna_EC.txt",sep="\t",quote=F,row.names=F,col.names=F)

circ_mrna_EC_degree.10.gwas<-intersect(intersect(as.character(circ_mrna_EC[,1]),gwas.circs_withcirc), as.data.frame(table(circ_mrna_EC[,1]))[which(as.data.frame(table(circ_mrna_EC[,1]))[,2]>100),1])
circ_mrna_EC_degree.10.gwas_net <- circ_mrna_EC[which(circ_mrna_EC[,1] %in% circ_mrna_EC_degree.10.gwas),]
length(unique(circ_mrna_EC_degree.10.gwas_net[,1])) #37
length(unique(circ_mrna_EC_degree.10.gwas_net[,2])) #1323

#triplets_EC
miRNA_circRNA_mRNA_sep2_5_circ<-miRNA_circRNA_mRNA_sep2_5
miRNA_circRNA_mRNA_sep2_5_circ[,1]<-paste0("circRNA_", miRNA_circRNA_mRNA_sep2_5[,1])
miRNA_circRNA_mRNA_sep2_5_circ[,2]<-paste0("circRNA_", miRNA_circRNA_mRNA_sep2_5[,2])

circ_mrna_EC_paste<-as.character(apply(circ_mrna_EC, 1, function(x){paste0(x[1],"_",x[2])}))
circ_mrna_EC_paste_triplets<-as.data.frame(miRNA_circRNA_mRNA_sep2_5_circ[which(miRNA_circRNA_mRNA_sep2_5_circ[,1] %in% circ_mrna_EC_paste),2])

triplets_EC<-cbind(sapply(strsplit(circ_mrna_EC_paste_triplets[,1],"_"), function(x){paste0(x[1],"_",x[2])}), sapply(strsplit(circ_mrna_EC_paste_triplets[,1],"_"), function(x){x[3]}), sapply(strsplit(circ_mrna_EC_paste_triplets[,1],"_"), function(x){x[4]}))

triplets_EC[,3]<-gsub("-5p","",triplets_EC[,3])
triplets_EC_knownmirna<-triplets_EC[which(triplets_EC[,3] %in% known_miRNAs),] #[1] 22539     3
dim(unique(triplets_EC_knownmirna[,1:2])) #[1] 11150     2
write.table(triplets_EC_knownmirna,"./triplets_EC_knownmirna.txt",sep="\t",quote=F,row.names=F,col.names=F)

triplets_EC_knownmirna_degree.10.gwas<-triplets_EC_knownmirna[which(triplets_EC_knownmirna[,1] %in% circ_mrna_EC_degree.10.gwas),]
dim(triplets_EC_knownmirna_degree.10.gwas) #[1] 21180     3
length(unique(triplets_EC_knownmirna_degree.10.gwas[,1])) #circs 35
length(unique(triplets_EC_knownmirna_degree.10.gwas[,2])) #genes 1252
length(unique(triplets_EC_knownmirna_degree.10.gwas[,3])) #mirnas 54
dim(unique(triplets_EC_knownmirna_degree.10.gwas[,1:2])) #[1] 10228     2

#common genes involved in 13 EC-specific programs (freq>5)
EC.geneset <-read.table("/pool1/zhangyan/projects/CAD/2/TableS14.txt", sep="\t", header=T)
EC.geneset.all<-EC.geneset[which(EC.geneset[,1] %in% c("K60_8","K60_50","K60_48","K60_47","K60_41","K60_39","K60_35","K60_31","K60_29","K60_28","K60_25","K60_2","K60_15")),2]
EC.geneset.all.freq<-as.data.frame(table(EC.geneset.all))

table(EC.geneset.all.freq[,2])
# Var1	Freq
# 1	1140
# 2	440
# 3	269
# 4	136
# 5	56
# 6	20
# 7	16
# 8	1
# 9	1

EC.geneset.all.freq.common<-as.character(EC.geneset.all.freq[which(EC.geneset.all.freq[,2]>5),1])
triplets_EC_knownmirna_commongenes.in.EC<-triplets_EC_knownmirna[which(triplets_EC_knownmirna[,2] %in% EC.geneset.all.freq.common),]
write.table(triplets_EC_knownmirna_commongenes.in.EC,"./triplets_EC_knownmirna_commongenes.in.EC.txt",sep="\t",quote=F,row.names=F,col.names=F)

length(unique(triplets_EC_knownmirna_commongenes.in.EC[,1])) #44
length(unique(triplets_EC_knownmirna_commongenes.in.EC[,2])) #26
length(unique(triplets_EC_knownmirna_commongenes.in.EC[,3])) #53

triplets_EC_knownmirna_commongenes.in.EC_net<-rbind(unique(triplets_EC_knownmirna_commongenes.in.EC[,c(1,3)]), unique(triplets_EC_knownmirna_commongenes.in.EC[,c(3,2)]))
write.table(triplets_EC_knownmirna_commongenes.in.EC_net,"./triplets_EC_knownmirna_commongenes.in.EC_net.txt",sep="\t",quote=F,row.names=F,col.names=F)


##circRNA_ZNF609
#gwas and known
intersect(unique(triplets_EC_knownmirna_commongenes.in.EC[,1]),intersect(gwas.circs_withcirc,known_circs_withcirc))
# [1] "circRNA_ZNF609"  "circRNA_HERPUD2"
intersect(unique(triplets_EC_knownmirna_commongenes.in.EC[,1]),gwas.circs_withcirc)
intersect(unique(triplets_EC_knownmirna_commongenes.in.EC[,1]),known_circs_withcirc)

triplets_EC_knownmirna_commongenes.in.EC[which(triplets_EC_knownmirna_commongenes.in.EC[,1] %in% "circRNA_ZNF609"),]
#  [1,] "circRNA_ZNF609" "SPARC"  "hsa-miR-15a"
#  [2,] "circRNA_ZNF609" "HSPG2"  "hsa-miR-15a"
#  [3,] "circRNA_ZNF609" "MGAT4A" "hsa-miR-15a"
#  [4,] "circRNA_ZNF609" "APLN"   "hsa-miR-15a"
#  [5,] "circRNA_ZNF609" "SPARC"  "hsa-miR-15b"
#  [6,] "circRNA_ZNF609" "HSPG2"  "hsa-miR-15b"
#  [7,] "circRNA_ZNF609" "APLN"   "hsa-miR-15b"
#  [8,] "circRNA_ZNF609" "APP"    "hsa-miR-15b"
#  [9,] "circRNA_ZNF609" "MGAT4A" "hsa-miR-16" 
# [10,] "circRNA_ZNF609" "SPARC"  "hsa-miR-16" 
# [11,] "circRNA_ZNF609" "APP"    "hsa-miR-16" 


##circRNA_ABCC1
#circRNA_ABCC1 with the highest degree
names(table(triplets_EC_knownmirna_commongenes.in.EC[,1]))[as.numeric(table(triplets_EC_knownmirna_commongenes.in.EC[,1]))>25]
# [1] "circRNA_ABCC1"    "circRNA_KIAA1586" "circRNA_PIGU"

triplets_EC_knownmirna_commongenes.in.EC_circABCC1<-unique(triplets_EC_knownmirna_commongenes.in.EC[which(triplets_EC_knownmirna_commongenes.in.EC[,1] %in% "circRNA_ABCC1"),])
dim(triplets_EC_knownmirna_commongenes.in.EC_circABCC1) #[1] 27  3
#  [1,] "circRNA_ABCC1" "MGAT4A"  "hsa-let-7b"  
#  [2,] "circRNA_ABCC1" "COL5A1"  "hsa-let-7b"  
#  [3,] "circRNA_ABCC1" "DUSP1"   "hsa-let-7b"  
#  [4,] "circRNA_ABCC1" "MGAT4A"  "hsa-let-7i"  
#  [5,] "circRNA_ABCC1" "APP"     "hsa-let-7i"  
#  [6,] "circRNA_ABCC1" "DUSP1"   "hsa-let-7i"  
#  [7,] "circRNA_ABCC1" "CALD1"   "hsa-let-7i"  
#  [8,] "circRNA_ABCC1" "MGAT4A"  "hsa-miR-15a" 
#  [9,] "circRNA_ABCC1" "APLN"    "hsa-miR-15a" 
# [10,] "circRNA_ABCC1" "MGAT4A"  "hsa-miR-16"  
# [11,] "circRNA_ABCC1" "APP"     "hsa-miR-16"  
# [12,] "circRNA_ABCC1" "CDH2"    "hsa-miR-18a" 
# [13,] "circRNA_ABCC1" "APP"     "hsa-miR-18a" 
# [14,] "circRNA_ABCC1" "MTUS1"   "hsa-miR-34a" 
# [15,] "circRNA_ABCC1" "MARCKS"  "hsa-miR-34a" 
# [16,] "circRNA_ABCC1" "MGAT4A"  "hsa-miR-34a" 
# [17,] "circRNA_ABCC1" "ZFP36L1" "hsa-miR-34a" 
# [18,] "circRNA_ABCC1" "MGAT4A"  "hsa-miR-34c" 
# [19,] "circRNA_ABCC1" "CDH2"    "hsa-miR-34c" 
# [20,] "circRNA_ABCC1" "MARCKS"  "hsa-miR-34c" 
# [21,] "circRNA_ABCC1" "CALD1"   "hsa-miR-34c" 
# [22,] "circRNA_ABCC1" "APP"     "hsa-miR-4306"
# [23,] "circRNA_ABCC1" "MTUS1"   "hsa-miR-98"  
# [24,] "circRNA_ABCC1" "APP"     "hsa-miR-98"  
# [25,] "circRNA_ABCC1" "MGAT4A"  "hsa-miR-98"  
# [26,] "circRNA_ABCC1" "DUSP1"   "hsa-miR-98"  
# [27,] "circRNA_ABCC1" "ZFP36L1" "hsa-miR-98"  
length(unique(triplets_EC_knownmirna_commongenes.in.EC_circABCC1[,1])) #1
length(unique(triplets_EC_knownmirna_commongenes.in.EC_circABCC1[,2])) #10
length(unique(triplets_EC_knownmirna_commongenes.in.EC_circABCC1[,3])) #9

DE_circRNAs_ABCC1<-DE_circRNAs[which(DE_circRNAs[,1] %in% circ_GPL22121[which(circ_GPL22121[,"GENE_SYMBOL"]=="ABCC1"),"ID"]),] #no diff expressed
write.table(DE_circRNAs_ABCC1,"./DE_circRNAs_ABCC1.txt",sep="\t",quote=F,row.names=F,col.names=T)


save.image(file="/pool1/zhangyan/projects/CAD/2/res.RData")
load(file="/pool1/zhangyan/projects/CAD/2/res.RData")


##circ_mrna_CAD_top200genes
setwd("/pool1/zhangyan/projects/CAD/2")
circ_mrna_CAD_degree_circ_mrna<-read.table("circ_mrna_CAD_degree_circ_mrna.txt", header=F, sep="\t")
#top 200 genes
circ_mrna_CAD_degree_top200genes<-circ_mrna_CAD_degree_circ_mrna[which(circ_mrna_CAD_degree_circ_mrna[,2]>26),1]

circ_mrna_CAD<-read.table("circ_mrna_CAD.txt", sep="\t", header=F)
circ_mrna_CAD_top200genes<-circ_mrna_CAD[which(circ_mrna_CAD[,2] %in% circ_mrna_CAD_degree_top200genes),]
write.table(circ_mrna_CAD_top200genes,"circ_mrna_CAD_top200genes.txt", sep="\t", quote=F,row.names=F, col.names=F)

