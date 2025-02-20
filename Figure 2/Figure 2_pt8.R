setwd("C:/MD_Anderson/4.PDAC_ST/6.SpatialInferCNV_by_patient_manual_anno/New_inferCNV_core_4/InferCNV_outputs_patient2_8")

setwd("C:/MD_Anderson/4.PDAC_ST_github/PDAC_Figures")
library(ape)
library(phylogram)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialInferCNV)
library(tidyverse)
library(ggtree)

#Use read.dendrogram() to import the dendogram file
input_path <- "C:/MD_Anderson/4.PDAC_ST/6.SpatialInferCNV_by_patient_manual_anno/New_inferCNV_core_4/InferCNV_outputs_patient2_8/" 

dendrogram <- read.dendrogram(file = paste0(input_path, "infercnv.observations_dendrogram.txt"))
mytree <- as.phylo(dendrogram)
my.subtrees <- subtrees(mytree)

#Generate a phylogenetic tree with node IDs for downstream manual clone assignment
png(paste0("./Figure 2_pt8_phylo.png"), width = 3500, height = 2500, res = 300)
plot(mytree, show.tip.label = FALSE)
nodelabels(text = 1:mytree$Nnode, node= 1:mytree$Nnode + Ntip(mytree))
dev.off()
	
#===================================================================
library(SpatialInferCNV)
library(tidyverse)

#Load the seurat object 
st <- load("C:/MD_Anderson/4.PDAC_ST/6.SpatialInferCNV_by_patient_manual_anno/New_inferCNV_core_4/ST_case_2_p85.Rdata")

Node1 <- SelectingSubTreeData(my.subtrees, 3)			
Node2 <- SelectingSubTreeData(my.subtrees, 698)
Node3 <- SelectingSubTreeData(my.subtrees, 1091)
Node4 <- SelectingSubTreeData(my.subtrees, 1842)
Node5 <- SelectingSubTreeData(my.subtrees, 2298)								  
Node6 <- SelectingSubTreeData(my.subtrees, 3113)											  				
Node7 <- SelectingSubTreeData(my.subtrees, 3395)											  		
Node8 <- SelectingSubTreeData(my.subtrees, 3665)
Node9 <- SelectingSubTreeData(my.subtrees, 4663)					
Node10 <- SelectingSubTreeData(my.subtrees, 5551)		
Node11 <- SelectingSubTreeData(my.subtrees, 6291)										
Node12 <- SelectingSubTreeData(my.subtrees, 6997)				
Node13 <- SelectingSubTreeData(my.subtrees, 7499)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2, Node3, Node4, Node5, Node6, Node7, Node8, Node9, Node10, Node11, Node12, Node13)

#Then drop the Histology column
colnames(MergedNodes) = c("Barcode", "Node")

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_3", "Clone A1",
                     ifelse(MergedNodes$Node == "Node_698" , "Clone A4",
					 ifelse(MergedNodes$Node == "Node_1091" , "Clone A3",
					 ifelse(MergedNodes$Node == "Node_1842" , "Clone A2", 					 
					 ifelse(MergedNodes$Node == "Node_2298" , "Clone B2", 
					 ifelse(MergedNodes$Node == "Node_3113" , "Clone B1", 					 
					 ifelse(MergedNodes$Node == "Node_3395" , "Clone C1",
					 ifelse(MergedNodes$Node == "Node_3665" , "Clone C2",					 
					 ifelse(MergedNodes$Node == "Node_4663" , "Clone D1",
					 ifelse(MergedNodes$Node == "Node_5551" , "Clone D2",				 
					 ifelse(MergedNodes$Node == "Node_6291" , "Clone E2",
					 ifelse(MergedNodes$Node == "Node_6997" , "Clone E1",
                     ifelse(MergedNodes$Node == "Node_7499" , "Clone E3", MergedNodes$Node)))))))))))))


#Include clone information in the metadata.
ST_case_2@meta.data$Histology_anno_clone <- ST_case_2@meta.data$Histology_anno <- "NT"
ST_case_2@meta.data$Histology_anno_clone[match(MergedNodes$Barcode, rownames(ST_case_2@meta.data))] <- MergedNodes$Node
ST_case_2@meta.data$Clone <- ST_case_2@meta.data$Histology_anno_clone 
ST_case_2@meta.data$Clone[which(ST_case_2@meta.data$CNV_cluster == "NT")] <- "NT"

#--------------------------------------------------
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(yarrr)
library(scales)

#Information update
ST_case_2@meta.data$orig.ident <- gsub("Sample85_", "Pt-8", ST_case_2@meta.data$orig.ident)
ST_case_2@meta.data$site[ST_case_2@meta.data$site == "Panc"] <- "Pri"
ST_case_2@meta.data$site[ST_case_2@meta.data$site == "Liver"] <- "LiM"
ST_case_2@meta.data$site[ST_case_2@meta.data$site == "Lung"] <- "LuM"
ST_case_2@meta.data$site[ST_case_2@meta.data$site == "PC"] <- "PerM"
ST_case_2@meta.data$site[ST_case_2@meta.data$orig.ident == "Pt-8A"] <- "Pri"
ST_case_2@meta.data$site[ST_case_2@meta.data$orig.ident == "Pt-8B"] <- "LiM-1"
ST_case_2@meta.data$site[ST_case_2@meta.data$orig.ident == "Pt-8C"] <- "LiM-2"
ST_case_2@meta.data$site[ST_case_2@meta.data$orig.ident == "Pt-8D"] <- "PerM"
names(ST_case_2@images) <- c(names(table(paste(ST_case_2@meta.data$orig.ident, ST_case_2@meta.data$site))))
Idents(ST_case_2) <- factor(ST_case_2@meta.data$Clone, levels = sort(names(table(ST_case_2@meta.data$Clone))) )

#Define a color panel for tumor clone
mycol = piratepal(palette = "google")
mycol5 <- c(piratepal("basel", trans = 0.3))
mycol6 <- c(piratepal("pony", trans = 0.3))
		
hue_pal()(9)
c0_color = colorRampPalette(c(mycol5[3], "white"))
c1_color = colorRampPalette(c("#ef3b2c", "white"))
c2_color = colorRampPalette(c(mycol5[8], "white"))
c3_color = colorRampPalette(c(mycol6[5], "white"))
c4_color = colorRampPalette(c("#41b6c4", "white"))

clone.color <- table(ST_case_2@meta.data$Clone)							
c0_col <- c0_color(5 + 2)[0:length(grep("Clone A", names(clone.color)))]
c1_col <- c1_color(5 + 2)[0:length(grep("Clone B", names(clone.color)))]
c2_col <- c2_color(5 + 2)[0:length(grep("Clone C", names(clone.color)))]
c3_col <- c3_color(5 + 2)[0:length(grep("Clone D", names(clone.color)))]
c4_col <- c4_color(3 + 2)[0:length(grep("Clone E", names(clone.color)))]
nt_col <- "white"
	
clone.color <- c(c0_col, c1_col, c2_col, c3_col, c4_col, nt_col)
names(clone.color) <- c(names(table(ST_case_2$Clone)))

#Define a color panel for tumor site
site.color <- c((brewer.pal(9, "Blues"))[8], brewer.pal(9, "Oranges")[8], brewer.pal(9, "RdPu")[6], brewer.pal(9, "YlGnBu")[5])
names(site.color) <- c("Pri", "LiM", "LuM", "PerM")
site.color[5:6] <- colorRampPalette(c(brewer.pal(9, "Oranges")[8], "white"))(3)[1:2]
names(site.color)[c(5,6)] <- c("LiM-1", "LiM-2")
site.color[7:8] <- colorRampPalette(c(brewer.pal(9, "YlGnBu")[5], "white"))(3)[1:2]
names(site.color)[c(7,8)] <- c("PerM-1", "PerM-2")

#Clustering & flip base on the clone & Manually verify the definition of the tumor clone.
allnode = mytree$tip.label
cluster = MergedNodes$Node
names(cluster) = MergedNodes$Barcode

groupInfo <- lapply(levels(as.factor(cluster)), function(cc, cluster, allnode){
					return(names(cluster)[cluster == cc])
					},	cluster, allnode)					 
names(groupInfo) <- names(table(cluster))

MPT = groupOTU(mytree, groupInfo)
gp1 <- ggtree(MPT, aes(color = group)) + theme_tree() + geom_tippoint(alpha = 0.7)
gp1 <- gp1 + scale_color_manual(values=clone.color[grep("Clone", names(clone.color))])
gp1 + geom_text(aes(label = node))
gp1 <- flip(gp1, 7909, 10203)
gp1 + geom_text(aes(label = node))
gp1 <- flip(gp1, 10204, 11300)
gp1 + geom_text(aes(label = node))
gp1 <- flip(gp1, 11301, 12568)
gp1 + geom_text(aes(label = node))
gp1 <- flip(gp1, 7910, 8604)
gp1 + geom_text(aes(label = node))
gp1 <- flip(gp1, 8605, 8997)
gp1 <- flip(gp1, 9749, 8998)
gp1 <- flip(gp1, 11302, 11572)
gp1 <- flip(gp1, 14198, 14903)

#Calculate the tissue composition for each subclone
row_anno_clean <- ST_case_2@meta.data[,c("site", "Clone")]
row_anno_clean$site <- factor(row_anno_clean$site, levels = c("Pri", "LiM", "LiM-1", "LiM-2", "LuM", "PerM", "PerM-1", "PerM-2"))
row_anno_clean <- row_anno_clean[grep("Clone", row_anno_clean$Clone),]
clone_composition <- table(row_anno_clean$site, row_anno_clean$Clone)
clone_composition <- t(clone_composition[,grep("Clone", colnames(clone_composition))])

clone_composition_df <- data.frame(matrix(as.vector(clone_composition), ncol = ncol(clone_composition)))
rownames(clone_composition_df) <- rownames(clone_composition)
colnames(clone_composition_df) <- colnames(clone_composition)
gp1 + geom_text(aes(label=node), hjust=0)

#Manually specify the parent node of each tumor clone to add the pie chart.
clone_composition_df$node <-  c(13654, 13291, 13182, 12790,  12097, 9267, 8023, 15294, 15032,  14343, 10877, 10732, 9897)

pie_df <- nodepie(clone_composition_df, cols=1:8, alpha=0.9, color = site.color ) 

gp_combo1 <- inset(gp1, pie_df, height=0.1, width = 0.1)
gp_combo2 <- inset(gp1, pie_df, height=0.1, width = 0.1, hjust = -0.035, vjust = 0.03)

dd = subset(gp_combo1[[1]], isTip)

#Get the order of spot for downstream heatmap
spot_order1 <- dd$label[order(dd$y, decreasing=TRUE)]

#Load the infercnv output
data1 = read.delim("C:/MD_Anderson/4.PDAC_ST/6.SpatialInferCNV_by_patient_manual_anno/New_inferCNV_core_4/InferCNV_outputs_patient2_8/infercnv.observations.txt", sep = " ", row.names = 1)
colnames(data1) <- gsub("\\.", "-", colnames(data1))
cnv <- t(data1)
data <- cnv[MergedNodes[grep("Clone", (MergedNodes$Node)),]$Barcode,]

#Remove the genes on chromosomes X and Y.
col_anno <- read.delim("C:/MD_Anderson/4.PDAC_ST/GeneOrderFile.txt", head = F)
col_anno <- col_anno[col_anno$V1 %in% colnames(cnv),]
col_anno <- col_anno[col_anno$V2 %in% paste0("chr", 1:22),]
rownames(col_anno) <- col_anno[,1]

col_anno_clean <- data.frame(col_anno[,1:2], row.names = rownames(col_anno))
intersect_gene <- intersect(colnames(cnv), rownames(col_anno_clean))
cnv_hdata_clean <- cnv[, intersect_gene]
col_anno_clean <- col_anno_clean[rownames(col_anno_clean) %in% intersect_gene,]
col_anno_clean <- data.frame(Chromosome = col_anno_clean[,2], row.names = rownames(col_anno_clean))
col_anno_clean$Chromosome <- factor(col_anno_clean$Chromosome, levels = paste0("chr", c(1:22)))

tb = table(col_anno_clean);
tb = tb[ paste0("chr", 1:22)]
seg.chr = cumsum(tb)

#Create color panel for heatmap
library(pheatmap)
library(pals) 
pals::cols25()[1:22] -> chr.color
names(chr.color) = paste0("chr", c(1:22))

aCol = list(Chromosome = chr.color, site = site.color, Clone = clone.color)
data_clean <- data[,rownames(col_anno_clean)]

#Generate CNV heatmap
mid_chr <- cumsum(tb)			
mid_tb <- round(tb/2)			
mid_chr_mid <- round(mid_chr/2)
mid_chr_mid[1] <- round(mid_chr[1]/2)
mid_chr_mid[2:length(mid_chr_mid)] <- round(mid_chr[-1] - mid_tb[-1])

data_clean_map <- data_clean

for (k in 1:length(colnames(data_clean_map))){
	colnames(data_clean_map)[k] <- " "
}
colnames(data_clean_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))

pheatmap(data_clean_map[spot_order1,], scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), angle_col = "90", 
         cluster_cols = F, gaps_row = NULL, gaps_col = seg.chr, annotation_colors = aCol, fontsize = 8,
		 annotation_row = row_anno_clean,
         annotation_col = NA, 
         cluster_rows = F, show_colnames = T, show_rownames = F, border_color = NA,
		 width = 10,
         height = 5,
		 cutree_col = 2,
		 file = paste0("Figure 2_pt8_heatmap_raw.png") 	)

#Merge the CNV from spot to clone level
dim(data)
dim(MergedNodes)

MergedNodes_clean <- MergedNodes[ MergedNodes$Barcode %in% rownames(data), ]
rownames(MergedNodes_clean) <- MergedNodes_clean$Barcode
length(rownames(data) == rownames(MergedNodes_clean))

data_merge <- data.frame(matrix(1, nrow = length(table(MergedNodes_clean$Node)), ncol = ncol(data)))
rownames(data_merge) <- names(table(MergedNodes_clean$Node))
colnames(data_merge) <- colnames(data)

for (i in 1:ncol(data_merge)){
	data_merge[,i] <- tapply(data[,i], MergedNodes_clean$Node, mean)
}
data_merge <- data_merge[,rownames(col_anno_clean)]

fig_thres <- c(max((data_merge - 1)), min((data_merge - 1)))
fig_thres <- min(abs(fig_thres))
data_merge[data_merge > 1 + fig_thres] <- 1 + fig_thres
data_merge[data_merge < 1 - fig_thres] <- 1 - fig_thres

#Load the infercnv output of reference
data2 = read.delim("C:/MD_Anderson/4.PDAC_ST/6.SpatialInferCNV_by_patient_manual_anno/New_inferCNV_core_4/InferCNV_outputs_patient2_8/infercnv.references.txt", sep = " ", row.names = 1)
colnames(data2) <- gsub("\\.", "-", colnames(data2))
cnv_ref <- t(data2)

#Alternatively, use a clean reference without any CNV. This node should be clean, as the root node needs to be specified for subsequent evolutionary analysis.
cnv_ref_clean <- matrix(1, nrow = nrow(data_merge), ncol = ncol(data_merge))
rownames(cnv_ref_clean) <- rownames(data_merge)
colnames(cnv_ref_clean) <- colnames(data_merge)
data_merge_ref <- rbind(data_merge, "NT" = apply(cnv_ref_clean, 2, mean))

#Remove those week CNV signal, here we used 1.02 at cutoff
data_merge_ref_clean <- data_merge_ref
data_merge_ref_clean[abs(log(data_merge_ref_clean, 2)) < log(1.02, 2)] <- 1
data_dist <- dist(data_merge_ref_clean)

#NJ tree construction at subclone level
library(ggtree)
library(phangorn)
treeNJ2  <- NJ(data_dist)
my_tree <- as.phylo(treeNJ2)
allnode = my_tree$tip.label
cluster = allnode
names(cluster) = allnode

groupInfo <- lapply(levels(as.factor(cluster)), function(cc, cluster, allnode){
					return(names(cluster)[cluster == cc])
					},	cluster, allnode)					 
names(groupInfo) <- names(table(cluster))

#Specify the root clone as NT.
treeNJ2 <- (root(treeNJ2, which(treeNJ2$tip.label == "NT")))
MPT5 = groupOTU(treeNJ2, groupInfo)

#Calculate the spot number
dot_size <- 3 * (log(as.numeric(table(ST_case_2@meta.data$Clone)[grep("Clone|NT", names(table(ST_case_2@meta.data$Clone)))]))) - 3
dot_size[dot_size > 18] <- 18
dot_size <- dot_size + 8

#Replace the NT color
clone.color["NT"] <- "darkgrey"
gp5 <- ggtree(MPT5, aes(color = group), size = 2) + theme_tree() #+ geom_tippoint(alpha = 0.7)
gp5 <- gp5 + scale_color_manual(values=clone.color[grep("Clone|NT", names(clone.color))])
gp5 <- gp5 + geom_tippoint(size = dot_size)
gp5$data$label <- gsub("Clone", "", gp5$data$label)
gp5 <- gp5 + geom_tiplab(hjust = 0.55, vjust = 0.45, color = "black", size = 14) +  theme(legend.position="none")
gp5 <- flip(gp5, 1, 15)
gp5 <- flip(gp5, 24, 23)

#Add NT spot number
clone_composition_df_ref <- rbind(clone_composition_df, "NT" = 0)
temp <- table(ST_case_2@meta.data$Clone, ST_case_2@meta.data$site)[c("NT"),]
clone_composition_df_ref["NT", match(names(temp), colnames(clone_composition_df_ref))] <- temp

gp1 + geom_text(aes(label=node), hjust=0)					 
clone_composition_df_ref$node <- c(1:nrow(clone_composition_df_ref))
pie_df <- nodepie(clone_composition_df_ref, cols=1:8, alpha=0.9, color = site.color ) 

gp_combo5 <- inset(gp5, pie_df, height=0.12, width = 0.2,  hjust = -0.2, vjust = 0.05) + xlim(NA, max(gp5$data$x) + 0.2)
dis_node <- gp5$data[gp5$data$isTip,]
gp_combo5a <- inset(gp5, pie_df, height=0.12, width = 0.2, hjust = c(-0.25 - max(dis_node$x) + dis_node$x),  vjust = 0.05) +  xlim(NA, max(dis_node$x) + 0.25)

pdf("Figure 2_pt8_part1.pdf", 6, 8)
print(gp_combo5a)
dev.off()

#CNV heatmap at subclone level
d5 = fortify(MPT5)
dd = subset(d5, isTip)
spot_order <- dd$label[order(dd$y, decreasing=TRUE)]
spot_order <- as.character(gp5$data$group[gp5$data$isTip == TRUE])[order(gp5$data$y[gp5$data$isTip == TRUE], decreasing = T)]

row_anno_clean_merge <- clone_composition_df[,1:2]
colnames(row_anno_clean_merge)[1:2] <- "Clone"
row_anno_clean_merge[,1:2] <- rownames(row_anno_clean_merge)
row_anno_clean_merge <- as.data.frame(row_anno_clean_merge[,1])
colnames(row_anno_clean_merge) <- "Clone"
rownames(row_anno_clean_merge) <- row_anno_clean_merge[,1]
row_anno_clean_merge[nrow(row_anno_clean_merge)+1,] <- "NT"
rownames(row_anno_clean_merge)[nrow(row_anno_clean_merge)] <- "NT"
data_merge <- data_merge[,rownames(col_anno_clean)]

data_merge_map <- data_merge_ref
colnames(data_merge_map)[1:length(data_merge_map)] <- " "
colnames(data_merge_map)[mid_chr_mid] <- gsub("chr", "", names(seg.chr))
data_merge_map[abs(log(data_merge_map, 2)) < log(1.02, 2)] <- 1

aCol2 <- aCol
aCol2$Clone[names(aCol2$Clone) == "NT"] <- "darkgrey"
aCol2$Clone <- aCol2$Clone[grep("Clone|NT", names(aCol2$Clone))]

pheatmap(data_merge_map[spot_order,], scale = 'none', color = colorRampPalette(c("navy", "white", "firebrick"))(50), angle_col = "90", 
         cluster_cols = F, gaps_row = NULL, gaps_col = seg.chr, annotation_colors = aCol2, fontsize = 8,
		 annotation_row = row_anno_clean_merge,
         annotation_col = NA, 
         cluster_rows = F, show_colnames = T, show_rownames = F, border_color = NA,
		 width = 6,
         height = 3.5,
		 cutree_col = 2, 
		 file = paste0("Figure 2_pt8_part2.png") 	)

# CNV score 
data1_clean <- data1
data1_clean[abs(data1_clean - 1) < 0.02] <- 1
ST_case_2_tumor_cnv <- apply(data1_clean, 2, sd)

ST_case_2_tumor <- subset(ST_case_2, idents = names(table(Idents(ST_case_2)))[ grep("Clone", names(table(Idents(ST_case_2))))  ])
ST_case_2_tumor_cnv <- ST_case_2_tumor_cnv[match(rownames(ST_case_2_tumor@meta.data), names(ST_case_2_tumor_cnv))]
length(names(ST_case_2_tumor_cnv) == rownames(ST_case_2_tumor@meta.data))
ST_case_2_tumor@meta.data$InferCNV_score <- ST_case_2_tumor_cnv

Idents(ST_case_2_tumor) <- factor(Idents(ST_case_2_tumor), levels = rev(spot_order))
p0 <- VlnPlot(ST_case_2_tumor, features = c("InferCNV_score"), pt.size=0) + theme(legend.position = "none") + xlab("")
p1 <- p0 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  scale_fill_manual(values = clone.color) + coord_flip() + ggtitle("")
p2 <- p1 + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) 
p3 <- p2 + theme(panel.border = element_blank(), axis.line=element_blank()) 
p4 <- p3 + theme(axis.title.y=element_blank(), axis.text.y=element_blank())
p5 <- p4 + theme(axis.ticks=element_blank())
		
pdf("Figure 2_pt8_part3.pdf", 1.5, 4)	 
print (p5)
dev.off()

clone.color["NT"] <- "white"

pdf("Figure 2_pt8_part4.pdf", 16, 5)
SpatialDimPlot(ST_case_2, label = FALSE, label.size = 3, alpha = 1, cols = clone.color, stroke = 0.01, pt.size.factor = 5000)
dev.off()


