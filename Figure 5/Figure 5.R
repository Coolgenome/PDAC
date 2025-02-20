##Figure 5
##Figure 5A, B, C, F, G are for illustratoin purposend were generated using the SpatialDimPlot function.

##Figure 5D Violin plot 
##------------------------------------------------------------------------
setwd("C:/MD_Anderson/4.PDAC_ST_github/PDAC_Figures")

library(ggplot2)
library(ggpubr)
library(Seurat)

load("Figure 5/Figure 5D_input.Rdata")

p1 <- ggviolin(tumor_neighbour, x= "Cell_type", y = "Percentage", color = "Cell_type", fill = "Cell_type", palette = mycol_cosmx, add = "boxplot",
			   add.params = list(fill = "white"), size = 0.35) + xlab("") + 
				ggtitle("Tumour neighbour proportion") + theme(plot.title = element_text(hjust = 0.3))
p1 <- p1 + stat_compare_means(comparisons = list(c("CAF", "Plasma cell")), method = "wilcox.test", label.y = 80) + NoLegend()
p1 <- p1 + ylim(0, 100) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


p2 <- ggviolin(tumor_neighbour, x= "Cell_type", y = "Percentage", color = "Cell_type", fill = "Cell_type", palette = mycol_cosmx, add = "boxplot",
			   add.params = list(fill = "white"), size = 0.55) + xlab("") + 
				ggtitle("Plasma cell neighbour proportion") + theme(plot.title = element_text(hjust = 0.3))
p2 <- p2 + stat_compare_means(comparisons = list(c("CAF", "Tumour cell")), method = "wilcox.test", label.y = 80) + NoLegend() 
p2 <- p2 + ylim(0, 100) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Figure 5D.pdf", 7, 4)
print (p1 + p2)
dev.off()

##Figure 5E Heatmap plot 
##------------------------------------------------------------------------
load("Figure 5/Figure 5E_input.Rdata")

library(pheatmap)

p1 <- pheatmap(scale(t(CXCL12_matrix)), cluster_rows = F, cluster_cols = T, color = colorRampPalette(c("blue", "#fff7fb", "#ffff33"))(50),  angle_col = "90", 
		annotation_col = FOV_annotation[,2:1], annotation_colors = aCol )

p2 <- pheatmap(scale(t(CXCR4_matrix)), cluster_rows = F, cluster_cols = T, color = colorRampPalette(c("blue", "#fff7fb", "#ffff33"))(50),  angle_col = "90", 
		annotation_col = FOV_annotation[,2:1], annotation_colors = aCol )


pdf("Figure 5E_1.pdf", 18, 3)
print (p1)
dev.off()

pdf("Figure 5E_2.pdf", 18, 3)
print (p2)
dev.off()