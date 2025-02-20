##Figure 1
##Figure 1A and 1B was a plot for illustratoin purpose, it's generated with biorender and illustrator

##Figure 1C UMAP plot 
##------------------------------------------------------------------------
setwd("C:/MD_Anderson/4.PDAC_ST_github/PDAC_Figures")

library(Seurat)
library(RColorBrewer)
library(yarrr)

load("PDAC_combined.Rdata")

site.color <- c((brewer.pal(9, "Blues"))[8], brewer.pal(9, "Oranges")[8], brewer.pal(9, "RdPu")[6], brewer.pal(9, "YlGnBu")[5])
names(site.color) <- c("Pri", "LiM", "LuM", "PerM")

p1 <- DimPlot(PDAC_combined, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin") 

pdf("Figure 1C.pdf", 6, 5.5)
print (p1)
dev.off()

##Figure 1D tree plot 
##------------------------------------------------------------------------
library(ggplot2)
library(ggtree)
library(ggpubr)

load("Figure 1/Figure 1D_input.Rdata")

Figure1D <- ggarrange(pt1_tree + ggtitle("Pt-1"), 
					  pt3_tree + ggtitle("Pt-3"), 
					  pt6_tree + ggtitle("Pt-6"), 
					  pt7_tree + ggtitle("Pt-7"), 
					  pt8_tree + ggtitle("Pt-8"), 
					  pt9_tree + ggtitle("Pt-9"), 
					  pt5_tree + ggtitle("Pt-5"), 
					  pt11_tree + ggtitle("Pt-11"), 
					  pt2_tree + ggtitle("Pt-2"),
					  pt4_tree + ggtitle("Pt-4"), 
					  pt13_tree + ggtitle("Pt-13"), 
					  pt12_tree + ggtitle("Pt-12"), 
					  pt10_tree + ggtitle("Pt-10"), 
					  ncol = 13, nrow = 1)
							
pdf("Figure 1D.pdf", 20, 4)
print (Figure1D)
dev.off()

##Figure 1E UMAP plot For each patient's merged data, generate a DimPlot using:
##------------------------------------------------------------------------
load("PDAC_Pt.Rdata")

DimPlot(Pt_1, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_2, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_3, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_4, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_5, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_6, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_7, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_8, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_9, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_10, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_11, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_12, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")
DimPlot(Pt_13, cols = site.color, label = FALSE, group.by = "Site") + ggtitle("Tissue origin")

#=========================================================================
##Below is Figure 1 relevant code for the Extended Data Figure.

##Extended Data Figure 1A & B UMAP plot 
##------------------------------------------------------------------------
PDAC_combined_tumor <- subset(PDAC_combined, idents = c("Classical", "Intermediate", "Basal"))

lineage.color <- c(piratepal("southpark", trans = .1)[1:3], "grey")
names(lineage.color) <- c("Classical", "Intermediate", "Basal", "Control")

treat.color <- c("#8dd3c7", "#bc80bd")
names(treat.color) <- c("Untreated", "Treated")

PDAC_combined_tumor$Lineage <- factor(PDAC_combined_tumor$Lineage, levels = c("Classical", "Intermediate", "Basal"))
PDAC_combined_tumor$Pt_id <- factor(PDAC_combined_tumor$Pt_id, levels = paste0("Pt-", 1:13))
PDAC_combined_tumor$Site <- factor(PDAC_combined_tumor$Site, levels = c("Pri", "LiM", "LuM", "PerM"))

p0 <- DimPlot(PDAC_combined_tumor, cols = lineage.color, label = TRUE, group.by = "Lineage") #+  theme(legend.position="none")
p1 <- DimPlot(PDAC_combined_tumor, group.by = "Pt_id") + ggtitle("Patient")
p2 <- DimPlot(PDAC_combined_tumor, cols = treat.color, group.by = "Treated")	+ ggtitle("Treatment")

pdf("Extended Data Figure 1 A & B.pdf", 12, 5)
print (p1 + p2)
dev.off()

##Extended Data Figure 4A UMAP plot 
##------------------------------------------------------------------------
p1 <- DimPlot(PDAC_combined_tumor, cols = lineage.color, group.by = "Lineage")	+ ggtitle("Lineage state")
p2 <- DimPlot(PDAC_combined_tumor, group.by = "Pt_id")	+ ggtitle("Patient")
p3 <- DimPlot(PDAC_combined_tumor, cols = site.color, group.by = "Site")	+ ggtitle("Site")

pdf("Extended Data Figure 4A.pdf", 16, 5)
print (p1 + p2 + p3)
dev.off()

