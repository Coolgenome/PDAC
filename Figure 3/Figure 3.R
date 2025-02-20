##Figure 3

##Figure 3A Bar plot 
##------------------------------------------------------------------------
setwd("C:/MD_Anderson/4.PDAC_ST_github/PDAC_Figures")

library(ggplot2)
library(ggpubr)

load("Figure 3/Figure 3A_input.Rdata")

p1 <- ggplot(ST_percentage_data, aes(x = Sample, y = count, group = Lineage, fill = Lineage)) + geom_bar(stat="identity") + scale_fill_manual(values = type.color)
p1 <- p1 + xlab("") + ylab("Count") + theme_bw() + ggtitle("")
p1 <- p1 + theme(legend.position = "none", axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) 
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p1 <- p1 + theme(legend.position = "right") + xlab("")

p2 <- ggplot(ST_percentage_data, aes(x = Sample, y = Percentage, group = Lineage, fill = Lineage)) + geom_bar(stat="identity") + scale_fill_manual(values = type.color)
p2 <- p2 + xlab("") + ylab("Percentage") + theme_bw() + ggtitle("")
p2 <- p2 + theme(legend.position = "none", axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) 
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2 <- p2 + theme(legend.position = "right") + xlab("")

pdf("Figure 3A.pdf", 12, 2.8)
print (p1)
print (p2)
dev.off()

##Extended Data Figure 4B Bar plot 
##------------------------------------------------------------------------
load("Figure 3/Extended Data Figure 4B_input.Rdata")

p3 <- ggplot(ST_percentage_data, aes(x = Sample, y = Percentage, group = Lineage, fill = Lineage)) + geom_bar(stat="identity") + scale_fill_manual(values = type.color)
p3 <- p3 + xlab("") + ylab("Percentage") + theme_bw() + ggtitle("")
p3 <- p3 + theme(legend.position = "none", axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black")) 
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p3 <- p3 + theme(legend.position = "right") + xlab("")
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = as.character(rep(site.color[c(1,2,4,3)], table(ST_percentage_data$Tissue)[c(4,1,3,2)]/3))))
p3 <- p3 + geom_vline(xintercept = as.numeric(cumsum((table(ST_percentage_data$Tissue)[c(4,1,3)])/3)) + 0.5, size = 2, color = "white")

pdf("Figure 3/Extended Data Figure 4B.pdf", 12, 4)
print (p3)
dev.off()


##Figure 3B Boxplot 
##------------------------------------------------------------------------
load("Figure 3/Figure 3B_input.Rdata")

T_violin_combo <- na.omit(T_violin_combo)

T_violin_combo_1 <- T_violin_combo[T_violin_combo$Lineage == "Classical", -2]
T_violin_combo_2 <- T_violin_combo[T_violin_combo$Lineage == "Intermediate", -2]
T_violin_combo_3 <- T_violin_combo[T_violin_combo$Lineage == "Basal", -2]

T_violin_combo_1$Tissue <- factor(T_violin_combo_1$Tissue, levels = c("Pri", "LiM", "PerM", "LuM"))
T_violin_combo_2$Tissue <- factor(T_violin_combo_2$Tissue, levels = c("Pri", "LiM", "PerM", "LuM"))
T_violin_combo_3$Tissue <- factor(T_violin_combo_3$Tissue, levels = c("Pri", "LiM", "PerM", "LuM"))

p1 <- ggpaired(T_violin_combo_1, x = "Tissue", y = "Percentage", id = "Pt", fill = "Tissue", line.color = "gray", line.size = 0.4, size = 2) +  scale_fill_manual(values = site.color) + xlab("Classical") + ylab("Percentage") + theme(legend.position = "none")
p2 <- ggpaired(T_violin_combo_2, x = "Tissue", y = "Percentage", id = "Pt", fill = "Tissue", line.color = "gray", line.size = 0.4, size = 2) +  scale_fill_manual(values = site.color) + xlab("Intermediate") + ylab("Percentage")+ theme(legend.position = "none")
p3 <- ggpaired(T_violin_combo_3, x = "Tissue", y = "Percentage", id = "Pt", fill = "Tissue", line.color = "gray", line.size = 0.4, size = 2) +  scale_fill_manual(values = site.color) + xlab("Basal") + ylab("Percentage")+ theme(legend.position = "none")

my_comparisons <- list( c("LuM", "PerM"), c("LuM", "LiM"), c("LuM", "Pri") )
p1 <- p1 + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater") )
p2 <- p2 + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater") )
p3 <- p3 + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "less") )

figure3B <- ggarrange(p1, p2, p3, ncol = 3)
figure3B

pdf("Figure 3/Figure 3B.pdf", 5, 2.8)
print (figure3B)
dev.off()

##Figure 3C D E Phylogenetic tree of tumor evolution
##------------------------------------------------------------------------
#The code for generating these figures can be found at the end of the section corresponding to the patient in Figure 2.
#Figure 2_pt5.R, Line 433
#Figure 2_pt9.R, Line 418
#Figure 2_pt13.R, Line 427

##Figure 3F Venny Plot 
##------------------------------------------------------------------------
library(ggvenn)
load("Figure 3/Figure 3F_input.Rdata")

x <- list("Pri Class" = rownames(ST_DEG_T1_1), "LiM Class" = rownames(ST_DEG_T2_1), "LuM Class" = rownames(ST_DEG_T3_1), "PerM Class" = rownames(ST_DEG_T4_1))
p1 <- ggvenn(x, fill_color = as.character(site.color), stroke_size = 0.5, set_name_size = 4)
x <- list("Pri Inter" = rownames(ST_DEG_T1_2), "LiM Inter" = rownames(ST_DEG_T2_2), "LuM Inter" = rownames(ST_DEG_T3_2), "PerM Inter" = rownames(ST_DEG_T4_2))
p2 <- ggvenn(x, fill_color = as.character(site.color), stroke_size = 0.5, set_name_size = 4)
x <- list("Pri Basal" = rownames(ST_DEG_T1_3), "LiM Basal" = rownames(ST_DEG_T2_3), "LuM Basal" = rownames(ST_DEG_T3_3), "PerM Basal" = rownames(ST_DEG_T4_3))
p3 <- ggvenn(x, fill_color = as.character(site.color), stroke_size = 0.5, set_name_size = 4)

pdf("Figure 3F.pdf", 6, 4.3)
print (p3)
dev.off()

pdf("Extended Data Figure 4 E & F.pdf", 6, 4.3)
print (p1)
print (p2)
dev.off()

##Figure 3G Heatmap Plot 
##------------------------------------------------------------------------
load("Figure 3/Figure 3G_input.Rdata")
library(pheatmap)

pdf("Figure 3G.pdf", 11.5, 3.5)
print (pheatmap( MP_lineage_cor, angle_col = 45, cluster_rows = F, cluster_col = F, display_numbers = MP_lineage_cor_text, fontsize = 15, color = colorRampPalette(c(type.color[1], "white", type.color[3]))(50)) )
dev.off()

##Figure 3H Heatmap Plot 
##------------------------------------------------------------------------
load("Figure 3/Figure 3h_input.Rdata")
library(pheatmap)

pdf("Figure 3H.pdf", 5, 1.7)
print (pheatmap( Gene_lineage_cor, angle_col = 45, cluster_rows = F, cluster_col = F, display_numbers = Gene_lineage_cor_text, color = colorRampPalette(c(type.color[1], "white", type.color[3]))(50)), fontsize = 25 )
dev.off()