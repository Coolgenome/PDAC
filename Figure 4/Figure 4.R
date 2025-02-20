##Figure 4
##Figure 4A, B, D, F, H, I are for illustratoin purposend were generated using the SpatialDimPlot function.

##Figure 4C Violin plot 
##------------------------------------------------------------------------
setwd("C:/MD_Anderson/4.PDAC_ST_github/PDAC_Figures")

library(ggplot2)
library(ggpubr)

load("Figure 4/Figure 4C_input.Rdata")

compare_1 <- list( c("Classical tumor bed", "Basal tumor bed"), c("Basal juxta", "Basal peri"), c("Basal juxta", "Basal tumor bed"), c("Basal tumor bed", "Basal peri"), c("Basal tumor bed", "Distant normal") )#, c("Int tumor bed", "Basal tumor bed"))

p1 <- ggviolin(figure4_c_input, x = "Region", y = "Value", fill = "Region", palette = (c("#bababa", "#756bb1",  "#fed98e", "#f03b20", "#f03b20",   "#f03b20"    )),
				add = "mean_sd", error.plot = "crossbar", add.params = list(fill = "white")) + xlab("") + ylab("") + theme(legend.position = "none") + ggtitle("myCAF score")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))	
p1 <- p1 + stat_compare_means(comparisons = compare_1, method = "wilcox.test") 	
p1 			

pdf("Figure 4C.pdf", 2.6, 5)
print (p1)
dev.off()

##Figure 4E Violin plot 
##------------------------------------------------------------------------
load("Figure 4/Figure 4E_input.Rdata")

compare_2 <- list( c("Peri", "Juxta"), c("Juxta", "Tumor bed"), c("Tumor bed", "Peri"), c("Tumor bed", "Distant normal") )

p1 <- ggviolin(figure4_e_input, x = "Region", y = "Value", fill = "Region", palette = (c("#bababa", "#756bb1",  "#fed98e", "#f03b20", "#f03b20",   "#f03b20"    )),
				add = "mean_sd", error.plot = "crossbar", add.params = list(fill = "white")) + xlab("") + ylab("") + theme(legend.position = "none") + ggtitle("apCAF score")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))	
p1 <- p1 + stat_compare_means(comparisons = compare_2, method = "wilcox.test") 	
p1 			

pdf("Figure 4E.pdf", 2.6, 5)
print (p1)
dev.off()
	
##Figure 4G Violin plot 
##------------------------------------------------------------------------
load("Figure 4/Figure 4G_input.Rdata")

compare_3 <- list( c("Classical tumor bed", "Basal tumor bed"), c("Basal juxta", "Basal peri"), c("Basal juxta", "Basal tumor bed"), c("Basal tumor bed", "Basal peri"), c("Basal tumor bed", "Distant normal") )#, c("Int tumor bed", "Basal tumor bed"))

p1 <- ggviolin(figure4_g_input, x = "Region", y = "Value", fill = "Region", palette = (c("#bababa", "#756bb1",  "#fed98e", "#f03b20", "#f03b20",   "#f03b20"    )),
				add = "mean_sd", error.plot = "crossbar", add.params = list(fill = "white")) + xlab("") + ylab("") + theme(legend.position = "none") + ggtitle("TAMs-4 Angio score")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))	
p1 <- p1 + stat_compare_means(comparisons = compare_3, method = "wilcox.test") 	
p1 			

pdf("Figure 4G.pdf", 2.6, 5)
print (p1)
dev.off()

##Figure 4J Violin plot 
##------------------------------------------------------------------------
load("Figure 4/Figure 4J_input.Rdata")

compare_4 <- list(c("TGFB1 - Tumor_classical", "TGFB1 - Tumor_basal"), c("TGFB1 - myCAF", "TGFB1 - Tumor_basal"), c("TGFB1 - myCAF", "TGFB1 - Tumor_basal"))
			
p1 <- ggviolin(figure4_j_input, x = "Pair", y = "Value", fill = "Pair", palette = c("#e5c494", "#377eb8", "#e41a1c"),
				add = "mean_sd", error.plot = "crossbar", add.params = list(fill = "white")) + xlab("") + ylab("") + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))	
p1 <- p1 + stat_compare_means(comparisons = compare_4, method = "wilcox.test") 
p1
					
pdf("Figure 4J.pdf", 3, 4.5)
print (p1)
dev.off()	

##Figure 4K Violin plot 
##------------------------------------------------------------------------
load("Figure 4/Figure 4K_input.Rdata")

compare_5 <- list(c("CXCL12 - Tumor_classical", "CXCL12 - Tumor_basal"), c("CXCL12 - Plasma", "CXCL12 - Tumor_classical"), c("CXCL12 - Plasma", "CXCL12 - Tumor_basal") )
			
p1 <- ggviolin(figure4_k_input, x = "Pair", y = "Value", fill = "Pair", palette = c("#7fc97f", "#377eb8", "#e41a1c"),
				add = "mean_sd", error.plot = "crossbar", add.params = list(fill = "white")) + xlab("") + ylab("") + theme(legend.position = "none")
p1 <- p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))	
p1 <- p1 + stat_compare_means(comparisons = compare_5, method = "wilcox.test") 		
p1 		

pdf("Figure 4K.pdf", 3, 4.5)
print (p1)
dev.off()	