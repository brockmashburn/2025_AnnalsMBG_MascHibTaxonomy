
## Packages
#method http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

library(factoextra)
library(FactoMineR)
#library(corrplot)
#library(RColorBrewer)
#library(plyr)
#library(dlookr)
#library(ggsci)
#library(scales)
#library(writexl)
#library(ggstatsplot)
library(mclust)
#
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(svglite)

################################################################################
### MORPHOMETRICS OF THE MASCARENE HIBISCUS
################################################################################

# establish working directory
setwd("/Users/brockmashburn/Documents/Work/Hibiscus/2_Mascarenes/1_NewSpecies/morphometrics/analysis")
# data for all my samples
data.all <- read.csv("MascData.csv")
data.all[c(1:72),c(1:20)]
data.all

################################################################################
## Set Theme
################################################################################

# make a theme: taken from Donoghue et al. 2022 - https://doi.org/10.1038/s41559-022-01823-x
# source code: https://github.com/eaton-lab/Oreinotinus-phylogeny/blob/main/Analyses/Morphology-ecomorphs-and-convergence/plot_igesge_mix.R
theme_clim <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
          # text = element_text(family = "Arial"),
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"),
          strip.text = element_text(size = 10, color = "black", face = "bold.italic"),
          strip.background = element_rect(color = "white", fill = "white", size = 1))
}

# load  function for geom_flat_violin
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

################################################################################
#### A: Univariate analysis
################################################################################

data.all$myID1
data.all$CurrentID

# add a column to df but change the id of the reunion fragilis sample
# to do a boxplot with fragilis vs. all boryanus s.l.
data.all$ID3 <- data.all$CurrentID
data.all[63,"ID3"] <- "boryanus"
data.all$ID3

## remove types from analysis and keep only three boryanus s.l. species
data.bory <- data.all[-c(14,25:44),]

#comparisons for labeling
comparisons1 <- list(c("boryanus", "fragilis"))

comparisons2 <- list(c("boryanus", "dargentii"), c("boryanus", "igneus"),
                       c("dargentii", "igneus")
                       )

### Figure 4B. Calyx Length
names(data.all)
# select only columns of interest and remove missing samples
calyxlobe.len <- na.omit(data.all[,c("Collection", "ID3", "calyx_lobe_L")])
calyxlobe.len$ID3 = factor(calyxlobe.len$ID3, levels=c("boryanus", "fragilis"))

# make the plot
plot4B <- 
  ggplot(data = calyxlobe.len, aes(ID3, calyx_lobe_L, fill = ID3)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = calyx_lobe_L, color = ID3), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Calyx Lobe Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 75)) +
  scale_fill_manual(values = c("lightgrey", "#009E73")) +
  scale_color_manual(values = c("lightgrey", "#009E73")) +
  stat_compare_means(comparisons = comparisons1, hide.ns = TRUE,
                     method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
plot4B

### Figure 4C. Corolla Length
names(data.all)
# select only columns of interest and remove missing samples
flower.len <- na.omit(data.all[,c("Collection", "ID3", "flower_L")])
flower.len$ID3 = factor(flower.len$ID3, levels=c("boryanus", "fragilis"))

# make the plot
plot4C <- 
  ggplot(data = flower.len, aes(ID3, flower_L, fill = ID3)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = flower_L, color = ID3), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Flower Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 75)) +
  scale_fill_manual(values = c("lightgrey", "#009E73")) +
  scale_color_manual(values = c("lightgrey", "#009E73")) +
  stat_compare_means(comparisons = comparisons1, hide.ns = TRUE,
                     method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
plot4C

### Figure4D. Pedicel Length
names(data.bory)
# select only columns of interest and remove missing samples
pedicel.len <- na.omit(data.bory[,c("Collection", "myID1", "pedicel_L")])

# make the plot
plot4D <- 
  ggplot(data = pedicel.len, aes(myID1, pedicel_L, fill = myID1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = pedicel_L, color = myID1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Pedicel Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 75)) +
  scale_fill_manual(values = c("#F0E442", "#0072B2", "#56B4E9")) +
  scale_color_manual(values = c("#F0E442", "#0072B2", "#56B4E9")) +
  stat_compare_means(comparisons = comparisons2, hide.ns = TRUE,
                     method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
plot4D

ggsave(plot4B, filename = "~/Documents/Work/Hibiscus/2_Mascarenes/1_NewSpecies/manuscript/figures/Fig4B.svg", 
       device = "svg", height = 6, width = 6, dpi = 600)
ggsave(plot4C, filename = "~/Documents/Work/Hibiscus/2_Mascarenes/1_NewSpecies/manuscript/figures/Fig4C.svg", 
       device = "svg", height = 6, width = 6, dpi = 600)
ggsave(plot4D, filename = "~/Documents/Work/Hibiscus/2_Mascarenes/1_NewSpecies/manuscript/figures/Fig4A.svg", 
       device = "svg", height = 6, width = 6, dpi = 600)

################################################################################
#### B: Multivariate analysis
################################################################################

#####
##### Model Based Clustering, all samples, 6 variables
names(data.all)
data.all$Number
data.6.norm <- scale(data.all[, -c(1:14)])
mc.6 <- Mclust(data.6.norm)
summary(mc.6)
mc.6$classification
mc.6$z

# save model based clustering results to a data frame
cluster.res <- data.frame(row.names=data.all$Number)
cluster.res$species <- data.all$myID1
cluster.res$RCC_mbc_k2 <- mc.all$classification

# Figure 4A
Fig4A <- fviz_mclust(mc.all, "classification",
                 show.clust.cent = FALSE, legend = "none", main = FALSE,
                 geom = c("point"), pointsize = 1,
                 repel= TRUE, shape = 20,
                 ggtheme = theme_clim(), palette = c("black","#009E73")) +
  theme(plot.subtitle = element_blank())
Fig4A

# save the plots
ggsave(Fig4A, filename = "~/Documents/Work/Hibiscus/2_Mascarenes/1_NewSpecies/manuscript/figures/HibMasc_clust1.tiff", 
       device = "tiff", height = 100, width = 100, units = "mm", dpi = 400)

### Figure S2A
S2A <- fviz_mclust(mc.all, "BIC", ggtheme = theme_clim())
S2A <- ggpubr::ggpar(S2A, title = "Model & Cluster Selection",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12)
S2A

## Figure S2B
S2B <- fviz_mclust(mc.all, "classification",
                 show.clust.cent = FALSE, legend = "none", main = FALSE,
                 geom = c("text", "point"), pointsize = 3,
                 repel= TRUE, shape = 20,
                 ggtheme = theme_clim(), palette = c("black","#009E73"))
S2B <- ggpubr::ggpar(S2B, title = "Clustering Results with Sample Numbers",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12,
                   legend = c(0.08, 0.14),
                   font.legend = 12,
                   ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank())
S2B

# perform a pca on the same dataset as the mb-clust
pca.6 <- prcomp(data.6.norm)
# show correlation between variables and dimension 1 (cos2 score)
pca.6$var$cos2

# Figure S2C 
S2C <- fviz_pca_var(pca.6, geom = c("arrow", "text"), repel = TRUE)
S2C <- ggpubr::ggpar(S2C, title = "Correlations of Variables to PC's",
                     font.main = c(18, "bold", "black"),      
                     font.x = 12, font.y = 12,
                     legend = "none",
                     ggtheme = theme_bw()) +
  theme(legend.title = element_blank())
S2C

# Figure S2D
S2D <- fviz_cos2(pca.6, choice = "var")
S2D <- ggpubr::ggpar(S2D, title = "Correlation Between Variables and PC1",
                     font.main = c(18, "bold", "black"),
                     font.x = 12, font.y = 12,
                     ggtheme = theme_bw()) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())
S2D

### put supplemental figs together into a figure
FigS2 <- ggarrange(S2A, S2B, S2C, S2D,
                         labels = c("A","B","C","D"), #labels given each panel
                         font.label = list(size = 24),
                         ncol = 2, nrow = 2, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
FigS2
# save the final plot
ggsave(FigS2, filename = "FigS2.pdf", height = 13, width = 13)

### save the clustering results as a csv
#write.csv(RCC.cluster.res, "RCC_cluster_results.csv")

#####
##### Model Based Clustering, bory only

names(data.all)
data.all$ID3
data.bory <- data.all[-c(25:44),]
data.bory$ID3
names(data.bory)
data.bory.norm <- scale(data.bory[,-c(1:14,21)])
mc.bory <- Mclust(data.bory.norm)
summary(mc.bory)
mc.NT$classification
# shows only one cluster
