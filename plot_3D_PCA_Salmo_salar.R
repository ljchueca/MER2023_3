########################################################################
######                                                            ######
######                      Plot a 3D PCA                         ######
######                                                            ######
########################################################################


library(ggplot2)
library(plot3D)

# set working directory
setwd("/Users/ljchueca/Documents/_Trabajos_departamento/2023_MER_Master_Course/Salmo_salar_PCA")

# read bamlist used with ANGSD
bams <- read.table("bamlist")[,1]

# extract sample names from bamlist
samples1 <- sub("/.*/", "", bams)
samples <- sub("_removed_duplicates.bam", "", samples1)

# read covariance matrix generated with PCAngsd
salsal_cov <- as.matrix(read.table("snps.ld_pruned.hwe_filter.cov"))

# append sample names as row and column names to covariance matrix
dimnames(salsal_cov) <- list(samples, samples)

# perform PCA
pca <- prcomp(salsal_cov, scale = TRUE)

# scree plot
eigenval <- pca$sdev^2
explained_var <- 100*(eigenval/sum(eigenval))
df1 <- data.frame(prin_comp = c(seq(1, length(eigenval))),
                  explained_var)
ggplot(df1, aes(prin_comp, explained_var)) +
  geom_col(fill = c(rep("steelblue", 3),
                    rep("grey40", length(eigenval)-3))) +
  xlab("Principal components") +
  ylab("Explained variance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# save scree plot in PDF
ggsave("scree_plot.pdf", width = 170, height = 92, units = "mm")

# create data frame
df2 <- as.data.frame(pca$x)

# add column for populations
df2$pop <- ""
df2$pop[c(grep("C18_102|C18_103", row.names(df2)))] <- "North America"
df2$pop[c(grep("C19_004|C19_008|C19_010|C19_011|C19_012|C19_013|C19_014|C19_015", row.names(df2)))] <- "Russia"
df2$pop[c(grep("C18_108|C18_133", row.names(df2)))] <- "Iceland"
df2$pop[c(grep("C18_0|C18_106|C18_104|C18_105|C18_110", row.names(df2)))] <- "Baltic Sea" #The pattern C18_0 joined all this sample names C18_019|C18_079|C18_093|C18_095
df2$pop[c(grep("C18_107|Cuni_ENN|Cuni_WIN|C19_002|C19_15|C19_216", row.names(df2)))] <- "Europe"
# change the order of groups
df2$pop <- factor(df2$pop,
                 levels = c("Russia", "North America", "Iceland", "Baltic Sea","Europe"))

# add column for subpopulations
df2$subpop <- ""
df2$subpop[c(grep("C18_102|C18_103", row.names(df2)))] <- "North America"
df2$subpop[c(grep("C19_004|C19_008|C19_010|C19_011|C19_012|C19_013|C19_014|C19_015", row.names(df2)))] <- "Russia"
df2$subpop[c(grep("C18_108|C18_133", row.names(df2)))] <- "Iceland"
df2$subpop[c(grep("C18_0|C18_106|C18_104|C18_105|C18_110", row.names(df2)))] <- "Baltic Sea"
df2$subpop[c(grep("C18_107|Cuni_ENN", row.names(df2)))] <- "Scotland"
df2$subpop[c(grep("Cuni_WIN|C19_216", row.names(df2)))] <- "Ireland"
df2$subpop[c(grep("C19_002|C19_155", row.names(df2)))] <- "Norway"
df2$subpop[c(grep("C19_152|C19_153|C19_154", row.names(df2)))] <- "Sweden"
# change the order of groups
df2$subpop <- factor(df2$subpop,
                    levels = c("Russia", 
                               "North America",
                               "Iceland",
                               "Baltic Sea", 
                               "Scotland", "Ireland", "Norway", "Sweden"))

# set color blind friendly palette. Information about colours in ggplot: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("Russia" = "#F72585", #56B4E9
               "North America" = "#B5179E", #D55E00
               "Iceland" = "#7209B7",  #CC79A7
               "Baltic Sea" = "#03045E",  #009E73
               "Scotland" = "#0077B6",  #006046
               "Ireland" = "#4361EE", #0096C7 #F0E442
               "Norway" = "#4895EF", #48CAE4 #0072B2
               "Sweden" = "#4CC9F0")  #ADE8F4 #E69F00

# set point shapes
shapes <- c(18, 15, 17, 8, 16) # Here you can find information about ggplot point shapes: http://www.sthda.com/english/wiki/ggplot2-point-shapes

# save PCA plot in PDF
pdf("Salsal_3dpca_plot.pdf", width = 8, height = 6)

# 3D PCA plot
layout(matrix(c(1, 1, 1, 0,
                1, 1, 1, 0,
                1, 1, 1, 2,
                1, 1, 1, 2),
              nrow = 4, ncol = 4, byrow = TRUE))
par(mar = c(0, 0, 0, 0))
scatter3D(df2$PC1, df2$PC2, df2$PC3, bty = "g", theta = 30, phi = 45,
          colvar = NULL, colkey = FALSE, col = cbPalette[as.factor(df2$subpop)],
          pch = shapes[as.factor(df2$pop)], cex = 2.5, cex.lab = 2,
          xlab = paste("PC1 (", round(explained_var[1], 2), "%)", sep = ""),
          ylab = paste("PC2 (", round(explained_var[2], 2), "%)", sep = ""),
          zlab = paste("PC3 (", round(explained_var[3], 2), "%)", sep = ""))
plot.new()
legend("left",
       c("Russia", 
         "North America",
         "Iceland",
         "Baltic Sea", 
         "Scotland", "Ireland", "Norway", "Sweden"),
       col = cbPalette,
       # pch = c(18, 15, 17, 8, 16, 16, 16, 16), # Should be the same order than in line 88
       pch = c(18, 15, 17, 8, rep(16, 4)), # Instead writing '16' four times for the four European subpopulations we write 'rep(16,4)'
       pt.cex = 2.5, cex = 1.5, y.intersp = 1.5, bty = "n")

invisible(dev.off())

