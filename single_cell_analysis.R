# This script uses the package Seurat to analyze single sequencing data.
# The data for this analysis can be found here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472

library(Seurat)
library(tidyverse)

# This is a very large tsv. It will take a while to import (24658 rows, 49300 columns).

retina_full <- read.table("~/Documents/GIT/gene_analysis/GSE63472_P14Retina_merged_digital_expression.txt",
                          fill = T,
                          sep = "\t",
                          header = T,
                          row.names = 1
)

# Most Seurat functions take a specific format so convert df to SeuratObject.
# Impose cutoffs to reduce noise and make the analysis faster.
# Column (min.genes) needs to have at least 900 non-zero values.

all_retinas <- CreateSeuratObject(
  raw.data = retina_full,
  min.cells = 1,
  min.genes = 900,
  project = "all_retinas"
)

# Normalize the data
all_retinas <- NormalizeData(object = all_retinas,
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4)

# Identify variable genes. Variable genes indicate underlying cellular differences.
# Consistent genes are equivalently expressed across all cell types and cannot be
# used to differentiate cells
# y.cutoff emperically determined based on number of varible genes paper reported
all_retinas <- FindVariableGenes(object = all_retinas,
                                 mean.function = ExpMean,
                                 dispersion.function = LogVMR,
                                 # x.low.cutoff = 0.0125,
                                 # x.high.cutoff = 3,
                                 y.cutoff = 1.7)

# How many variable genes?
length(x = all_retinas@var.genes)

# Running out of memory so save gene names and
# remove original df once seurate object has been created
genes <- row.names(retina_full)
rm(retina_full)

# Use linear regression to identify and remove unwanted sources of variation
# Only using number of genes as variables to regress since no mito genes were aligned
all_retinas <- ScaleData(object = all_retinas, vars.to.regress = "nUMI")

# Use PCA to identify groups of genes that explain the greatest amounts of variation
# across the entire data set. Only use genes defined as variable previously
all_retinas <- RunPCA(object = all_retinas,
                      pc.genes = all_retinas@var.genes,
                      pcs.compute = 40)

# How genes within a single PC scores along its axis
# VizPCA(object = all_retinas, pcs.use = 1:2)
# Compare two PCs
# PCAPlot(object = all_retinas, dim.1 = 1, dim.2 = 2)

# Heatmap to visualize genes that are driving source of heterogenaity within a PC
# Number of PCs to display can be throttled, try changing pc.use to 1:n
# PCHeatmap(object = all_retinas,
#           pc.use = 1:12,
#           cells.use = 500,
#           do.balanced = TRUE,
#           label.columns = FALSE,
#           use.full = FALSE)

# A very important step is determining how many PCs to include. Some may be significant,
# others may be spurrious.
# JackStraw works to identify low p-values based on PCA score and iterates with
# random variations in data to identify PCs most affected by random variation.
# I don't understand this very well, go back and understand more.
# This is a heavy function and takes a while to run.
all_retinas <- JackStraw(object = all_retinas, num.pc = 40, num.replicate = 1000)

# Plot result of jackstraw analysis. Anything above the dotted line are significant PCs.
# Also pay attention to p-value (number next to PC number)
JackStrawPlot(object = all_retinas, PCs = 1:40)
# All 20 (max) appear to be significant

# Cluster the cells based on PCA score. Graph will not print.
all_retinas <- FindClusters(object = all_retinas,
                            reduction.type = "pca",
                            dims.use = 1:32,
                            resolution = 0.6,
                            print.output = 0,
                            save.SNN = F)

# Perform dimensional reduction using TSNE
all_retinas <- RunTSNE(object = all_retinas, dims.use = 1:32, do.fast = TRUE)

# End of main processing, save off so can start from this point in future
save(all_retinas, file = "~/Documents/GIT/gene_analysis/all_retinas.Rdata")

load("~/Documents/GIT/gene_analysis/all_retinas.Rdata")

# Plot graphical rendering of TSNE dimensional reduction.
# Individual clusters will be color coded and numbered.
TSNEPlot(object = all_retinas, do.label = TRUE, pt.size = 0.5)

# Don't know what cell type each cluster corresponds to.
# Identify gene identifiers (markers) of each cluster and take the top 10 to help identify cell types
all_retinas.markers <- FindAllMarkers(object = all_retinas, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
all_retinas_cluster_markers <- all_retinas.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
View(all_retinas_cluster_markers)

# Current list of cluster cell types, just from own knowledge or googling.
# 0 = Rods (RHO)
# 1 = Cones (OPN1SW/MW)
# 2 = ON Bipolars (PCP2, TRPM1)
# 3 = MG (RLBP)
# 4 = OFF Bipolar
# 5 = 
# 6 = ON Bipolar
# 7 =
# 8 =
# 9 =
# 11 =
# 12 =
# 13 =
# 14 =
# 15 =
# 16 =
# 17 =
# 18 =
# 19 =
# 20 =
# 21 =
# 22 =
# 23 =
# 24 =
# 25 =
# 26 =

# Finding markers takes a while so save off to load again in future
save(all_retinas.markers, file = "~/R_analysis/single_cell_RNAseq/all_retinas_markers")

# Find clusters where certain genes are expressed.
# These are known markers based on previous research.
rlbp <- grep("RLBP", genes, value = T) # Muller glia
lhx1 <- grep("LHX1", genes, value = T) # Horizontal cells
vsx2 <- grep("VSX", genes, value = T) # Bipolars (vsx1 and 2)
rho <- grep("RHO", genes, value = T) # Rhodopsin (rods? its widespread)
gabr <- grep("GABR", genes, value = T) # Gaba A receptors
opn1 <- grep("OPN1", genes, value = T) # Cones
atoh7 <- grep("ATOH7", genes, value = T)
pou <- grep("POU", genes, value = T)
map1 <- grep("MAP1", genes, value = T)
map2 <- grep("MAP2", genes, value = T)
pou4f <- grep("POU4F", genes, value = T)
prox1 <- grep("PROX1", genes, value = T)
bhlh <- grep("BHLH", genes, value = T)
nrn1 <- grep("NRN1", genes, value = T) # RGC?
gap43 <- grep("GAP43", genes, value = T) # RGC?
vsnl1 <- grep("VSNL1", genes, value = T) # RGC?
rgc <- grep("TUBB3", genes, value = T) # RGC
p27 <- grep("CDKN1B", genes, value = T)
ascl1 <- grep("ASCL1", genes, value = T)
ccn <- grep("CCND3", genes, value = T)
otx <- grep("OTX", genes, value = T)

FeaturePlot(object = all_retinas, features.plot = opn1, cols.use = c("yellow", "blue"), reduction.use = "tsne")

# In mice there are 2 cone types, but here they appear clustered into one (cluster 1)
# Isolate cells belonging to just that cluster and rerun analysis to see if can parse cell types

# In mice there are many different bipolar cell types. Here they are spread over 8 clusters
# (2, 4, 6, 9, 10, 13, 15, 23). Collect cells belonging to these clusters and rerun analysis to uncover known or novel cell subtypes

# Cones = 1
# BPs = 23, 10, 2, 9, 6, 4, 15, 13

cone_cells <- names(all_retinas@ident[all_retinas@ident == 1])

bp_cells <- names(all_retinas@ident[all_retinas@ident %in% c(23, 10, 2, 9, 6, 4, 15, 13)])

# If want to clean out stored objects run this
rm(all_retinas, all_retinas.markers, all_retinas_cluster_markers)

retina_full <- read.table("~/Documents/GIT/gene_analysis/GSE63472_P14Retina_merged_digital_expression.txt",
                          fill = T,
                          sep = "\t",
                          header = T,
                          row.names = 1
)

cones_only <- retina_full[,cone_cells] 

bp_only <- retina_full[,bp_cells]

all_cones <- CreateSeuratObject(
  raw.data = cones_only,
  min.cells = 1,
  min.genes = 900,
  project = "all_bps"
)

all_cones <- NormalizeData(object = all_cones,
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4)

all_cones <- FindVariableGenes(object = all_cones,
                             mean.function = ExpMean,
                             dispersion.function = LogVMR,
                             y.cutoff = 1.7)

all_cones <- ScaleData(object = all_cones, vars.to.regress = "nUMI")

all_cones <- RunPCA(object = all_cones,
                  pc.genes = all_cones@var.genes,
                  pcs.compute = 40)

all_cones <- JackStraw(object = all_cones, num.pc = 40, num.replicate = 1000)

JackStrawPlot(object = all_cones, PCs = 1:40)

all_cones <- FindClusters(object = all_cones,
                        reduction.type = "pca",
                        dims.use = 1:32,
                        resolution = 0.6,
                        print.output = 0,
                        save.SNN = F)

all_cones <- RunTSNE(object = all_cones, dims.use = 1:32, do.fast = TRUE)

TSNEPlot(object = all_cones, do.label = TRUE, pt.size = 0.5)

all_cones.markers <- FindAllMarkers(object = all_cones, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
all_cones_cluster_markers <- all_cones.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
View(all_cones_cluster_markers)

FeaturePlot(object = all_cones, features.plot = opn1, cols.use = c("yellow", "blue"), reduction.use = "tsne")



all_bps <- CreateSeuratObject(
  raw.data = bp_only,
  min.cells = 1,
  min.genes = 900,
  project = "all_bps"
)

all_bps <- NormalizeData(object = all_bps,
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4)

all_bps <- FindVariableGenes(object = all_bps,
                                 mean.function = ExpMean,
                                 dispersion.function = LogVMR,
                                 y.cutoff = 1.7)

all_bps <- ScaleData(object = all_bps, vars.to.regress = "nUMI")

all_bps <- RunPCA(object = all_bps,
                      pc.genes = all_bps@var.genes,
                      pcs.compute = 40)

all_bps <- JackStraw(object = all_bps, num.pc = 40, num.replicate = 1000)

JackStrawPlot(object = all_bps, PCs = 1:40)

all_bps <- FindClusters(object = all_bps,
                            reduction.type = "pca",
                            dims.use = 1:32,
                            resolution = 0.6,
                            print.output = 0,
                            save.SNN = F)

all_bps <- RunTSNE(object = all_bps, dims.use = 1:32, do.fast = TRUE)

TSNEPlot(object = all_bps, do.label = TRUE, pt.size = 0.5)

all_bps.markers <- FindAllMarkers(object = all_bps, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
all_bps_cluster_markers <- all_bps.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
View(all_bps_cluster_markers)

FeaturePlot(object = all_bps, features.plot = "PCP2", cols.use = c("yellow", "blue"), reduction.use = "tsne")



mg_cells <- names(all_retinas@ident[all_retinas@ident == 3])
mg_only <- retina_full[,mg_cells]

all_mg <- CreateSeuratObject(
  raw.data = mg_only,
  min.cells = 1,
  min.genes = 900,
  project = "all_bps"
)

all_mg <- NormalizeData(object = all_mg,
                           normalization.method = "LogNormalize",
                           scale.factor = 1e4)

all_mg <- FindVariableGenes(object = all_mg,
                               mean.function = ExpMean,
                               dispersion.function = LogVMR,
                               y.cutoff = 1.7)

all_mg <- ScaleData(object = all_mg, vars.to.regress = "nUMI")

all_mg <- RunPCA(object = all_mg,
                    pc.genes = all_mg@var.genes,
                    pcs.compute = 40)

all_mg <- JackStraw(object = all_mg, num.pc = 40, num.replicate = 1000)

JackStrawPlot(object = all_mg, PCs = 1:40)

all_mg <- FindClusters(object = all_mg,
                          reduction.type = "pca",
                          dims.use = 1:32,
                          resolution = 0.6,
                          print.output = 0,
                          save.SNN = F)

all_mg <- RunTSNE(object = all_mg, dims.use = 1:32, do.fast = TRUE)

TSNEPlot(object = all_mg, do.label = TRUE, pt.size = 0.5)

all_mg.markers <- FindAllMarkers(object = all_mg, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
all_mg_cluster_markers <- all_mg.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
View(all_mg_cluster_markers)