#Reproducing the Robinson Pipeline Catalyst Workflow for Practice 03/15/2024#
#Link to Pipeline: It is important to refer to this pipeline when using your actual data
https://www.bioconductor.org/packages/release/workflows/vignettes/cytofWorkflow/inst/doc/cytofWorkflow.html 

#Set working directory 
getwd()
setwd("C:/Users/Kayla/OneDrive/Documents/Documents/GitHub/The-Robinson-Pipeline-for-Cytof-Analysis-in-R")

#Download and load Packages

install.packages("readxl")
install.packages("HDCytoData")
install.packages("CATALYST")
install.packages("SingleCellExperiment")
install.packages("devtools")
install.packages("flowCore")
install.packages("cowplot")

library(readxl)
library(HDCytoData)
library(CATALYST)
library(SingleCellExperiment)
library(devtools)
library(flowCore)
library(cowplot)
library()

#Downloading the data
url <- "https://zenodo.org/records/10039274/files"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))

#library(HDCytoData) is important in order to create your flowset
fs <- Bodenmiller_BCR_XL_flowSet()

#__________If there is an issue creating your flowset, you may have an issue with HDCytoData package. Use this: 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
#Then
BiocManager::install("HDCytoData")
#
install.packages(prepData)
BiocManager::install("CATALYST")
#___________________________________________________
#Set variables for each of your files
panel <- "PBMC8_panel_v3.xlsx"
#Read the excel file 
download.file(file.path(url, panel), destfile = panel, mode = "wb")
#Load your data frame you just downloaded
panel <- read_excel(panel)
head(data.frame(panel))

#spot check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs)) #If issue, there was most likely a problem loading your excel file possibly due to mismatch errors

#specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

#construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

#_______________If for some reason you get an error, try constructing SCE this way:
md$sample_id <- factor(md$sample_id,
                       levels = md$sample_id[order(md$condition)])

sce <- prepData(fs, panel, md, features = panel$fcs_colname)
#_______________________________________________________________
#Check that all panel column names are in the flowset
all(panel$fcs_colname %in% colnames(fs))

#Ask how many cells are in your single cell experiment..you should see all of your conditions
n_cells(sce)


#Here is where you plot all of your desired data...I only listed commonly used ones

#Histogram
plotCounts(sce, group_by = "sample_id", color_by = "condition")
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
p

#DiagnosticHeatmap
plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "RdYlBu")))

#MDS Plot
pbMDS(sce, color_by = "condition", label_by = "sample_id")

plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "RdYlBu")))

#NRS Score
plotNRS(sce, features = "type", color_by = "panel$antigen")

#Clustering...this is important to do before plotting heatmap
#Make sure you keep note of what seeds you use because they are unique to what you generate and important for replication
set.seed(2322)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 10, seed = 7435)

#Ploting the Heatmap
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta10",  scale = c("last"),
                bars = TRUE, perc = TRUE)


plotMultiHeatmap(sce,
                 hm1 = "type", hm2 = "pS6", k = "meta10",
                 row_anno = FALSE, bars = TRUE, perc = TRUE)

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(2322)
sce <- runDR(sce, "UMAP", cells = 5e3, features = "type")
sce <- runDR(sce, "TSNE", cells = 5e3, features = "type")

#__________________Make More rows within Cluster Expression
ClusterExpression <- plotClusterExprs(sce, k = "meta10", features = "type")

ClusterExpression$facet$params$nrow <- 4

ClusterExpression

??plotClusterExprs

plotDR(sce, "UMAP", color_by = "meta10")
plotDR(sce, "UMAP", color_by = "meta10", facet_by = "condition")
plotDR(sce, "UMAP", color_by = "CD1c", facet_by = "sample_id")
#_______________________________________________________________

# facet by sample
plotDR(sce, "UMAP", color_by = "meta10", facet_by = "sample_id")
plotDR(sce, "UMAP", facet_by = "panel$antigen" color_by = "meta10")
plotDR(sce, "UMAP", color_by = c(), facet_by = "condition", ncol = 4)

#TSNE Plot
#COWPLOT package is important here
p1 <- plotDR(sce, "TSNE", color_by = "meta10") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta10")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))


####MISC####

md <- "PBMC8_metadataMATCHED.xlsx"
panel <- "PBMC8_panel_v3.xlsx"

md <- read_xlsx(md)
panel <- read_xlsx(panel)

head(data.frame(panel))
head(data.frame(md))

#IF issues with creating your flowset, there is likely a probelm with recalling the files. A fix that worked for me was to download the FCS files directly onto my laptop/within my directory and call directly to that. 
# This line didn't work () fs <- read.flowSet("C:/Users/Kayla/CyTOF Practice 20230117/FCSFiles") 
#Try this fix:
fcs_directory <- "C:/Users/Kayla/CyTOF Practice 20230117/FCSFilesMATCHED"
fs <- read.flowSet(dir(fcs_directory, pattern = "*\\.fcs", full.names = TRUE),
                   name.keyword = "FILENAME",
                   transformation = FALSE, truncate_max_range = FALSE)

#Reorder sample Ids
factor(md$sample_id, levels = unique(md$sample_id))
unique(md$sample_id)
factor(md$condition, levels = c("PBMC", "LCL"))
factor(md$condition, levels = unique(md$condition))


