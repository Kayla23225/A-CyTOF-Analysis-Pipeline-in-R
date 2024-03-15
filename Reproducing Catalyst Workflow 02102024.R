#Reproducing Catalyst Workflow for Practice#
#####
#Begin#

#Set working directory 
getwd()
setwd("C:/Users/Kayla/OneDrive/Documents/CyTOF Data Analysis")


#Download Packages
#If issue with HDCytoData package, use this: 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
#Then
BiocManager::install("HDCytoData")
#
install.packages(prepData)
BiocManager::install("CATALYST")

#Downloading the data
library(readxl)
url <- "https://zenodo.org/records/10039274/files"
md <- "PBMC8_metadata.xlsx"
download.file(file.path(url, md), destfile = md, mode = "wb")
md <- read_excel(md)
head(data.frame(md))

library(HDCytoData)
fs <- Bodenmiller_BCR_XL_flowSet()

##
#Set variables for each of your files
panel <- "PBMC8_panel_v3.xlsx"
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel)
head(data.frame(panel))

# spot check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs))


# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)




md <- "PBMC8_metadataMATCHED.xlsx"
panel <- "PBMC8_panel_v3.xlsx"

md <- read_xlsx(md)
panel <- read_xlsx(panel)

head(data.frame(panel))
head(data.frame(md))

#Reorder sample Ids
factor(md$sample_id, levels = unique(md$sample_id))
unique(md$sample_id)
factor(md$condition, levels = c("PBMC", "LCL"))
factor(md$condition, levels = unique(md$condition))


#Make a variable for the FCS files into a flowset
# This line didn't work (permission denied?) fs <- read.flowSet("C:/Users/Kayla/CyTOF Practice 20230117/FCSFiles") 

fcs_directory <- "C:/Users/Kayla/CyTOF Practice 20230117/FCSFilesMATCHED"
fs <- read.flowSet(dir(fcs_directory, pattern = "*\\.fcs", full.names = TRUE),
                   name.keyword = "FILENAME",
                   transformation = FALSE, truncate_max_range = FALSE)

#Check that all panel column names are in the flowset
all(panel$fcs_colname %in% colnames(fs))

# construct SingleCellExperiment
md$sample_id <- factor(md$sample_id,
                       levels = md$sample_id[order(md$condition)])

sce <- prepData(fs, panel, md, features = panel$fcs_colname)

#Ask how many cells are in your single cell experiment
n_cells(sce)

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

#Clustering!!!
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

#Make More rows within Cluster Expression
ClusterExpression <- plotClusterExprs(sce, k = "meta10", features = "type")

ClusterExpression$facet$params$nrow <- 4

ClusterExpression

??plotClusterExprs


plotDR(sce, "UMAP", color_by = "meta10")
plotDR(sce, "UMAP", color_by = "meta10", facet_by = "condition")
plotDR(sce, "UMAP", color_by = "CD1c", facet_by = "sample_id")

# facet by sample
plotDR(sce, "UMAP", color_by = "meta10", facet_by = "sample_id")
plotDR(sce, "UMAP", facet_by = "panel$antigen" color_by = "meta10")
plotDR(sce, "UMAP", color_by = c(), facet_by = "condition", ncol = 4)


#TSNE Plot
#NEed to download COWPLOT first
p1 <- plotDR(sce, "TSNE", color_by = "meta10") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta10")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))



####Reload####
library(readxl)
library(HDCytoData)
library(CATALYST)
library(SingleCellExperiment)
library(devtools)
library(flowCore)
library(cowplot)
library()



