#Adapted from JhuangLab: "https://github.com/JhuangLab/CytoTree"
#
#Script reflects steps taken to create major figures and analyze pre gated 
#samples to find unique B cell populations and compare differences between 
#healthy control and individuals with recent onset T1D.
#
#Load packages
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(CytoML)
library(LSD)
library(CytoTree)
library(ggplot2)
library(ggthemes)
library(umap)

#Upload raw FCS files and merge. FCS path must be modified based on download directory
fcs.path <- "Desktop/RawFCSfiles/Recent onsets or Healthy controls"
fcs.files <- list.files(fcs.path, pattern = '.fcs$', full = TRUE)
fcs.data <- runExprsMerge(fcs.files, comp = F, transformMethod = "arcsinh", mergeMethod = "all")

#rename channels with marker names
recol <- c(`Pr141Di<141Pr_CD22>` = "CD22",
           `Nd143Di<143Nd_CD5>` = "CD5",
           `Nd144Di<144Nd_CD69>` = "CD69",
           `Nd145Di<145Nd_CD95>` = "CD95",
           `Nd146Di<146Nd_IgD>` = "IgD",
           `Sm147Di<147Sm_CD20>` = "CD20",
           `Nd148Di<148Nd_CD38>` = "CD38",
           `Sm149Di<149Sm_CD24>` = "CD24",
           `Nd150Di<150Nd_CD86>` = "CD86",
           `Sm152Di<152Sm_CD36>` = "CD36",
           `Sm154Di<154Sm_CD40>` = "CD40",
           `Gd155Di<155Gd_CD27>` = "CD27",
           `Gd156Di<156Gd_CXCR3>` = "CXCR3",
           `Gd158Di<158Gd_CD10>` = "CD10",
           `Gd160Di<160Gd_antiFITC>` = "Tetanus binding",
           `Dy161Di<161Dy_CD80>` = "CD80",
           `Dy162Di<162Dy_CD11c>` = "CD11c",
           `Dy163Di<163Dy_CD72>` = "CD72",
           `Dy164Di<164Dy_CD71>` = "CD71",
           `Ho165Di<165Ho_CXCR5_EQ>` = "CXCR5",
           `Er166Di<166Er_CD23>` = "CD23",
           `Er168Di<168Er_CD138>` = "CD138",
           `Tm169Di<169Tm_IgG>` = "IgG",
           `Er170Di<170Er_antiBiotin>` = "Insulin binding",
           `Yb171Di<171Yb_PDL1>` = "PDL1",
           `Yb172Di<172Yb_IgM>` = "IgM",
           `Yb173Di<173Yb_CD21>` = "CD21",
           `Lu175Di<175Lu_CCR4_EQ>` = "CCR4",
           `Yb176Di<176Yb_PD1>` = "PD1")
colnames(fcs.data)[match(names(recol), colnames(fcs.data))] = recol
fcs.exp <- fcs.data[, recol]

#Group markers needed for dimensional reduction and clustering
markers <- c("CD22", "CD5", "CD69", "CD95", "IgD", "CD20", "CD38", "CD24", "CD86", "CD36", "CD40",
             "CD27", "CXCR3", "CD10", "CD80", "CD11c", "CD72", "CD71", "CXCR5", "CD23", "CD138", 
             "IgG", "PDL1", "IgM", "CD21", "CCR4", "PD1")

#create Cyt object
cyt <- createCYT(raw.data = fcs.exp, markers = markers,
                 normalization.method = "none")
cyt

#Starting with SOM clustering, xdim multiplied by ydim = total number of clusters
set.seed(1)
cyt <- runCluster(cyt, cluster.method = "som", xdim = 6, ydim = 6)
cyt <- processingCluster(cyt, k = 15)

#Phenograph clustering can be adjusted by changing K value.
cyt <- runCluster(cyt, cluster.method = "phenograph", k = 15)

#Adjust UMAP settings
umap.defaults
custom.settings = umap.defaults
custom.settings$min_dist = 0.5
custom.settings
cyt <- runUMAP(cyt, umap.config = custom.settings)

#Plot UMAP with clusters or by marker by changing color.by
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id", 
       alpha = 1, main = "tSNE", category = "categorical", show.cluser.id = T)

#optionally can also run tSNE
cyt <-runTSNE(cyt, check_duplicates = F)

#MST plots can be generated using the following functions:
cyt <- buildTree(cyt, dim.type = "raw")
cyt <- buildTree(cyt, dim.type = "umap", dim.use = 1:2)
cyt <- buildTree(cyt, dim.type = "tsne", dim.use = 1:2)

#plot MST
cyt@meta.data$branch.id <- paste0("B", cyt@meta.data$branch.id)
plotTree(cyt, color.by = "branch.id", show.node.name = T, cex.size = 2)

#clusters can be identified by heatmap
plotClusterHeatmap(cyt)

#modif branch id to refelect identity of clusters. Example:
branch.id <- cyt@meta.data$branch.id
branch.id[cyt@meta.data$cluster.id %in% c(1,5,8)] = "Immature B cells"
cyt@meta.data$branch.id <- branch.id

#re-plot to reflect cluster identity
plotTree(cyt, color.by = "branch.id", show.node.name = F, cex.size = 2)
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "branch.id", 
       alpha = 1, main = "UMAP", category = "categorical", show.cluser.id = T)

