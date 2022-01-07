library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
#Read data
input_dir_PREV_1 <- "C:/Users/Spectre/Desktop/PREV_1_filtered_feature_bc_matrix"
PREV1.data <- Read10X(input_dir_PREV_1)
PREV_1<- CreateSeuratObject(counts = PREV1.data, project = "10X")
PREV_1@meta.data[,1]<-rep("PREV")

input_dir_PSTV_1 <- "C:/Users/Spectre/Desktop/PSTV_1_filtered_feature_bc_matrix"
PSTV1.data <- Read10X(input_dir_PSTV_1)
PSTV_1<- CreateSeuratObject(counts = PSTV1.data, project = "10X")
PSTV_1@meta.data[,1]<-rep("PSTV")

#Merge
AGGR <- merge(x = PREV_1, 
              y = PSTV_1, 
              add.cell.id = c("PREV_1", "PSTV_1"))

#Space release
rm(PREV_1)
rm(PREV1.data)
rm(PSTV_1)
rm(PSTV1.data)



#Get the gene list of mito. genes
AGGR[["percent.mt"]] <- PercentageFeatureSet(AGGR, pattern = "^MT-")
VlnPlot(AGGR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(AGGR, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AGGR, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

AGGR <- subset(AGGR, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#Normalization
AGGR <- NormalizeData(AGGR , normalization.method = "LogNormalize", scale.factor = 10000)
AGGR<- FindVariableFeatures(AGGR, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AGGR), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AGGR)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Dimensionality Reduction
#Scaling the data
all.genes <- rownames(AGGR)
AGGR <- ScaleData(AGGR, features = all.genes)
#Perform linear dimensional reduction
AGGR <- RunPCA(AGGR, features = VariableFeatures(object = AGGR))
#Examine and visualize PCA results a few different ways
print(AGGR[["pca"]], dims = 1:14, nfeatures = 5)
VizDimLoadings(AGGR, dims = 1:10, reduction = "pca")
DimPlot(AGGR, reduction = "pca")
DimHeatmap(AGGR, dims = 10, cells = 500, balanced = TRUE)

AGGR <- JackStraw(AGGR, num.replicate = 100)
AGGR <- ScoreJackStraw(AGGR, dims = 1:20)
JackStrawPlot(AGGR, dims = 1:20)
ElbowPlot(AGGR)

AGGR <- FindNeighbors(AGGR, dims = 1:9)
AGGR <- FindClusters(AGGR, resolution = 1)

AGGR <- RunTSNE(AGGR , dims = 1:13)
head(AGGR@reductions$tsne@cell.embeddings)
write.csv(AGGR@reductions$umap@cell.embeddings,file="UMAP.csv")
p2 <- DimPlot(AGGR, reduction = "tsne")
p2


#Find marker genes of clusters
cluster1.markers <- FindMarkers(AGGR, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
AGGR.markers <- FindAllMarkers(AGGR, only.pos = TRUE, min.pct = 0.25)
?FindMarkers
top5 <- AGGR.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top5,file="top5_markers_clusters.txt",sep="\t")
#DoHeatmap(AGGR, features = top5$gene)

#Cell annotation
library(SingleR)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(AGGR, slot="data")
clusters <- AGGR@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.table(celltype,file="AGGR_celltype.txt",sep="\t")

T_cell_CD4_Naive=c(0,3)
B_cell_Naive=c(10)
T_cell_CD8=c(4,11,14,16)
NK_cell=c(1,2,6,18)
Monocyte_CD14P=c(7,9)
Monocyte_CD16P=c(17)
#Monocyte_CD14P=c(3)
NK_cell_IL2=c(8)
B_cell_immature=c(12,13)
#B_cell_Memory=c(12)
T_cell_CD4_central_memory=c(5)
Platelets=c(15)

current.cluster.ids <- c(T_cell_CD4_Naive,
                         B_cell_Naive,
                         T_cell_CD8,
                         NK_cell,
                         Monocyte_CD14P,
                         Monocyte_CD16P,
                         NK_cell_IL2,
                         B_cell_immature,
                         T_cell_CD4_central_memory,
                         Platelets
)

new.cluster.ids <- c(rep("T_cell_CD4_Naive",length(T_cell_CD4_Naive)),
                     rep("B_cell_Naive",length(B_cell_Naive)),
                     rep("T_cell_CD8",length(T_cell_CD8)),
                     rep("NK_cell",length(NK_cell)),
                     rep("Monocyte_CD14P",length(Monocyte_CD14P)),
                     rep("Monocyte_CD16P",length(Monocyte_CD16P)),
                     rep("B_cell_immature",length(B_cell_immature)),
                     rep("T_cell_CD4_central_memory",length(T_cell_CD4_central_memory)),
                     rep("NK_cell_IL2",length(NK_cell_IL2)),
                     rep("Platelets",length(Platelets))
)

AGGR@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AGGR@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
plot1=row.names(AGGR@meta.data)
tsne_1<-DimPlot(AGGR, reduction = "tsne", group.by = "celltype", pt.size=0.5,cells = plot1)
tsne_2<-DimPlot(AGGR, reduction = "tsne", group.by = "orig.ident", pt.size=0.5,cells = plot1)

marker.genes<-row.names.data.frame(AGGR.markers)
marker.genes=bitr(marker.genes,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
cluster.go<-enrichGO(gene=marker.genes[,"ENTREZID"],keyType = "ENTREZID",OrgDb=org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = TRUE)

write.table(AGGR@active.ident,file="AGGR_CLUSTER_BARCODE.txt",sep="\t")
