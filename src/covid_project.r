suppressMessages({library(Seurat)
                  library(Matrix)
                  library(dplyr)
                  library(ggplot2)
                  library(data.table)
                 library(stringr)
                 library(tibble)})

setwd('/data/brian/courses/singlecell')

nCoV.integrated <- readRDS(file = 'nCoV_mine_integrated_celltype_final_marker.rds')

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

saveRDS(nCoV.integrated, file = "nCoV_mine_integrated_celltype.rds")

##change to integrated assay
DefaultAssay(nCoV.integrated) <- "integrated"
dpi = 300
png(file="qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2) 

###tsne and umap
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')

write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()
saveRDS(nCoV.integrated, file = "nCoV_mine_integrated_celltype_final.rds")

#Draw Heatmap
dpi <- 300
hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
tt1 = DoHeatmap(object = subset(nCoV.integrated, downsample = 500), features = top10$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)

markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18')
hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) -> top30
var.genes = c(nCoV.integrated@assays$RNA@var.features,top30$gene,markers)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
saveRDS(nCoV.integrated, file = "nCoV_mine_integrated_celltype_final_marker.rds")

dpi <- 300
png(file="tsne_markers.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap_markers.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

####marker expression
dpi = 300
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF',
            'FCGR3A','TREM2','KRT18','HBB')
#markers = c('HBB')
#markers = c('DCN') --> No DCN
png(file="violin_marker.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
print(VlnPlot(object = nCoV.integrated, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
dev.off()
png(file="umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()
for(marker in markers){
  png(file=paste("violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = nCoV.integrated, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()
  png(file=paste("umap_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = nCoV.integrated, features = marker,cols = c("lightgrey","#ff0000")))
  dev.off()
}

markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF',
            'FCGR3A','TREM2','KRT18','HBB')

pdf(file="marker_heatmap.pdf", width = 10, height = 8)
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6))
print(pp)