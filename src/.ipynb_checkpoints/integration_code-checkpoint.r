suppressMessages({library(Seurat)
                  library(Matrix)
                  library(dplyr)
                  library(ggplot2)
                  library(data.table)
                 library(stringr)})
#library(Seurat)
#library(Matrix)
#library(dplyr)
#library(ggplot2)
#library(data.table)

samples <- as.list(list.files(path = 'covid_balf/singlecellstudy/data', pattern='*.h5', recursive= FALSE))
moderate <- samples[c(1,2,4)]
severe <- samples[c(3,5,6,10,11,12)]
healthy <- samples[seq(7,9,1)]

Get_paths<- function(listing){
    datadir<- list()
    for (i in 1:length(listing)){
        datadir[i]<- paste('/data/brian/courses/singlecell/covid_balf/singlecellstudy/data/',listing[i], sep="")}
    return(datadir)
}

Get_patient<- function(pat){
    patients <- list()
    for (i in pat){
        patients[i] <- strsplit(i[[1]], "[_]")[[1]][3]
    }
    return (patients)
}

Create_Seurat<- function(datas){
    datas<- Get_paths(datas)
    ndata<- data.table(ID = levels(as.factor(unlist(datas))))
    patients <- Get_patient(datas)
    for (i in datas){
            df.data <- data.frame()
            df <- data.frame()
            ndata$data <- sapply(ndata$ID, function(i, df.data, df){
                df.data <- Seurat::Read10X_h5(i ,use.names = TRUE, unique.features = TRUE)
                df <- CreateSeuratObject(counts = df.data , project = unlist(patients[i]), 
                                         min.cells = 3, min.features = 200)})
        for (j in 1:dim(ndata)[1]){
            ndata$data[[j]][['percent.mito']] <- PercentageFeatureSet(ndata$data[[j]], pattern = "^MT-")
        }
        return (ndata)
        }
}

moderate.df <- Create_Seurat(moderate)

severe.df<- Create_Seurat(severe)

healthy.df<- Create_Seurat(healthy) # I did not include GSM3660650

dpi <- 300
Preprocess_vis <- function(df){
    for (i in df[['data']]){
        print(dim(i))
        png(file=paste(as.character(unique(unlist(i@meta.data$orig.ident))),"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
        print(VlnPlot(object = i, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
        dev.off()
        png(file=paste(as.character(unique(unlist(i@meta.data$orig.ident))),"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
        print(FeatureScatter(object = i, feature1 = "nCount_RNA", feature2 = "percent.mito"))
        dev.off()
        png(paste(as.character(unique(unlist(i@meta.data$orig.ident))),"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
        print(FeatureScatter(object = i, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
        dev.off()
    }
}

Preprocess_vis(healthy.df)
Preprocess_vis(moderate.df)
Preprocess_vis(severe.df)

Preprocess_filter <- function(df){
    nFeature_RNA_low = 200
    nFeature_RNA_high = 6000
    nCount_RNA_req = 1000
    percent.mito.req = 10
    df_list <- list()
    for (i in df[['data']]){
        nCoV.seurat.filter <- subset(x = i, subset = nFeature_RNA > nFeature_RNA_low
                                     & nFeature_RNA < nFeature_RNA_high 
                                     & nCount_RNA > nCount_RNA_req & percent.mito < percent.mito.req)
        nCov.seurat.filter.norm <- NormalizeData(nCoV.seurat.filter, verbose = FALSE)
        nCov.seurat.filter.norm.selection <- FindVariableFeatures(nCov.seurat.filter.norm, selection.method = "vst", 
                                                                  nfeatures = 2000,verbose = FALSE)
        df_list <- append(df_list, nCov.seurat.filter.norm.selection)
    }
    return (df_list)}

healthy.df.filtered <- Preprocess_filter(healthy.df) # I did not include GSM3660650
moderate.df.filtered<- Preprocess_filter(moderate.df)
severe.df.filtered<- Preprocess_filter(severe.df)

metadata_change <- function(df, string){
    df.metadata <- df@meta.data 
    df.metadata$ID <- rownames(df.metadata)
    df.metadata <- df.metadata %>%
        mutate(ID = str_replace(ID, "-1", string))
    df <- AddMetaData(object = df,
                     metadata = df.metadata)
    rownames(df@meta.data) <- df@meta.data$ID 
    df@meta.data$ID <- NULL
    return (df)
}
healthy.df.filtered[[1]]<- metadata_change(healthy.df.filtered[[1]], 
                                            string = "_1")
healthy.df.filtered[[2]] <- metadata_change(healthy.df.filtered[[2]], 
                                            string = "_2")
healthy.df.filtered[[3]] <- metadata_change(healthy.df.filtered[[3]], 
                                            string = "_3")
moderate.df.filtered[[1]] <- metadata_change(moderate.df.filtered[[1]], 
                                            string = "_5")
moderate.df.filtered[[2]] <- metadata_change(moderate.df.filtered[[2]], 
                                            string = "_6")
moderate.df.filtered[[3]] <- metadata_change(moderate.df.filtered[[3]], 
                                            string = "_7")
severe.df.filtered[[1]] <- metadata_change(severe.df.filtered[[1]], 
                                            string = "_8")
severe.df.filtered[[2]] <- metadata_change(severe.df.filtered[[2]], 
                                            string = "_9")
severe.df.filtered[[3]] <- metadata_change(severe.df.filtered[[3]], 
                                            string = "_10")
severe.df.filtered[[4]] <- metadata_change(severe.df.filtered[[4]], 
                                            string = "_11")
severe.df.filtered[[5]] <- metadata_change(severe.df.filtered[[5]], 
                                            string = "_12")
severe.df.filtered[[6]] <- metadata_change(severe.df.filtered[[6]], 
                                            string = "_13")

all <- c(healthy.df.filtered, moderate.df.filtered, severe.df.filtered)

nCoV <- FindIntegrationAnchors(object.list = all, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))
saveRDS(nCoV.integrated, file = "nCoV_mine.rds")
