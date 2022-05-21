library(reticulate)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)
args <- commandArgs(trailingOnly = TRUE)
sc <- import("scanpy")

adata <- sc$read(args[[1]])

data <- t(adata$X)
rownames(data) <- rownames(adata$var)
# colnames(data) = rownames(adata$obs)

obs <- rownames(adata$obs)
obs_l = strsplit(obs,'_')
obs_name <- c()
i = 1
for (ob in obs_l)
{
    n = length(ob)
    batch = paste(ob[1:n-1],collapse= "")
    barcode = ob[n]
    obs_name[i] <- paste(batch,barcode,sep= "_")
    i = i+1
    # obs_name <- c(obs_name,new_name)
}
colnames(data) = obs_name

anno <- read.csv("~/code/batch/GSE158055_cell_annotation.csv")
cellName_l <- as.list(anno['cellName'])[[1]]
obs_l = strsplit(cellName_l,'_')
obs_name <- c()
i = 1
for (ob in obs_l)
{   
    n = length(ob)
    batch = paste(ob[1:n-1],collapse= "")
    barcode = ob[n]
    obs_name[i] <- paste(batch,barcode,sep= "_")
    i = i+1
    # obs_name <- c(obs_name,new_name)
}
rownames(anno) = obs_name
h_meta <- anno[colnames(data),]
colnames(h_meta) <- c('cellName','orig.ident','celltype','majorType')
h_meta <- h_meta[,c('orig.ident','celltype','majorType')]
merge <- CreateSeuratObject(data,meta = h_meta)

# QC
merge[["percent.mt"]] <- PercentageFeatureSet(merge,pattern = '^MT')
merge[["percent.rb"]] <- PercentageFeatureSet(merge,pattern = '^RP[SL]')
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(merge))
merge[["percent.HB"]]<-PercentageFeatureSet(merge, features=HB.genes)

minGene=500
maxGene=4000
maxUMI=15000
pctMT=10
pctHB=1

merge <- subset(merge, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
merge <- NormalizeData(merge) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()

saveRDS(merge,"~/code/batch/data/origin.rds")
