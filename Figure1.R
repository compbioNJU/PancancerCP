options(stringsAsFactors = F)
library(cowplot)
library(ggpubr)
library(ggforce)
library(pals)
library(scales)
library(ggthemes)
library(Seurat)
library(paletteer)
# library(plot1cell)
library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(xlsx)

library(pals)
library(ggthemes)
mycols <- unique(c(tableau_color_pal("Jewel Bright")(9),
                   tableau_color_pal("Winter")(8),
                   tableau_color_pal('Tableau 20')(20), 
                   tableau_color_pal('Classic 20')(20),
                   tableau_color_pal("Green-Orange-Teal")(12),
                   tableau_color_pal("Red-Blue-Brown")(12),
                   tableau_color_pal("Nuriel Stone")(9),
                   tableau_color_pal("Summer")(8),
                   tableau_color_pal("Color Blind")(10),
                   tableau_color_pal("Classic Color Blind")(10),
                   tableau_color_pal("Classic Green-Orange 12")(12),
                   tableau_color_pal("Classic Blue-Red 12")(12)
))


setwd("~/works/metastasis")
seuratObj <- readRDS('metastasis.final.rds')

# seuratObj@meta.data <- seuratObj@meta.data %>% mutate(
#   major_cluster=as.character(major_cluster),
#   primary_cluster=as.character(primary_cluster),
#   subtype=as.character(subtype)) %>% mutate(
#     major_cluster=ifelse(major_cluster=="Kupffer", "Macrophage", major_cluster),
#     primary_cluster=ifelse(primary_cluster=="Kupffer", "Macrophage", primary_cluster),
#     subtype=ifelse(subtype=="Kupffer_VSIG4", "Mac_VSIG4", subtype),
#     subtype=ifelse(subtype=="Kupffer_CCL18", "Mac_CCL18", subtype),
#     subtype=ifelse(subtype=="Epi_MMP7", "Endocrine_MMP7", subtype)
#   ) %>% droplevels()

refGenes <- seuratObj@misc$geneName

cellmeta <- seuratObj@meta.data %>% 
  select(primary_cluster, major_cluster, subtype) %>% 
  group_by(primary_cluster, major_cluster, subtype) %>% 
  summarise(n=n()) %>% data.frame(row.names = NULL) %>% 
  arrange(primary_cluster, major_cluster, -n) %>% 
  mutate(primary_cluster=factor(primary_cluster),
         major_cluster=factor(major_cluster, levels=unique(major_cluster)),
         subtype=factor(subtype, levels=unique(subtype)))
rownames(cellmeta) <- cellmeta$subtype

seuratObj@meta.data <- seuratObj@meta.data %>% 
  mutate(primary_cluster=factor(primary_cluster, levels=levels(cellmeta$primary_cluster))) %>%
  mutate(major_cluster=factor(major_cluster, levels=levels(cellmeta$major_cluster))) %>%
  mutate(subtype=factor(subtype, levels=levels(cellmeta$subtype))) 

mmks <- levels(cellmeta$subtype) %>% gsub("^[^_]*_","",x=.) %>% 
  strsplit("_") %>% unlist() %>% unique()

DEGs <- seuratObj@misc$markerGenes %>% 
  filter(avg_log2FC > log2(1.6)) %>% 
  pull(gene) %>% unique()
topmks <- refGenes[DEGs]

# output directories
data_dir <- "data/"
fig_dir <- 'figures/'
dir.create(data_dir, showWarnings = F)
dir.create(fig_dir, showWarnings = F)

# output directories
data_dir <- "data/"
fig_dir <- 'figures/'
dir.create(data_dir, showWarnings = F)
dir.create(fig_dir, showWarnings = F)

samples <- read.xlsx('finalmetadata.xlsx', sheetIndex = 1, header = T) %>%
  dplyr::rename(tissue=Tissue, study=Project, id=SampleID, sex=Sex, ref=Ref) %>%
  dplyr::select(tissue, study, id, sex, ref) %>% 
  mutate(group = gsub(".*-","",id)) %>%
  mutate(group = gsub("_[0-9]","",group)) %>% 
  mutate(patient = gsub("-.*","",id))
rownames(samples) <- samples$id

gcol <- setNames(palettes_d$awtools$spalette[c(-2,-5)], names(table(samples$group)))

seuratObj@meta.data <- seuratObj@meta.data %>% 
  mutate(sample=sampleID) %>% 
  mutate(samples[as.character(sample), c('study','patient','group','id','tissue','sex')])


table(cellmeta$primary_cluster)
majors <- table(cellmeta$primary_cluster) |> names()
subtypes <- table(cellmeta$major_cluster) |> names()
studies <- table(seuratObj$study) |> names()
tissues <- table(seuratObj$tissue) |> names()
groupx <- table(seuratObj$group) |> names()

groupcol <- setNames(palettes_d$ggthemes$Classic_10[1:length(groupx)], groupx)
studycol <- setNames(c(palettes_d$ggthemes$calc,
                       palettes_d$ggthemes$gdoc)[1:length(studies)], studies)
tissuecol <- setNames(c(palettes_d$ggthemes$stata_s2color,
                        palettes_d$ggthemes$colorblind)[1:length(tissues)], tissues)

mcol <- setNames(c(palettes_d$ggthemes$Tableau_20,
                   palettes_d$ggthemes$calc,
                   palettes_d$ggthemes$gdoc,
                   palettes_d$ggthemes$stata_s2color,
                   palettes_d$ggthemes$colorblind)[1:length(majors)],
                 majors)
scol <- setNames(c(palettes_d$ggthemes$Green_Orange_Teal,
                   palettes_d$ggthemes$Red_Blue_Brown,
                   palettes_d$ggthemes$Color_Blind,
                   palettes_d$ggthemes$Miller_Stone
)[1:length(subtypes)], subtypes)

pcols <- c(palettes_d$yarrr$google, palettes_d$RColorBrewer$Dark2,
           palettes_d$yarrr$appletv, palettes_d$ggthemes$Nuriel_Stone,
           palettes_d$khroma$bright, palettes_d$yarrr$basel,palettes_d$suffrager$classic,
           palettes_d$rcartocolor$Bold, palettes_d$suffrager$hanwell, 
           palettes_d$miscpalettes$brightPastel, palettes_d$yarrr$up,
           palettes_d$ggthemes$Classic_Cyclic, palettes_d$yarrr$espresso, 
           palettes_d$yarrr$decision, palettes_d$yarrr$nemo, mycols) %>% 
  substr(1,7) %>% toupper() %>% unique()

xcol <- setNames(pcols[1:nlevels(cellmeta$subtype)], levels(cellmeta$subtype))

# zcol <- setNames(xcol[cellmeta %>% pull(subtype)], 
#                  cellmeta %>% pull(cluster))

DimPlot(seuratObj, label = F, label.box = F, raster = T, 
        group.by = "subtype", cols = xcol) + NoLegend()
DimPlot(seuratObj, label = F, label.box = F, raster = T, 
        group.by = "major_cluster", cols = scol) + NoLegend()
DimPlot(seuratObj, label = F, label.box = F, raster = T, 
        group.by = "primary_cluster", cols = mcol) + NoLegend()


if(!is.null(names(rownames(seuratObj@assays$RNA@data)))){
  rownames(seuratObj@assays$RNA@data) <- names(rownames(seuratObj@assays$RNA@data))
}

if(!any(is.na(refGenes[rownames(seuratObj@assays$RNA@data)]))){
  rownames(seuratObj@assays$RNA@data) <- refGenes[rownames(seuratObj@assays$RNA@data)]
}
if(!any(is.na(refGenes[rownames(seuratObj@assays$RNA@scale.data)]))){
  rownames(seuratObj@assays$RNA@scale.data) <- refGenes[rownames(seuratObj@assays$RNA@scale.data)]
}
if(!any(is.na(refGenes[rownames(seuratObj@assays$SCT@data)]))){
  rownames(seuratObj@assays$SCT@data) <- refGenes[rownames(seuratObj@assays$SCT@data)]
}
if(!any(is.na(refGenes[rownames(seuratObj@assays$SCT@scale.data)]))){
  rownames(seuratObj@assays$SCT@scale.data) <- refGenes[rownames(seuratObj@assays$SCT@scale.data)]
}

avgExp0 <- AverageExpression(seuratObj, group.by = 'subtype', assays = 'SCT')$SCT
scor <- cor(avgExp0[topmks,], method = 's')
avgExp1 <- AverageExpression(seuratObj, group.by = 'primary_cluster', assays = 'SCT')$SCT
mcor <- cor(avgExp1[topmks,], method = 's')

annotation_colors <- list(primary=mcol, major=scol, subtype=xcol)
annotation_row <- cellmeta %>% mutate(primary=primary_cluster, 
                                      major=major_cluster,
                                      subtype=subtype) %>%
  dplyr::select(primary, major, subtype)


library(pheatmap)
pdf(paste0(fig_dir, 'Fig1.cor.pdf'), height = 12, width = 12)
phtm <- pheatmap(mcor, cellwidth = 10, cellheight = 10, fontsize = 10)
pheatmap(scor, cellwidth = 5, cellheight = 5, fontsize = 5, 
         annotation_row = annotation_row, 
         annotation_col = annotation_row, 
         annotation_colors = annotation_colors,
         border_color=NA)
majors <- phtm$tree_row$labels[phtm$tree_row$order] %>% as.character()
subs <- c()
for(x in majors){
  xx <- cellmeta %>% filter(primary_cluster==x) %>% pull(subtype) %>% as.character()
  if(length(xx) > 1){
    phtm <- pheatmap(scor[xx,xx], cellwidth = 10, cellheight = 10, fontsize = 10, main = x)
    subs <- c(subs, phtm$tree_row$labels[phtm$tree_row$order] %>% as.character())
  }else{
    subs <- c(subs, xx)
  }
}
n <- length(subs)
# subs <- subs[c(1:(n-7),(n-5):n,n-6)] ##### TODO ********** ######
pheatmap(scor[subs, subs], cellwidth = 5, cellheight = 5, fontsize = 5, 
         cluster_rows = F, cluster_cols = F, border_color = NA, 
         annotation_row = annotation_row, 
         annotation_col = annotation_row, 
         annotation_colors = annotation_colors)

dev.off()


seuratObj$primary_cluster <- factor(seuratObj$primary_cluster, levels = majors)

cellmeta <- cellmeta[subs, ] %>% mutate(cid=sprintf("c%03d",1:nrow(.)))
icol <- setNames(xcol[cellmeta %>% pull(subtype) %>% as.character()], 
                 cellmeta %>% pull(cid))
seuratObj$cid <- cellmeta[seuratObj$subtype %>% as.character(),'cid']

icor <- scor
rownames(icor) <- cellmeta[rownames(scor),'cid']
colnames(icor) <- cellmeta[colnames(scor),'cid']


library(dendextend)
getdend <- function(x){
  xx <- cellmeta %>% filter(as.character(primary_cluster) %in% x) %>% 
    pull(cid) %>% as.character()
  phtm <- pheatmap(icor[xx,xx], silent = T)
  phtm$tree_row$labels[phtm$tree_row$order] %>% as.character()
  
  dend <- phtm$tree_row %>% as.dendrogram()
  temp_col <- icol[xx][order.dendrogram(dend)]
  temp_col <- factor(temp_col, unique(temp_col))
  dend %<>% highlight_branches_col %>%
    color_branches(clusters=as.numeric(temp_col), col=levels(temp_col)) %>%
    ## set("nodes_pch", 19) %>%
    set("labels_colors", as.character(temp_col))
  dend
}

getdend("Plasma") %>% plot()
nodePar <- list(lab.cex = 0.6, pch = 19, cex = 0.7, col = "blue")


dend <- merge(
  merge(
    merge(getdend('Neutrophil'), 
          merge(getdend('DC'), 
                merge(getdend('Macrophage'),getdend('Monocyte')))), 
    merge(merge(getdend(c('Mast','Plasma')),getdend('B')),
          merge(getdend('CD4T'), merge(getdend('CD8T'),getdend('NK'))))
  ),
  merge(
    merge(merge(getdend(c('Photoreceptor','Endothelial')), merge(getdend('Fibroblast'), getdend('Pericyte'))),
          merge(getdend('Melanocyte'), merge(getdend('Endocrine'), getdend('Epithelial')))),
    getdend(c('Acinar','Hepatocyte'))
  )
) 

nodes <- get_leaves_attr(dend, "label")

seuratObj$minor <- seuratObj$cid <- factor(seuratObj$cid, levels = nodes)

meta <- cellmeta %>% mutate(minor=factor(cid,levels=nodes)) %>% arrange(minor)
rownames(meta) <- meta$minor 

p1 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "minor", cols = icol) + ggtitle(NULL) + NoAxes()
p2 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "subtype", cols = xcol) + ggtitle(NULL) + NoAxes()
p3 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "primary_cluster", cols = mcol) + ggtitle(NULL) + NoAxes()
p4 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "major_cluster", cols = scol) + ggtitle(NULL) + NoAxes()
p5 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "group", cols = groupcol) + ggtitle(NULL) + NoAxes()
p6 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "id") + ggtitle(NULL) + NoAxes()
samplecol <- setNames(ggplot_build(p6)$data[[1]][,1] %>% unique(), 
                      unique(seuratObj$id)) 
p7 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "patient") + ggtitle(NULL) + NoAxes()
p8 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "tissue", cols = tissuecol) + ggtitle(NULL) + NoAxes()
p9 <- DimPlot(seuratObj, raster = T, # label = T, label.box = T, 
              group.by = "study", cols = studycol) + ggtitle(NULL) + NoAxes()

pdf(paste0(fig_dir, 'Fig1.UMAP.pdf'), height = 8.27, width = 8.27)
print(p1 + NoLegend())
print(p2 + NoLegend())
print(p3 + NoLegend())
print(p4 + NoLegend())
print(p5 + NoLegend())
print(p6 + NoLegend())
print(p7 + NoLegend())
print(p8 + NoLegend())
print(p9 + NoLegend())
print(as_ggplot(get_legend(p1)))
print(as_ggplot(get_legend(p2)))
print(as_ggplot(get_legend(p3)) + as_ggplot(get_legend(p4)))
print(as_ggplot(get_legend(p5)))
print(as_ggplot(get_legend(p6)))
print(as_ggplot(get_legend(p7)))
print(as_ggplot(get_legend(p8)))
print(as_ggplot(get_legend(p8)))
dev.off()
rm(p1,p2,p3,p4,p5,p6,p7,p8,p9)


# cols <- setNames(mycols[1:nlevels(seuratObj$subtype)], levels(seuratObj$subtype))

pdf(paste0(fig_dir, "Fig1.UMAP-subtype.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'subtype', raster = F,
             cols = xcol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-minor.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'minor', raster = F,
             cols = icol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-primary_cluster.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'primary_cluster', raster = F,
             cols = mcol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-major_cluster.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'major_cluster', raster = F,
             cols = scol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-group.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'group', raster = F,
             cols = groupcol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-id.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'id', raster = F,
             pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-patient.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'patient', raster = F,
             pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-tissue.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'tissue', raster = F,
             cols = tissuecol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
pdf(paste0(fig_dir, "Fig1.UMAP-study.pdf"), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'study', raster = F, 
             # split.by = 'study', ncol = 4,
             cols = studycol, pt.size = 0.1, order = T) + NoAxes() 
print(p + ggtitle(NULL) + NoLegend())
print(as_ggplot(get_legend(p)))
dev.off()
rm(p)


library(pheatmap)
annotation_col <- meta %>% select(primary_cluster, major_cluster, cid) %>% 
  rename(major=primary_cluster, subset=major_cluster)
ann_colors <- list(major=mcol, subset=scol, cid=icol)

pdf(paste0(fig_dir, "Fig1.dendrogram.pdf"), height=10, width=10)
ids <- rev(nodes)
pheatmap(icor[ids, ids], cluster_rows = F, cluster_cols = F, border_color = NA,
         annotation_col=annotation_col, annotation_row=annotation_col, 
         labels_row = meta[ids,'cid'],
         annotation_colors = ann_colors, 
         cellwidth = 5, cellheight = 5, fontsize = 5)
op <- par(mar=c(2,1,1,1), mfrow=c(1,4))
dend %>% set("leaves_pch",19) %>% set("leaves_col",icol[nodes]) %>% 
  set("branches_lwd",1.5) %>% set("labels",meta[nodes,'cid']) %>% 
  set("labels_colors", icol[nodes]) %>% plot(horiz=T, xpd = T)
par(op)
op <- par(mfrow=c(4,1))
print(pal.bands(icol[rev(nodes)]))
print(pal.bands(xcol[meta[rev(nodes),'subtype']]))
print(pal.bands(scol[meta[nodes,'major_cluster'] |> unique()]))
print(pal.bands(mcol[meta[nodes,'primary_cluster'] |> unique()]))
par(op)
dev.off()


top10mks <- seuratObj@misc$top20Marker %>% dplyr::group_by(cluster) %>% 
  dplyr::top_n(n=10, wt=avg_log2FC) %>% pull(gene) %>% unique()

clusterExp0 <- AverageExpression(seuratObj, group.by = 'cid', 
                                 slot = 'data', assays = 'SCT')$SCT
mksx <- unique(c(refGenes[top10mks], mmks))

markerExp <- t(apply(clusterExp0[mksx, ], 1, scale))
colnames(markerExp) <- as.character(colnames(clusterExp0))

cord <- nodes
geneOrder <- data.frame(cluster=factor(colnames(markerExp)[apply(markerExp, 1, which.max)], levels = cord))
geneOrder$gene <- unlist(lapply(levels(geneOrder$cluster), function(i){
  order(markerExp[which(as.character(geneOrder$cluster)==i),i], decreasing = T)
}))
markerExp <- markerExp[order(geneOrder$cluster, geneOrder$gene), ]

pltd <- markerExp[,cord]
pltd[pltd > 3] <- 3
pltd[pltd < -3] <- -3
rownames(pltd) <- NULL 

# mapal <- colorRampPalette(c(rep('white',5),brewer.reds(5)))(256)
mapal <- colorRampPalette(c(rev(brewer.blues(8))[-(1:3)],brewer.ylorrd(5)))(256)

phtm <- ComplexHeatmap::pheatmap(pltd, use_raster=F, border_color = NA, fontsize = 5, color = mapal, 
                                 labels_row = NULL, cellwidth = 4, cluster_cols = F, cluster_rows = F,
                                 annotation_colors=ann_colors, annotation_col=annotation_col)
ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=match(mmks, rownames(markerExp)), 
                                                                  labels=mmks, labels_gp=grid::gpar(fontsize = 8)))
mk1 <- phtm + ha

pdf(paste0(fig_dir, "Fig1.marker-heatmap.pdf"), height=15, width=12)
print(mk1)
dev.off()
rm(mk1,phtm,ha)

gord <- meta[cord,] %>% pull(subtype) %>% gsub(".*_","",x=.)
gordx <- setNames(refGenes[match(gord,refGenes)],cord)
gordx <- gordx[!duplicated(gordx)]
Idents(seuratObj) <- factor(seuratObj$cid, levels = cord)
p <- DotPlot(seuratObj, features = gordx %>% unique()) + 
  scale_color_gradientn(colours=rev(brewer.rdylbu(20)), guide = "colourbar") +
  theme_minimal_grid() + theme(axis.text.x=element_text(angle=90, hjust=1))

pdf(paste0(fig_dir, "Fig1.marker.pdf"), height=25, width=28)
print(p)
op <- par(mfrow=c(4,1))
print(pal.bands(setNames(icol[names(gordx)],gordx)))
print(pal.bands(setNames(icol[names(gordx)],refGenes[gordx])))
par(op)
rm(p)
dev.off()


cellstat0 <- table(seuratObj$minor)
cellstat1 <- table(seuratObj$minor, seuratObj$group)
cellstat2 <- table(seuratObj$minor, seuratObj$sample)
cellstat3 <- table(seuratObj$minor, seuratObj$patient)
cellstat4 <- table(seuratObj$minor, seuratObj$tissue)
cellstat1 <- t(t(cellstat1) / colSums(cellstat1)) 
cellstat2 <- t(t(cellstat2) / colSums(cellstat2)) 
cellstat3 <- t(t(cellstat3) / colSums(cellstat3)) 
cellstat4 <- t(t(cellstat4) / colSums(cellstat4)) 

library(philentropy)
library(foreach) 
JSscore <- function(mat){
  JS <- foreach (i=1:nrow(mat), .combine=c) %dopar% {
    x <- rbind(mat[i,]/sum(mat[i,]), rep(1/ncol(mat), ncol(mat)))
    suppressMessages(JSD(x, est.prob="empirical"))
  }
  setNames(JS, rownames(mat))
}

js1 <- JSscore(cellstat1)
js2 <- JSscore(cellstat2)
js3 <- JSscore(cellstat3)
js4 <- JSscore(cellstat4)

p0 <- cellstat0 %>% as.data.frame() %>% 
  mutate(cluster=factor(Var1, levels = nodes), label=cluster) %>% 
  ggplot(aes(fill=cluster, x=log10(Freq), y=label)) + 
  geom_bar(stat = 'identity') + xlab("Cell number (log10)") + 
  scale_fill_manual(values = icol) + theme_minimal_vgrid() + NoLegend()
p1 <- cellstat1 %>% as.data.frame() %>% 
  mutate(cluster=factor(Var1, levels = nodes), label=cluster) %>% 
  mutate(group=factor(Var2)) %>% 
  ggplot(aes(fill=group, x=Freq, y=label)) + 
  geom_bar(position="fill", stat="identity") + xlab("Cell proportion") + 
  scale_fill_manual(values = groupcol) + theme_minimal_vgrid() + NoLegend()
p2 <- cellstat2 %>% as.data.frame() %>% 
  mutate(cluster=factor(Var1, levels = nodes), label=cluster) %>% 
  mutate(group=factor(Var2)) %>% 
  ggplot(aes(fill=group, x=Freq, y=label)) + 
  geom_bar(position="fill", stat="identity")  + xlab("Cell proportion") + 
  theme_minimal_vgrid() + NoLegend()

p3 <- data.frame(cluster=names(js1), score=js1, row.names = names(js1)) %>% 
  mutate(cluster=factor(cluster, levels = nodes),label=cluster) %>% 
  ggplot(aes(fill=cluster, x=score, y=label)) + 
  geom_bar(stat = 'identity') + xlab("JSD/group") + 
  scale_fill_manual(values = icol) + theme_minimal_vgrid() + NoLegend()
p4 <- data.frame(cluster=names(js2), score=js2, row.names = names(js2)) %>% 
  mutate(cluster=factor(cluster, levels = nodes),label=cluster) %>% 
  ggplot(aes(fill=cluster, x=score, y=label)) + 
  geom_bar(stat = 'identity') + xlab("JSD/sample") + 
  scale_fill_manual(values = icol) + theme_minimal_vgrid() + NoLegend()
p5 <- data.frame(cluster=names(js3), score=js3, row.names = names(js3)) %>% 
  mutate(cluster=factor(cluster, levels = nodes),label=cluster) %>% 
  ggplot(aes(fill=cluster, x=score, y=label)) + 
  geom_bar(stat = 'identity') + xlab("JSD/patient") + 
  scale_fill_manual(values = icol) + theme_minimal_vgrid() + NoLegend()
p6 <- data.frame(cluster=names(js4), score=js4, row.names = names(js4)) %>% 
  mutate(cluster=factor(cluster, levels = nodes),label=cluster) %>% 
  ggplot(aes(fill=cluster, x=score, y=label)) + 
  geom_bar(stat = 'identity') + xlab("JSD/tissue") + 
  scale_fill_manual(values = icol) + theme_minimal_vgrid() + NoLegend()


pdf(paste0(fig_dir, "Fig1.barplot.pdf"), height=15, width=6)
print(p0)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()
rm(p0,p1,p2,p3,p4,p5,p6)
rm(js1,js2,js3,js4,cellstat0,cellstat1,cellstat2,cellstat3,cellstat4)


######### hdWGCNA 
seuratObj$tissueGroup <- sprintf("%s-%s",seuratObj$tissue,seuratObj$group)
for(x in unique(seuratObj$tissueGroup)){
  cells <- names(seuratObj$tissueGroup)[seuratObj$tissueGroup %in% x]
  obj <- subset(seuratObj, cells = cells)
  saveRDS(obj, sprintf("%s.seuratObj.rds", x))
}
rm(obj)



######### CellPhoneDB 
## data 
counts_matrix <- GetAssayData(seuratObj, assay = 'RNA', slot="data")
if(!is.null(names(rownames(counts_matrix)))){
  rownames(counts_matrix) <- rownames(counts_matrix) %>% names()
}
cell_mata <- data.frame(Cell=names(seuratObj$major_cluster), 
                        cell_type=as.character(seuratObj$major_cluster), 
                        row.names = names(seuratObj$major_cluster),
                        stringsAsFactors = F)
for(i in unique(seuratObj$sample)){
  cids <- Cells(seuratObj)[seuratObj$sample == i]
  out <- counts_matrix[,cids] %>% as.matrix()
  write.table(out, gzfile(sprintf("cellphonedb/%s_counts.txt.gz",i)), sep="\t", quote=F)
  rm(out)
  out <- cell_mata[cids, ]
  write.table(out, gzfile(sprintf("cellphonedb/%s_meta.txt.gz",i)), sep="\t", quote=F, row.names = F)
}
for(i in unique(seuratObj$tissueGroup)){
  cids <- Cells(seuratObj)[seuratObj$tissueGroup == i]
  out <- counts_matrix[,cids] %>% as.matrix()
  write.table(out, gzfile(sprintf("cellphonedb/%s_counts.txt.gz",i)), sep="\t", quote=F)
  rm(out)
  out <- cell_mata[cids, ]
  write.table(out, gzfile(sprintf("cellphonedb/%s_meta.txt.gz",i)), sep="\t", quote=F, row.names = F)
}
rm(counts_matrix)





#### NMF 
nmfZscore <- NULL 
for(x in list.files('nmf/res2/', pattern = '*.Zscore.txt', full.names = T)){
  s <- read.delim(x, check.names = F, header = T)
  if(is.null(nmfZscore)){
    nmfZscore <- s
  }else{
    nmfZscore <- nmfZscore %>% full_join(s, by = 'gene')
  }
}
rownames(nmfZscore) <- nmfZscore$gene 
nmfZscore <- nmfZscore %>% dplyr::select(-gene)

nmfZscore <- nmfZscore[rowSums(is.na(nmfZscore))/ncol(nmfZscore)<0.1, ]
nmfZscore_scaled <- apply(nmfZscore, 2, scale) %>% 
  magrittr::set_rownames(rownames(nmfZscore))

zscorp <- Hmisc::rcorr(nmfZscore_scaled, type = 'pearson')
zscors <- Hmisc::rcorr(nmfZscore_scaled, type = 'spearman')

# zscorp <- cor(nmfZscore_scaled, use ='pairwise.complete.obs', method = 'p')
# zscors <- cor(nmfZscore_scaled, use ='pairwise.complete.obs', method = 's')

library(factoextra)
set.seed(31)
# function to compute total within-cluster sum of squares
p1 <- fviz_nbclust(zscorp$r, kmeans, method = "wss", k.max = 30) + 
  theme_minimal() + ggtitle("the Elbow Method - p")
p2 <- fviz_nbclust(zscors$r, kmeans, method = "wss", k.max = 30) + 
  theme_minimal() + ggtitle("the Elbow Method - s")
p <- matrix(-log10(p.adjust(zscorp$P)), 
            nrow = nrow(zscorp$P), 
            ncol = ncol(zscorp$P),
            dimnames = list(rownames(zscorp$P),colnames(zscorp$P)))
p[is.infinite(p)] <- max(p[is.finite(p)])
p[is.na(p)] <- 0
p3 <- fviz_nbclust(p, kmeans, method = "wss", k.max = 30) + 
  theme_minimal() + ggtitle("the Elbow Method - p")
p <- matrix(-log10(p.adjust(zscors$P)), 
            nrow = nrow(zscors$P), 
            ncol = ncol(zscors$P),
            dimnames = list(rownames(zscors$P),colnames(zscors$P)))
p[is.infinite(p)] = max(p[is.finite(p)])
p[is.na(p)] <- 0
p4 <- fviz_nbclust(p, kmeans, method = "wss", k.max = 30) + 
  theme_minimal() + ggtitle("the Elbow Method - s")

pdf("Fig2.nmf_metaprogram_nbclust.pdf", height=5, width=8.27)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


k <- 16
pltd <- zscorp$r ## sqrt(abs(zscorp))*sign(zscorp)
pltd[zscorp$P > 0.01] <- 0
pltd[pltd > 0.5] <- 0.5

annotation_row <- data.frame(id=rownames(pltd), row.names=rownames(pltd)) %>% 
  mutate(sample=gsub("_[0-9]+$","",x=id)) %>% 
  mutate(tissue=samples[sample,'tissue']) %>%
  mutate(group=samples[sample,'group']) %>%
  dplyr::select(group,tissue,sample)
ann_colors <- list(sample=samplecol, tissue=tissuecol, group=gcol)

# pltd[pltd < -0.5] <- -0.5
# cols <- brewer.rdylbu(100) %>% rev()
# cols <- c(rep(cols[1],round(100*sum(pltd<0)/length(pltd))), cols)
# phtm <- pheatmap(pltd, cutree_rows = k, cutree_cols = k) ## color = cols, 

pltd[pltd < 0] <- 0
phtm <- pheatmap(pltd, cutree_rows = k, cutree_cols = k)
metaProgram <- cutree(phtm$tree_row, k = k)[phtm$tree_row$order]
metaProgram <- setNames(sprintf("MP%02d",metaProgram), names(metaProgram))
mpcol <- setNames(palettes_d$awtools$bpalette[c(1:k)], unique(metaProgram))



annotation_row <- annotation_row %>% mutate(meta=metaProgram[rownames(.)])
ann_colors[['meta']] <- mpcol

dev.off()

pdf("Fig2.nmf_metaprogram_heatmap.pdf", height=15, width=15)
ComplexHeatmap::pheatmap(pltd, cutree_rows = k, cutree_cols = k, 
                         annotation_row = annotation_row, 
                         annotation_col = annotation_row,
                         annotation_colors = ann_colors, 
                         cellwidth = 1, cellheight = 1, 
                         fontsize_row = 1, fontsize_col = 1)
pheatmap(pltd, cutree_rows = k, cutree_cols = k, 
                 annotation_row = annotation_row, 
                 annotation_col = annotation_row,
                 annotation_colors = ann_colors, 
                 cellwidth = 1, cellheight = 1, 
                 fontsize_row = 1, fontsize_col = 1)
dev.off()



cell2program <- NULL 
for(x in list.files('nmf/res2', pattern = '*.usage.norm.txt', full.names = T)){
  m <- read.delim(x, check.names = F, header = T)
  o <- data.frame(cell=rownames(m), program=colnames(m)[apply(m,1,which.max)],
             value=apply(m,1,function(x){x[which.max(x)]}))
  cell2program <- rbind(cell2program, o)
}
cell2program <- cell2program %>% mutate(metaP=metaProgram[program]) %>% 
  mutate(cluster=seuratObj@meta.data[rownames(.),'cid'] %>% droplevels())
mpscore <- matrix(0, ncol = length(mpcol), nrow = ncol(malObj), 
                  dimnames = list(Cells(malObj),names(mpcol)))
malObj@meta.data <- cbind(malObj@meta.data, mpscore[rownames(malObj@meta.data),])


metaP2cls <- table(cell2program$cluster, cell2program$metaP)
p1 <- metaP2cls %>% as.data.frame() %>% 
  dplyr::rename(cluster=Var1, metaP=Var2, cell=Freq) %>% 
  mutate(metaP = factor(metaP, levels = unique(metaProgram) %>% rev())) %>% 
  ggplot(aes(y=metaP, fill=cluster, x=cell)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = icol[levels(cell2program$cluster)]) +
  theme(legend.position = 'none')
metaP2cls <- (100 * metaP2cls / rowSums(metaP2cls))
metaP2cls <- metaP2cls %>% as.data.frame() %>% 
  dplyr::rename(cluster=Var1, metaP=Var2, fract=Freq) %>% 
  mutate(metaP = factor(metaP, levels = unique(metaProgram) %>% rev()))
p2 <- metaP2cls %>% ggplot(aes(y=metaP, fill=cluster, x=fract)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = icol[levels(cell2program$cluster)]) +
  theme(legend.position = 'none')

pdf("Fig2.metaP.cellstat.pdf", height=8.27, width=8.27)
print(p1)
print(p2)
print(pal.bands(icol[levels(cell2program$cluster)]))
dev.off()


library(Nebulosa)
pdf("Fig2.metaP.score.pdf", height=8, width=5)
for(m in names(mpcol)){
  x <- cell2program %>% dplyr::filter(metaP==m)
  malObj@meta.data[rownames(x),m] <- x$value
  p <- FeaturePlot(malObj, features = m, raster = F, 
                   pt.size = 1, order = T, 
                   cols = c('grey',mpcol[m])) + NoAxes()
  p2 <- plot_density(malObj, m) + NoAxes()
  print(p + p2)
}
dev.off()


malObj$metaP <- 'NA'
malObj$metaP[rownames(cell2program)] <- cell2program$metaP
malObj$metaP <- factor(malObj$metaP, levels = c('NA', names(mpcol)))

pdf("Fig2.malObj.metaP.pdf", height=8.27, width=8.27)
p <- DimPlot(malObj, group.by = 'metaP', 
             order = T, cols = c('NA'='grey',mpcol))
print(p + NoLegend() + NoAxes() + ggtitle(NULL))
print(as_ggplot(get_legend(p)))
dev.off()


metaP2cls <- table(cell2program$cluster, cell2program$metaP)
x <- (100 * metaP2cls / rowSums(metaP2cls))
metaP2cls[x < 10] <- NA
# install.packages("ggalluvial")
library(ggalluvial)
p <- ggplot(data = data.frame(metaP2cls) %>% dplyr::filter(!is.na(Freq)),
       aes(axis2 = factor(Var2, levels=unique(metaProgram)), 
           axis1 = Var1, y = Freq)) +
  geom_alluvium(aes(fill = c(Var1))) +
  geom_stratum(aes(fill = Var1)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c(icol,mpcol)) +
  theme_void() 

pdf("Fig2.malObj.cell2metaP.pdf", height=8.27, width=5)
print(p + NoLegend() + NoAxes() + ggtitle(NULL))
print(as_ggplot(get_legend(p)))
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)
out <- NULL 
for(m in unique(metaProgram)){
  p <- names(metaProgram)[metaProgram == m]
  s <- nmfZscore_scaled[,p]
  gs <- apply(s, 1, median, na.rm=T) %>% sort(decreasing = T)
  ego <- gseGO(geneList     = gs,
               OrgDb        = org.Hs.eg.db,
               keyType      = 'SYMBOL', 
               ont          = "BP",
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.01,
               verbose      = FALSE)
  ego <- ego@result %>% mutate(meta=m)
  out <- rbind(out, ego)
}

toppath <- out %>% group_by(meta) %>% top_n(10, -p.adjust)
library(tidyr)
pltd <- out %>% dplyr::filter(ID %in% toppath$ID) %>%
  pivot_wider(id_cols = meta, values_from = p.adjust, 
                             names_from = Description) %>% as.data.frame()
pltd[is.na(pltd)] <- 1
rownames(pltd) <- pltd$meta 
pltd <- pltd %>% dplyr::select(-meta) %>% 
  log10() %>% abs() %>% sqrt() %>% t()

pdf("Fig2.nmf_metaprogram_pathway.pdf", height=12, width=8.27)
pheatmap(pltd[,unique(metaProgram)], cluster_cols = F, border_color = NA,
         cellheight = 1, cellwidth = 5, fontsize = 1)
pheatmap(pltd[,unique(metaProgram)], cluster_cols = T, border_color = NA,
         cellheight = 1, cellwidth = 5, fontsize = 1)
dev.off()



library(msigdbr)
msigdbr_show_species()
msigdbr_species()

# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

library(clusterProfiler)
out2 <- NULL 
for(m in unique(metaProgram)){
  p <- names(metaProgram)[metaProgram == m]
  s <- nmfZscore_scaled[,p]
  gs <- apply(s, 1, median, na.rm=T) %>% sort(decreasing = T)
  em <- GSEA(gs, TERM2GENE = m_t2g)
  em <- em@result %>% mutate(meta=m)
  out2 <- rbind(out2, em)
}

library(tidyr)
pltd <- out2 %>% pivot_wider(id_cols = meta, values_from = p.adjust, 
                            names_from = ID) %>% as.data.frame()
pltd[is.na(pltd)] <- 1
rownames(pltd) <- pltd$meta 
pltd <- pltd %>% dplyr::select(-meta) %>% 
  log10() %>% abs() %>% sqrt() %>% t()

pdf("Fig2.nmf_metaprogram_hallmark.pdf", height=12, width=8.27)
pheatmap(pltd[,unique(metaProgram)], cluster_cols = F)
pheatmap(pltd[,unique(metaProgram)], cluster_cols = T)
dev.off()




###Check and see the meta data info on your Seurat object
colnames(seuratObj@meta.data)  

###Prepare data for ploting
Idents(seuratObj) <- 'primary_cluster'
circ_data <- plot1cell::prepare_circlize_data(seuratObj, scale = 0.8)
set.seed(1234)
cluster_colors <- rand_color(length(levels(seuratObj)))
group_colors <- rand_color(length(names(table(seuratObj$group))))
rep_colors <- rand_color(length(names(table(seuratObj$primary_cluster))))

###plot and save figures
png(filename =  'circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "major_cluster",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()
