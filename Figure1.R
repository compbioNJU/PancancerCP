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