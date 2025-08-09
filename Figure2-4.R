library(CytoTRACE)
cell1 <- seuratObj@meta.data %>% filter(primary_cluster %in% c('Epithelial') & 
                                          group %in% c("PT", "HM") &
                                          patient %in% sampaired & 
                                          (!tissue %in% 'Liver')) 
cell2 <- seuratObj@meta.data %>% filter(primary_cluster %in% c('Epithelial') & 
                                          group %in% c("PN", "PT") &
                                          #patient %in% sampaired & 
                                          (tissue %in% 'Liver')) 
cells <- cell1 %>% droplevels()
# cells <- cells %>% filter(subtype %in% (table(cells$subtype)[table(cells$subtype) > 20] %>% rownames())) %>% droplevels()
cells <- cells %>% slice_sample(n = 1000, by = subtype)
obj <- subset(seuratObj, cells=cells %>% rownames())
obj@meta.data <- droplevels(obj@meta.data)
obj <- RunUMAP(obj, dims = 1:50, n.components = 3)
DimPlot(obj) 

cyto_obj <- CytoTRACE(obj@assays$RNA@data %>% as.matrix(), 
                     ncores = 8, subsamplesize = 2000)
saveRDS(cyto_obj, 'epi.CytoTRACE.Obj.rds')
obj$CytoTRACE <- cyto_obj$CytoTRACE

pheno_data <- setNames(as.character(obj$subtype), Cells(obj))
plotCytoTRACE(cyto_obj, phenotype = pheno_data, 
              emb = Embeddings(obj, reduction = 'umap'))


pdf(paste0(fig_dir, 'Fig4.CytoTRACE.pdf'), height=6, width=6)
p0 <- DimPlot(obj, group.by = 'subtype', pt.size = 0.5, cols = xcol)
p1 <- DimPlot(obj, group.by = 'group', pt.size = 0.5, cols = groupcol)
p2 <- FeaturePlot(obj, features = "CytoTRACE", pt.size = 0.5, order = T) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))
print(p0 + NoLegend() + NoAxes() + ggtitle(NULL))
print(p1 + NoLegend() + NoAxes() + ggtitle(NULL))
print(p2 + NoLegend() + NoAxes() + ggtitle(NULL))
print(as_ggplot(get_legend(p0)) + 
        as_ggplot(get_legend(p1)) + 
        as_ggplot(get_legend(p2)))

ord <- obj@meta.data %>% group_by(subtype) %>% summarise(m=median(CytoTRACE)) %>%
  arrange(m) %>% pull(subtype) %>% as.character() %>% rev()
p1 <- obj@meta.data %>% mutate(subtype = factor(subtype, levels = ord)) %>%
  ggplot(aes(x=CytoTRACE, y=subtype, color=subtype)) +
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = xcol) +
  theme_bw()
p2 <- obj@meta.data %>% mutate(subtype = factor(subtype, levels = ord)) %>%
  ggplot(aes(x=-CytoTRACE, y=subtype, color=subtype)) +
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = xcol) +
  theme_bw()
p3 <- obj@meta.data %>% mutate(group = factor(group, levels = names(groupcol))) %>%
  ggplot(aes(x=CytoTRACE, y=group, color=group)) +
  geom_boxplot(outlier.color = NA) + scale_color_manual(values = groupcol) +
  facet_wrap(vars(tissue), scales = 'free_y', ncol = 1) + 
  theme_bw()
print((p2 + NoLegend()) + 
        (p3 + NoLegend()))
print((p1 + NoLegend()) + 
        (p2 + NoLegend()))
print(as_ggplot(get_legend(p1)) + as_ggplot(get_legend(p2)))
dev.off()


pltd <- obj@meta.data %>% mutate(rank=rank(-CytoTRACE)) %>%
  mutate(cnv=cnvScore[rownames(.)])
pltd$index <- cut(pltd$rank, breaks = seq(0, nrow(pltd), by=100))

pdf(paste0(fig_dir, 'Fig4.CytoTRACE.bin.pdf'), height=8.27, width=8.27)
x <- table(pltd$index, pltd$subtype)
x <- x / rowSums(x)
rownames(x) <- NULL
# barplot(t(x), col=xcol[colnames(x)], border = NA, axes = F, main="subtype")
x[x > 0.45] <- 0.45
pheatmap(t(x)[rev(ord),], cluster_rows = T, cluster_cols = F, 
         cellheight = 10, border_color = NA,
         color = colorRampPalette(rev(brewer.pal(11,"Spectral")))(255))
pheatmap(t(x)[rev(ord),], cluster_rows = F, cluster_cols = F, 
         cellheight = 10, border_color = NA,
         color = colorRampPalette(rev(brewer.pal(11,"Spectral")))(255))

HM <- CNV <- NULL
for(i in unique(pltd$tissue) %>% sort()){
  px <- pltd %>% filter(tissue %in% i)
  o <- px %>% group_by(index) %>% 
    summarise(cnv=median(cnv, na.rm = T)) %>% 
    as.data.frame() %>% na.omit()
  rownames(o) <- o$index
  o <- o[levels(pltd$index), ] %>% dplyr::select(cnv)
  colnames(o) <- i
  rownames(o) <- levels(pltd$index)
  if(is.null(CNV)){
    CNV <- o
  }else{
    CNV <- cbind(CNV, o)
  }
  
  x <- table(px$index, px$group)
  # x <- t(t(x) / colSums(x))
  x <- x / rowSums(x)
  rownames(x) <- NULL
  o <- data.frame(x[,'HM'])
  colnames(o) <- i
  if(is.null(HM)){
    HM <- o
  }else{
    HM <- cbind(HM, o)
  }
}
HM <- t(HM) / apply(HM, 2, max, na.rm=T)
pheatmap(HM, cluster_rows = T, cluster_cols = F, 
         cellheight = 10, border_color = NA,
         color = colorRampPalette(rev(brewer.pal(11,"Spectral")))(255))
# pheatmap(t(CNV), cluster_rows = T, cluster_cols = F, 
#          cellheight = 10, border_color = NA,
#          color = colorRampPalette(rev(brewer.pal(11,"Spectral")))(255))
dev.off()



library(monocle)

data <- obj@assays$RNA@counts
pd <- new('AnnotatedDataFrame', data = obj@meta.data)
fData <- data.frame(gene_short_name=obj@misc$geneName[row.names(data)], 
                    gene_id=row.names(data), 
                    row.names=row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle CellDataSet class
cds <- monocle::newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      ## lowerDetectionLimit = 0.5,
                      expressionFamily = VGAM::negbinomial.size());
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) +
                     2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) -
                     2*sd(log10(pData(cds)$Total_mRNAs)))
cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
             pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)

topMarkers <- obj@misc$markerGenes %>% dplyr::filter(p_val_adj < 0.05 & pct.1>0.25 & pct.2<pct.1 & !grepl("^MT-", name))
ordering_genes <- unique(topMarkers$gene)
cds <- monocle::setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, 
                       ## auto_param_selection = F, 
                       method = 'DDRTree') 
cds <- orderCells(cds)

root_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    R_counts <- pData(cds) %>% group_by(State) %>% 
      summarise(cyto=median(CytoTRACE)) %>% arrange(-cyto)
    return(as.numeric(R_counts[1,'State']))
  } else {
    return(1)
  }
}
cds <- monocle::orderCells(cds, root_state = root_state(cds))

plot_cell_trajectory(cds, color_by = "CytoTRACE")
plot_cell_trajectory(cds, color_by = "Pseudotime")
monocle::plot_cell_trajectory(cds, cell_size=0.5, color_by="subtype") + 
  scale_color_manual(values=xcol)
monocle::plot_cell_trajectory(cds, cell_size=0.5, color_by="tissue") + 
  scale_color_manual(values=tissuecol)

cds$CytoTime <- rescale(cds$CytoTRACE, to = range(cds$Pseudotime))

obj$Pseudotime <- pData(cds)[Cells(obj), 'Pseudotime'] 
obj$State <- pData(cds)[Cells(obj), 'State']

Idents(obj) <- 'subtype'
mkx <- FindAllMarkers(object = obj)
obj@misc$allmarkers <- mkx
saveRDS(obj, 'epi.seuratObj.rds')
saveRDS(cds, 'epi.monocle2.rds')

pdf(paste0(fig_dir, 'Fig4.CytoTRACE_Pseudotime.pdf'), height=6, width=6)
p0 <- DimPlot(obj, group.by = 'subtype', pt.size = 0.5, cols = xcol)
px <- DimPlot(obj, group.by = 'tissue', pt.size = 0.5, cols = tissuecol)
p1 <- DimPlot(obj, group.by = 'State', pt.size = 0.5, cols = brewer.dark2(9))
p2 <- FeaturePlot(obj, features = "CytoTRACE", pt.size = 0.5) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")))
p3 <- FeaturePlot(obj, features = "Pseudotime", pt.size = 0.5) + 
  scale_colour_gradientn(colours = rev(brewer.brbg(10)))
p4 <- pData(cds) %>% ggplot(aes(x=CytoTRACE, y=Pseudotime)) + 
  geom_point(aes(color=CytoTRACE)) + 
  geom_smooth(color='red', method = lm) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) + 
  theme_pubr() + NoLegend()
p5 <- pData(cds) %>% ggplot(aes(x=CytoTRACE, y=Pseudotime)) + 
  geom_point(aes(color=Pseudotime)) + 
  geom_smooth(color='red', method = lm) + 
  scale_colour_gradientn(colours = rev(brewer.brbg(10))) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) + 
  theme_pubr() + NoLegend()
p6 <- obj@meta.data %>% mutate(rank=rank(-CytoTRACE)) %>%
  mutate(index=cut(rank, breaks = seq(0, nrow(.), by=20))) %>%
  group_by(index) %>% summarise(x=median(CytoTRACE), 
                                y=median(Pseudotime, na.rm = T)) %>%
  ggplot(aes(x=x, y=y)) + geom_point(aes(color=x)) + 
  geom_smooth(color='black', method = lm) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) + 
  theme_pubr() + NoLegend()
p7 <- obj@meta.data %>% mutate(rank=rank(-CytoTRACE)) %>%
  mutate(index=cut(rank, breaks = seq(0, nrow(.), by=20))) %>%
  group_by(index) %>% summarise(x=median(CytoTRACE), 
                                y=median(Pseudotime, na.rm = T)) %>%
  ggplot(aes(x=x, y=y)) + geom_point(aes(color=y)) + 
  geom_smooth(color='black', method = lm) + 
  scale_colour_gradientn(colours = rev(brewer.brbg(10))) +
  stat_regline_equation(aes(label=paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) + 
  theme_pubr() + NoLegend()
print(p6 + p7 + ggplot() + ggplot())
print(p4 + p5 + ggplot() + ggplot())
print(p0 + NoLegend() + NoAxes() + ggtitle(NULL))
print(px + NoLegend() + NoAxes() + ggtitle(NULL))
print(p1 + NoLegend() + NoAxes() + ggtitle(NULL))
print(p2 + NoLegend() + NoAxes() + ggtitle(NULL))
print(p3 + NoLegend() + NoAxes() + ggtitle(NULL))
print(as_ggplot(get_legend(p0)))
print(as_ggplot(get_legend(p1)) + 
        as_ggplot(get_legend(p2)) + 
        as_ggplot(get_legend(p3)))
dev.off()

diff_test_x <- monocle::differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(CytoTime)")
saveRDS(diff_test_x, 'epi.pseudotime_diff_test.rds')
sig_gene_x <- row.names(subset(diff_test_x, qval < 0.001))
sig_genes <- refGenes[intersect(sig_gene_x, DEGs)] 


glist <- list()
for(tx in unique(cell1$tissue)){
  cells <- cell1 %>% filter(tissue %in% tx) %>% droplevels()
  # cells <- cells %>% filter(subtype %in% (table(cells$subtype)[table(cells$subtype) > 20] %>% rownames())) %>% droplevels()
  cells <- cells %>% slice_sample(n = 1000, by = subtype)
  objx <- subset(seuratObj, cells=cells %>% rownames())
  objx@meta.data <- droplevels(objx@meta.data)
  
  cyto_objx <- CytoTRACE(objx@assays$RNA@data %>% as.matrix(), 
                         ncores = 8, subsamplesize = 2000)
  objx$CytoTRACE <- cyto_objx$CytoTRACE
  
  
  library(monocle)
  data <- objx@assays$RNA@counts
  pd <- new('AnnotatedDataFrame', data = objx@meta.data)
  fData <- data.frame(gene_short_name=objx@misc$geneName[row.names(data)], 
                      gene_id=row.names(data), 
                      row.names=row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle CellDataSet class
  cdsx <- monocle::newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  ## lowerDetectionLimit = 0.5,
                                  expressionFamily = VGAM::negbinomial.size());
  cdsx <- estimateSizeFactors(cdsx)
  cdsx <- estimateDispersions(cdsx)
  cdsx <- detectGenes(cdsx, min_expr = 0.1)
  
  pData(cdsx)$Total_mRNAs <- Matrix::colSums(exprs(cdsx))
  upper_bound <- 10^(mean(log10(pData(cdsx)$Total_mRNAs)) +
                       2*sd(log10(pData(cdsx)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(cdsx)$Total_mRNAs)) -
                       2*sd(log10(pData(cdsx)$Total_mRNAs)))
  cdsx <- cdsx[,pData(cdsx)$Total_mRNAs > lower_bound &
                 pData(cdsx)$Total_mRNAs < upper_bound]
  cdsx <- detectGenes(cdsx, min_expr = 0.1)
  
  cdsx <- monocle::setOrderingFilter(cdsx, ordering_genes)
  
  cdsx <- reduceDimension(cdsx, max_components = 2, 
                          ## auto_param_selection = F, 
                          method = 'DDRTree') 
  cdsx <- orderCells(cdsx)
  
  root_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
      R_counts <- pData(cds) %>% group_by(State) %>% 
        summarise(cyto=median(CytoTRACE)) %>% arrange(-cyto)
      return(as.numeric(R_counts[1,'State']))
    } else {
      return(1)
    }
  }
  cdsx <- monocle::orderCells(cdsx, root_state = root_state(cdsx))
  
  cdsx$CytoTime <- rescale(cdsx$CytoTRACE, to = range(cdsx$Pseudotime))
  
  objx$Pseudotime <- pData(cdsx)[Cells(objx), 'Pseudotime'] 
  objx$State <- pData(cdsx)[Cells(objx), 'State']
  
  diff_test_x <- monocle::differentialGeneTest(cdsx, fullModelFormulaStr = "~sm.ns(CytoTime)")
  sig_genes <- refGenes[intersect(row.names(subset(diff_test_x, qval < 0.001)), DEGs)] 
  glist[[tx]] <- sig_genes
  saveRDS(objx, sprintf('%s.seuratObj.rds', tx))
  saveRDS(cdsx, sprintf('%s.monocle2.rds', tx))
  saveRDS(diff_test_x, sprintf('%s.diff_test.rds', tx))
}


library(mgcv)
library(scales)
# generalized additive model for logistic regression
pst <- pData(cds) %>% dplyr::select(CytoTRACE) %>% arrange(CytoTRACE)
dat <- GetAssayData(seuratObj)[sig_genes, rownames(pst)] %>% data.frame()

newdat <- data.frame(p.time=seq(from=min(pst$CytoTRACE, na.rm=T), to=max(pst$CytoTRACE, na.rm=T), length.out=500))
fit <- lapply(seq(1:nrow(dat)),function(x){
  df <- data.frame(value=as.numeric(dat[x,]), p.time=pst$CytoTRACE)
  mod <- gam(value~s(p.time, bs="cr"), data=df)
  pred <- predict(mod, newdat, type="response")
  zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
  return(zscore)
})
names(fit) <- rownames(dat)
fit <- do.call(rbind, fit)
fit <- t(apply(fit, 1, function(x){rescale(x, c(-1,1))}))

# reformat output
row.o <- apply(fit, 1, which.max)
fit <- fit[order(row.o, decreasing=F),]
TFs <- read.delim("hg.TF.txt", row.names = 1)
epiGs <- apply(avgExp1, 1, function(x){colnames(avgExp1)[which.max(x)]}) 
epiGs <- epiGs[epiGs == 'Epithelial'] %>% names() 
tfg <- intersect(epiGs, rownames(fit)[rownames(fit) %in% 
                        c(driverGs,refGenes[rownames(TFs)])])
saveRDS(fit, 'epi.CytoTRACE_heatmap.rds')

## x <- glist %>% unlist() %>% table()
## pltd <- fit[rownames(fit) %in% intersect(driverGs, names(x[x>3])), ]
pltd <- fit[rownames(fit) %in% c(
  refGenes[rownames(cnvpltd)], refGenes[rownames(TFs)]
), ]
tfg <- intersect(rownames(pltd), mksx)

# plot
mapal2 <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
pdf(paste0(fig_dir, 'Fig4.CytoTRACE_heatmap.pdf'), height=8.27, width=8.27)
ComplexHeatmap::pheatmap(pltd, color = mapal2, cluster_cols = F, 
                         fontsize_row = 3, fontsize = 8,
                         border_color = NA,
                         cluster_rows = F, show_colnames = F, 
                         main = sprintf("Genes = %s", nrow(pltd)))
phtm <- ComplexHeatmap::pheatmap(pltd, color = mapal2, cluster_cols = F, border_color = NA,
                                 fontsize_row = 1, fontsize = 8, show_rownames = F,
                                 cluster_rows = F, show_colnames = F, use_raster=T,
                                 main = sprintf("Genes = %s", nrow(pltd)))
ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=match(tfg, rownames(pltd)), 
                                                                  labels=tfg))
print(phtm+ha)
dev.off()


##### TODO *******
grep('CCL2|CTHRC1|YBX3', rownames(pltd))

gx <- rownames(pltd)
cytoglist <- list(
  C1=gx[1:145],
  C2=gx[146:185],
  C3=gx[186:293],
  C4=gx[294:length(gx)]
)
gotable2 <- NULL
for(k in names(cytoglist)){
  ego <- enrichGO(gene          = cytoglist[[k]],
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  gotable2 <- rbind(gotable2, ego@result %>% mutate(cluster=k))
}

pdf(paste0(fig_dir, 'Fig4.CytoTRACE_heatmap.GO.pdf'), height=8.27, width=7.8)
print(plotGO(gotable2))
dev.off()



gs <- setNames(rep(0, nrow(fit)), rownames(fit))
x <- glist %>% unlist %>% table()
x <- x[rownames(x) %in% names(gs)]
gs[names(x)] <- x



library(Nebulosa)
plot_density(obj, "ENSG00000072364", reduction = 'umap')

pheatmap(fit[intersect(rownames(fit), driverGs), ], 
         cluster_cols = F, cluster_rows = T, cutree_rows = 3)
