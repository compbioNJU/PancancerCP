library(dplyr)
library(tidyverse)  
library(ggsankey)   
library(ggplot2)    
library(cols4all)   

IntOGen <- read.delim('IntOGen-DriverGenes.txt', row.names = 1)
COSMIC <- read.delim('Census_allTue Jan 23 04 12 40 2024.txt', row.names = 1)
driverGs <- union(rownames(IntOGen), rownames(COSMIC))

metadata <- readRDS('../metadata.rds')
load('../color.RData')

fs <- scan("infercnv.observations.list", what=character())
out <- NULL 
for(x in fs){
  o <- read.delim(x, sep = ' ', check.names = F)
  o <- data.frame(gene=rownames(o), o)
  if(is.null(out)){
    out <- o
  }else{
    out <- out %>% full_join(o, by = 'gene')
  }
  print(dim(out))
}
rownames(out) <- out$gene
out <- out %>% select(-gene)
out[is.na(out)] <- 1
out <- out[, colnames(out) %in% (cellmeta %>% filter(primary_cluster == 'Epithelial') %>% rownames())]
saveRDS(out, 'infercnv.observations.rds')

if(F){
  library(irlba)
  npcs <- 100
  pca.results <- irlba(A = as.matrix(x = t(out)), nv = npcs)
  feature.loadings <- pca.results$v
  
  sdev <- pca.results$d/sqrt(max(1, ncol(out) - 1))
  cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  
  rownames(x = feature.loadings) <- rownames(x = out)
  colnames(x = feature.loadings) <- paste0('PC_', 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = out)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  
  saveRDS(cell.embeddings, 'infercnv.cell.embeddings.rds')
  saveRDS(feature.loadings, 'infercnv.feature.loadings.rds')
  
  library(uwot)
  set.seed(seed = 2023)
  umap.output <- umap(X = cell.embeddings,
                      n_threads = 10, n_neighbors = 100,
                      n_components = 2, metric = 'cosine',
                      n_epochs = NULL, learning_rate = 1.0,
                      min_dist = 0.8, spread = 1.0,
                      set_op_mix_ratio = 1.0, local_connectivity = 1,
                      repulsion_strength = 1, negative_sample_rate = 5,
                      a = NULL, b = NULL, fast_sgd = FALSE,
                      verbose = TRUE, ret_model = FALSE
  )
  
  saveRDS(umap.output, 'infercnv.umap.rds')
}

library(Seurat)
cnvPca <- RunPCA(as.matrix(x = out), npcs = 50)
cnvUmap <- RunUMAP(cnvPca@cell.embeddings, seed.use = 2023, 
                   n.neighbors = 100, min.dist = 0.8)
saveRDS(cnvPca, 'infercnv.pca.rds')
saveRDS(cnvUmap, 'infercnv.umap.rds')

library(dplyr)
cnvNeighbor <- FindNeighbors(cnvPca@cell.embeddings)
cnvClusters <- FindClusters(cnvNeighbor$snn, resolution = 1)
cnvClusters$res.1 <- factor(cnvClusters$res.1, 
                            levels = table(cnvClusters$res.1) %>% 
                              sort(decreasing = T) %>% names())
saveRDS(cnvClusters, 'infercnv.cluster.rds')

cnvmeta <- cnvUmap@cell.embeddings %>% as.data.frame() %>% 
  mutate(cluster=cnvClusters[rownames(.),1]) %>%
  mutate(metadata[rownames(.),]) %>% droplevels()
saveRDS(cnvmeta, 'infercnv.cellmeta.rds')

SingleDimPlot(cnvmeta, dims = c('UMAP_1','UMAP_2'), 
              col.by = 'cluster', pt.size=0.25, 
              raster = F, shape.by = NULL,
              alpha.by = NULL, cols = pcols)

library(Seurat)
cnvObj <- CreateSeuratObject(counts = out, meta.data = cnvmeta)
saveRDS(cnvObj, 'infercnv.seuratObj.rds')
cnvObj <- FindVariableFeatures(object = cnvObj)
cnvObj <- ScaleData(object = cnvObj)
cnvObj <- RunPCA(object = cnvObj, npcs = 50)
cnvObj <- FindNeighbors(object = cnvObj)
cnvObj <- FindClusters(object = cnvObj, resolution = 1)
# cnvObj <- RunTSNE(object = cnvObj, dims = 1:30)
cnvObj <- RunUMAP(object = cnvObj, n.neighbors = 100, 
                     dims = 1:30, min.dist = 0.8)
x <- cnvUmap@cell.embeddings[rownames(cnvObj@reductions$umap@cell.embeddings), ]
colnames(x) <- colnames(cnvObj@reductions$umap@cell.embeddings)
cnvObj@reductions$umap@cell.embeddings <- x
DimPlot(cnvObj)


FeaturePlot(cnvObj, features = 'ENSG00000283627', 
            cols = c('grey', 'orangered1', 'red'))


load('../color.RData')
samplex <- metadata$sample %>% sort() %>% unique()
samplecol <- setNames(colorRampPalette(palettes_d$ggthemes$Classic_Cyclic)(length(samplex)), samplex)


pdf(paste0(fig_dir, 'Fig3.inferCNV.DimPlot.pdf'), height=8.27, width=8.27)
p1 <- DimPlot(cnvObj, group.by = 'cluster', cols = cnvcol, raster = F) + NoAxes() 
p2 <- DimPlot(cnvObj, group.by = 'subtype', cols = xcol[levels(cnvObj$subtype)], raster = F) + NoAxes() 
p3 <- DimPlot(cnvObj, group.by = 'sample', cols = samplecol, raster = F) + NoAxes() 
p4 <- DimPlot(cnvObj, group.by = 'group', cols = groupcol, raster = F) + NoAxes() 
p5 <- DimPlot(cnvObj, group.by = 'tissue', cols = tissuecol, raster = F) + NoAxes() 

print(DimPlot(cnvObj, group.by = 'cluster', label.box = T,
              cols = cnvcol, raster = F, label = T) + 
        NoAxes() + ggtitle(NULL) + NoLegend())
print(p1 + ggtitle(NULL) + NoLegend())
ggarrange(p2 + ggtitle(NULL) + NoLegend(),
          p3 + ggtitle(NULL) + NoLegend(),
          p4 + ggtitle(NULL) + NoLegend(),
          p5 + ggtitle(NULL) + NoLegend(),
          ncol = 2, nrow = 2)
print(as_ggplot(get_legend(p1)))
print(as_ggplot(get_legend(p2)))
print(as_ggplot(get_legend(p3)))
print(as_ggplot(get_legend(p4)))
print(as_ggplot(get_legend(p5)))
dev.off()

png(paste0(fig_dir, 'Fig3.inferCNV.DimPlot1.png'), height=800, width=800)
print(p1 + ggtitle(NULL) + NoLegend())
dev.off()
png(paste0(fig_dir, 'Fig3.inferCNV.DimPlot2.png'), height=800, width=800)
ggarrange(p2 + ggtitle(NULL) + NoLegend(),
          p3 + ggtitle(NULL) + NoLegend(),
          p4 + ggtitle(NULL) + NoLegend(),
          p5 + ggtitle(NULL) + NoLegend(),
          ncol = 2, nrow = 2)
dev.off()

rm(p1,p2,p3,p4,p5)
nodex <- c(levels(cnvmeta$subtype),
           levels(cnvmeta$cluster),
           unique(cnvmeta$sample), 
           unique(cnvmeta$group), 
           unique(cnvmeta$tissue))
cnvcol <- setNames(mycols[1:nlevels(cnvmeta$cluster)], 
                   levels(cnvmeta$cluster))
ncols <- c(xcol, cnvcol, samplecol,
           groupcol, tissuecol)
showlabs <- c(unique(cnvmeta$group), 
              unique(cnvmeta$tissue))

library(ggplot2)
library(ggsankey)
pdf(paste0(fig_dir, 'Fig3.infercnv.sankey.pdf'), width = 15, height = 8.27)
cnvmeta %>% dplyr::select(subtype, cluster, sample, group, tissue) %>%
  make_long(subtype, cluster, sample, group, tissue) %>% 
  mutate(label=ifelse(node %in% showlabs, node, NA)) %>%
  mutate(node = factor(node, levels=nodex)) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node, 
             label = label,
             fill = node)) +
  geom_sankey(
    flow.alpha = 0.9,
  ) +  
  geom_sankey_label(size = 3, color = 'white', fill = NA) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = ncols) + 
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Cancer cells")
dev.off()

library(entropy)
p1 <- apply(cnvmeta %>% ## filter(patient %in% sampaired) %>% 
        dplyr::select(cluster, sample) %>% 
        table(), 2, entropy, unit='log2') %>% 
  as.data.frame() %>% dplyr::rename(Entropy='.') %>% 
  dplyr::mutate(samples[rownames(.),c('group','tissue')]) %>% 
  mutate(group=factor(group, levels = c("NN","PN","PT","HM"))) %>%
  ggplot(aes(y=group, x=Entropy, color=group)) + 
  geom_boxplot() + ## geom_jitter() + 
  scale_color_manual(values = groupcol) 

p2 <- apply(cnvmeta %>% ## filter(patient %in% sampaired) %>% 
              dplyr::select(cluster, sample) %>% 
              table(), 2, entropy, unit='log2') %>% 
  as.data.frame() %>% dplyr::rename(Entropy='.') %>% 
  dplyr::mutate(samples[rownames(.),c('group','tissue')]) %>% 
  ggplot(aes(x=tissue, y=Entropy, color=tissue)) + 
  geom_boxplot() + ## geom_jitter() + 
  scale_color_manual(values = tissuecol) +
  stat_compare_means()

library(philentropy)
p3 <- apply(cnvmeta %>% ## filter(patient %in% sampaired) %>% 
        dplyr::select(cluster, sample) %>% 
        table(), 2, function(x){
          x <- rbind(x/sum(x), rep(1/length(x), length(x)))
          suppressMessages(JSD(x, est.prob="empirical"))
        }) %>% 
  as.data.frame() %>% dplyr::rename(JSD='.') %>% 
  dplyr::mutate(samples[rownames(.),c('group','tissue')]) %>% 
  mutate(group=factor(group, levels = c("NN","PN","PT","HM"))) %>%
  ggplot(aes(y=group, x=JSD, color=group)) + 
  geom_boxplot() + ## geom_jitter() + 
  scale_color_manual(values = groupcol)

p4 <- apply(cnvmeta %>% ## filter(patient %in% sampaired) %>% 
              dplyr::select(cluster, sample) %>% 
              table(), 2, function(x){
                x <- rbind(x/sum(x), rep(1/length(x), length(x)))
                suppressMessages(JSD(x, est.prob="empirical"))
              }) %>% 
  as.data.frame() %>% dplyr::rename(JSD='.') %>% 
  dplyr::mutate(samples[rownames(.),c('group','tissue')]) %>% 
  ggplot(aes(x=tissue, y=JSD, color=tissue)) + 
  geom_boxplot() + ## geom_jitter() + 
  scale_color_manual(values = tissuecol) +
  stat_compare_means()

pdf(paste0(fig_dir, 'Fig3.infercnv.JSD.pdf'), width = 8.27, height = 6)
print(ggarrange(p1, p3, nrow = 2, ncol = 2))
print(ggarrange(p2, p4, nrow = 2, ncol = 1))
dev.off()

genepos <- read.delim('../gene_ordering_file.txt', row.names = 1, header = F)
genepos <- genepos %>% mutate(chrom=factor(V2, levels=unique(V2))) 
chromsize <- genepos %>% group_by(chrom) %>% 
  summarise(size=as.numeric(max(V3))) %>% as.data.frame()
rownames(chromsize) <- chromsize$chrom
chromsize <- setNames(c(0, cumsum(chromsize$size) %>% head(-1)), rownames(chromsize))
## chrom chrompos abspos
genepos <- genepos %>% mutate(chrompos=round(mean(c(V3,V4)))) %>%
  mutate(abspos=chromsize[chrom] + chrompos) %>% 
  select(chrom, chrompos, abspos)

out <- readRDS('infercnv.observations.rds')

cnvgenes <- intersect(rownames(genepos), rownames(out))

chr <- as.numeric(genepos[cnvgenes, 'chrom']) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR)

my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=9, name="RdBu")))(n = 999)
col_breaks = c(seq(-1,-0.5,length=50), 
               seq(-0.5,-0.3,length=150),
               seq(-0.3,0.3,length=600), 
               seq(0.3,0.5,length=150),
               seq(0.5, 1,length=50))

## https://github.com/navinlabcode/copykat
for(i in levels(cnvmeta$cluster)){
  print(i)
  cnvdat <- cbind(
    genepos[cnvgenes, ],
    out[cnvgenes, cnvmeta %>% filter(cluster == i) %>% rownames()] 
  )
  
  png(sprintf('figures/infercnv/inferCNV.cluster%s.png', i), height=500, width=800)
  cid <- 4:ncol(cnvdat)
  if(length(cid) > 500){
    set.seed(i)
    cid <- cid %>% sample(500) %>% sort()
  }
  pdat <- apply(cnvdat[,cid], 2, scale) %>% t()
  copykat::heatmap.3(pdat, dendrogram="r", 
                     distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"), 
                     hclustfun = function(x) hclust(x, method="ward.D2"),
                     ColSideColors=chr1, Colv=NA, Rowv=TRUE, labRow="", labCol="",
                     notecol="black", col=my_palette, breaks=col_breaks, key=TRUE,
                     keysize=1, density.info="none", trace="none",
                     cexRow=0.1, cexCol=0.1, cex.main=1, cex.lab=0.1,
                     symm=F, symkey=F, symbreaks=T, cex=1, cex.main=4)
  dev.off()
}


sams <- cnvmeta %>% select(sample, tissue, group) %>% unique()
rownames(sams) <- sams$sample
cnvSample <- do.call('cbind', lapply(rownames(sams), function(i){
  cid <- cnvmeta %>% filter(sample == i) %>% rownames()
  out[, cid] %>% rowMeans()
}))
colnames(cnvSample) <- rownames(sams)

cnvhGenes <- rownames(cnvSample)[apply(cnvSample, 1, sd) > 0.025]

cnvScore <- colMeans(abs(out[cnvhGenes, ] - 1))
# cnvScore[cnvScore > 0.05] <- 0.5
cnvScorePlot <- data.frame(cnvScore) %>% mutate(cnvmeta[rownames(.),])


htm <- pheatmap(cor(cnvSample), silent = T)
sam2 <- sams[htm$tree_row$labels[htm$tree_row$order], ] %>% 
  mutate(sample=factor(sample, levels = sample), group = factor(group, levels = unique(group))) %>% 
  arrange(group, sample)
cnvcor <- cor(cnvSample)

pheatmap(cnvcor[rownames(sam2), rownames(sam2)], border_color = NA, fontsize = 5,
         cluster_rows = F, cluster_cols = F,
         annotation_row = sams %>% select(group, tissue))

pdf(paste0(fig_dir, 'Fig3.Epithelial.DimPlot.pdf'), height=8.27, width=8.27)
p <- DimPlot(seuratObj, group.by = 'primary_cluster', cols.highlight = mcol['Epithelial'], raster = F,
        cells.highlight=metadata %>% filter(primary_cluster %in% 'Epithelial') %>% rownames()) +
  NoLegend() + NoAxes()
print(p)
dev.off()
rm(p)



library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

cnvScorePlot <- cnvScorePlot %>% 
  mutate(sample=factor(sample, levels = rownames(sam2))) 

samstat <- samples$patient %>% table
sampaired <- names(samstat[samstat>1])

cnvplt <- cnvScorePlot %>% filter(patient %in% sampaired) %>%
  filter(!(tissue == "Pancreas" & group == "PN")) %>%
  mutate(group = factor(group, levels = c("PN","PT","HM"))) %>%
  mutate(tissue=factor(tissue, levels = c("Eye","Intestine","Stomach","Pancreas","Liver"))) 

p1 <- cnvplt %>% ggplot(aes(y = cnvScore, x = tissue, fill = group)) +
  geom_boxplot(outlier.colour = NA) + 
  scale_fill_manual(values = groupcol) 

plist <- lapply(levels(cnvplt$tissue), function(x){
  cnvplt %>% filter(tissue == x) %>%
    ggplot(aes(y = cnvScore, x = group, fill = group)) +
    geom_boxplot(outlier.colour = NA) + 
    scale_fill_manual(values = groupcol) +
    stat_compare_means() + theme_minimal()
})
p2 <- ggarrange(plotlist = plist, nrow = 2, ncol = 3, 
          common.legend = T, align = 'hv')

pdf(paste0(fig_dir, 'Fig3.inferCNV.score.pdf'), height=8.27, width=8.27)
print(ggarrange(p1, nrow = 2))
print(p2)
dev.off()

CHR2 <- c('white','darkred')[(cnvgenes %in% cnvhGenes)+1]
chr2 <- cbind(CHR2,CHR)

for(i in sampaired){
  print(i)
  cnvdat <- cbind(
    genepos[cnvgenes, ],
    out[cnvgenes, cnvmeta %>% filter(patient == i) %>% rownames()] 
  )
  
  png(sprintf('figures/infercnv/inferCNV.patient-%s.png', i), height=500, width=800)
  cid <- 4:ncol(cnvdat)
  if(length(cid) > 500){
    set.seed(2024)
    cid <- cid %>% sample(500) %>% sort()
  }
  cells <- rbind(groupcol[cnvmeta[colnames(cnvdat[,cid]), 'group']] %>% as.character())
  pdat <- apply(cnvdat[,cid], 2, scale) %>% t()
  copykat::heatmap.3(pdat, dendrogram="r", main = i,
                     distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"), 
                     hclustfun = function(x) hclust(x, method="ward.D2"),
                     ColSideColors=chr2, RowSideColors=cells, Colv=NA, Rowv=TRUE, labRow="", labCol="",
                     notecol="black", col=my_palette, breaks=col_breaks, key=TRUE,
                     keysize=1, density.info="none", trace="none",
                     cexRow=0.1, cexCol=0.1, cex.main=1, cex.lab=0.1,
                     symm=F, symkey=F, symbreaks=T, cex=1, cex.main=4)
  dev.off()
}


cid <- cnvmeta %>% filter(group %in% c('PT','HM') & tissue != 'Liver') %>% rownames()
cnvtest <- cbind(t(out[cnvhGenes, cid]-1), 
                 cnvmeta[cid, ] %>% select(group, tissue))

cnvm <- cnvtest %>% group_by(group, tissue) %>%
  summarise_at(cnvhGenes, mean, na.rm = TRUE)
cnvfc <- cnvm %>% group_by(tissue) %>% 
  summarise(across(!group, ~.x[group == 'HM']-.x[group == 'PT'])) %>%
  as.data.frame()
rownames(cnvfc) <- cnvfc$tissue
cnvfc <- cnvfc %>% select(-tissue)

cnvp <- cnvtest %>% group_by(tissue) %>% 
  summarise(across(!group, ~wilcox.test(.x[group == 'PT'], .x[group == 'HM'])$p.value)) %>%
  as.data.frame()
rownames(cnvp) <- cnvp$tissue
cnvp <- cnvp %>% select(-tissue) 


cnvpltd <- cnvfc[, apply(abs(cnvfc) > 0.09, 2, any)] %>% t()
phtm <- pheatmap::pheatmap(cnvpltd, silent = T)
k <- 4
annotation_row <- cutree(phtm$tree_row, k = k) %>% as.data.frame() %>%
  dplyr::rename(cluster='.') %>% mutate(cluster = sprintf("K%s", cluster))
cnvt <- data.frame(cluster=annotation_row, gene = rownames(annotation_row),
                   value=apply(cnvpltd, 1, function(x){max(abs(x))})[rownames(annotation_row)]) %>%
  group_by(cluster) %>% top_n(5, value) %>% pull(gene)
cnvt <- c(intersect(rownames(tfs), rownames(cnvpltd)), cnvt) %>% unique()
cnvt <- c(intersect(names(refGenes)[match(driverGs, refGenes)],rownames(cnvpltd)), cnvt) %>% unique()

ann_colors <- list(
  cluster = setNames(tableau_color_pal("Nuriel Stone")(k), sprintf("K%s", 1:k))
)
library(pals)
mapal <- my_palette ## brewer.piyg(256) %>% rev()
p1 <- ComplexHeatmap::pheatmap(cnvpltd, use_raster=F, border_color = NA, fontsize = 10, 
                                 labels_row = NULL, cellwidth = 10, color = mapal, cutree_rows = k,
                                 annotation_colors=ann_colors, annotation_row=annotation_row)
p2 <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=match(cnvt, rownames(cnvpltd)), 
                                                                  labels=refGenes[cnvt], 
                                                                  labels_gp=grid::gpar(fontsize = 8)))
library(org.Hs.eg.db)
library(clusterProfiler)
cnvglist <- split(rownames(annotation_row), annotation_row$cluster)
gotable <- NULL
for(k in names(cnvglist)){
  ego <- enrichGO(gene          = cnvglist[[k]],
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  gotable <- rbind(gotable, ego@result %>% mutate(cluster=k))
}

## pink_material
plotGO <- function(res, top_n=10, scale_size=c(1,6), color=paletteer_d("ggsci::teal_material")[1:6]){
  res$score <- -log10(res$p.adjust)
  res$Description <- stringr::str_trunc(res$Description, 100, "right")
  res <- res %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = score, n = top_n, with_ties=FALSE)
  
  mat <- res %>% 
    dplyr::select(-ID, -GeneRatio, -BgRatio,  -pvalue, -p.adjust, -qvalue, -geneID, -score) %>%  
    tidyr::pivot_wider(names_from = cluster, values_from = Count) %>% data.frame() 
  row.names(mat) <- mat$Description  # put gene in `row`
  mat <- mat[,-1] #drop gene column as now in rows
  mat[is.na(mat)] <- 0
  clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
  res$Description <- factor(res$Description, levels = clust$labels[clust$order])
  p <- ggplot2::ggplot(res, aes(x=cluster, y=Description)) + 
    ggplot2::geom_point(aes(size= Count, fill = score),shape=21,alpha=0.9) +
    ggplot2::scale_fill_gradientn(colours =color) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    ggplot2::scale_size_continuous(range = scale_size) + 
    ggplot2::ylab('') + ggplot2::xlab("") +
    ggplot2::labs(fill = "-log10(p.adjust)", size="Count") 
  p
}


pdf(paste0(fig_dir, 'Fig3.inferCNV.topGenes.pdf'), height=8.27, width=8.27)
print(p1 + p2)
print(plotGO(gotable))
dev.off()


cnvd <- cnvm %>% select(-group) %>% t %>% as.data.frame()
colnames(cnvd) <- cnvm %>% pull(group)
cnvd$fc <- cnvd$HM - cnvd$PT
cnvd$fdr <- p.adjust(unlist(cnvp)[rownames(cnvd)], method = 'fdr')
cnvd$name <- refGenes[rownames(cnvd)]

library(ggpubr)
sp <- ggscatter(cnvd, x = "PT", y = "HM", color = "lightgray")
sp + stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  gradient_fill("YlOrRd") + geom_abline()

ggscatter(cnvd, x = "PT", y = "HM", size = 0.3,
          palette = "jco",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(method = "spearman")

table(cnvmeta$cluster)
cnvdat <- cbind(
  genepos[cnvgenes, ],
  out[cnvgenes, slice_sample(cnvmeta %>% filter(patient %in% sampaired), 
                             n=100, by = patient) %>% rownames()] 
)
pdat <- apply(cnvdat[,4:ncol(cnvdat)], 2, scale) %>% t()

cells <- rbind(
  cluster=cnvcol[cnvmeta[rownames(pdat), 'cluster']],
  sample=samplecol[cnvmeta[rownames(pdat), 'sample']],
  tissue=tissuecol[cnvmeta[rownames(pdat), 'tissue']],
  group=groupcol[cnvmeta[rownames(pdat), 'group']]
)
colnames(cells) <- rownames(pdat)

hcc <- hclust(parallelDist::parDist(pdat, threads = 8, 
                                    method = "euclidean"), 
              method = "ward.D2")

png(paste0(fig_dir, "Fig3.inferCNV.all.png"), height=800, width=1000)
copykat::heatmap.3(pdat, dendrogram="r", main = 'Patients with paired samples', 
                   distfun = function(x) parallelDist::parDist(x, threads = 8, method = "euclidean"), 
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   ColSideColors=chr2, RowSideColors=cells, Colv=NA, Rowv=TRUE, labRow="", labCol="",
                   notecol="black", col=my_palette, breaks=col_breaks, key=TRUE,
                   keysize=1, density.info="none", trace="none",
                   cexRow=0.1, cexCol=0.1, cex.main=1, cex.lab=0.1,
                   symm=F, symkey=F, symbreaks=T, cex=1, cex.main=4)
dev.off()

cid <- hcc$labels[hcc$order]
cid <- cnvmeta[cid, ] %>% mutate(group=factor(group, levels = unique(group))) %>%
  arrange(group) %>% rownames(.)
cellcols <- (cells %>% t %>% as.data.frame())[cid, ]

png(paste0(fig_dir, "Fig3.inferCNV.all2.png"), height=600, width=1000)
copykat::heatmap.3(pdat[cid,], dendrogram="none", main = 'Patients with paired samples', 
                   # distfun = function(x) parallelDist::parDist(x, threads = 8, method = "euclidean"), 
                   # hclustfun = function(x) hclust(x, method="ward.D2"),
                   ColSideColors=chr2, RowSideColors=cells[,cid], Colv=NA, Rowv=NA, labRow="", labCol="",
                   notecol="black", col=my_palette, breaks=col_breaks, key=FALSE,
                   keysize=1, density.info="none", trace="none",
                   cexRow=0.1, cexCol=0.1, cex.main=1, cex.lab=0.1,
                   symm=F, symkey=F, symbreaks=T, cex=1, cex.main=4)
dev.off()


pdf(paste0(fig_dir, "Fig3.inferCNV.colbar.pdf"), height=6, width=8.27)
op <- par(mar=c(1,1,1,1), mfrow=c(6,1))
pals::pal.bands(cellcols$cluster)
pals::pal.bands(cellcols$sample)
pals::pal.bands(cellcols$tissue)
pals::pal.bands(cellcols$group)
pals::pal.bands(chr2[,1])
pals::pal.bands(chr2[,2])
par(op)
dev.off()


cnvScorePlot %>% ggplot(aes(cnvScore, y = sample,
                            fill = 0.5 - abs(0.5 - stat(ecdf))), color=NA) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_distiller(palette = "PiYG", name='CNV score') +
  theme_ridges() + xlim(c(0,0.04))




pdf(paste0(fig_dir, 'Fig3.inferCNV.all.pdf'), height=5, width=8.27)
pltd <- pdat[cid, ]
pltd[pltd > 5] <- 5
pltd[pltd < -5] <- -5
ComplexHeatmap::pheatmap(pltd, border_color = NA, color = my_palette,
                         cluster_rows = F, cluster_cols = F, 
                         show_rownames = F, show_colnames = F)
dev.off()