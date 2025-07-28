library(dplyr)
# 相关R包下载与载入：
library(tidyverse)  # 载入tidyverse包，包括多个数据处理和绘图包
library(ggsankey)   # 载入ggsankey包，用于创建桑基图
library(ggplot2)    # 载入ggplot2包，用于创建绘图
library(cols4all)   # 载入cols4all包，用于定制颜色

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
# cnvObj <- CreateSeuratObject(counts = out)
saveRDS(cnvObj, 'infercnv.seuratObj.rds')
# Idents(cnvObj) <- 'cluster'
# all.markers <- FindAllMarkers(object = cnvObj)
# 
# cnvObj@misc$cluster.marker <- all.markers
# saveRDS(cnvObj, 'infercnv.seuratObj.rds')

# cnvObj <- NormalizeData(object = cnvObj, normalization.method = 'LogNormalize')
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


# library(Nebulosa)
# plot_density(cnvObj, "ENSG00000283627", reduction = 'umap')
FeaturePlot(cnvObj, features = 'ENSG00000283627', 
            cols = c('grey', 'orangered1', 'red'))


load('../color.RData')
samples <- metadata$sample %>% sort() %>% unique()
scol <- setNames(colorRampPalette(palettes_d$ggthemes$Classic_Cyclic)(length(samples)), samples)

# 指定因子，调整显示顺序：
# 将df数据框中的'node'列转换为因子，并指定显示顺序
nodes <- c(levels(cnvmeta$subtype),
           levels(cnvmeta$cluster),
           unique(cnvmeta$sample), 
           unique(cnvmeta$group), 
           unique(cnvmeta$tissue))
cnvcol <- setNames(mycols[1:nlevels(cnvmeta$cluster)], 
                   levels(cnvmeta$cluster))
ncols <- c(xcol, cnvcol, scol,
           groupcol, tissuecol)
showlabs <- c(unique(cnvmeta$group), 
              unique(cnvmeta$tissue))

library(ggplot2)
library(ggsankey)
pdf('infercnv.sankey.pdf', width = 15, height = 8.27)
cnvmeta %>% select(subtype, cluster, sample, group, tissue) %>%
  make_long(subtype, cluster, sample, group, tissue) %>% 
  mutate(label=ifelse(node %in% showlabs, node, NA)) %>%
  mutate(node = factor(node, levels=nodes)) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node, 
             label = label,
             fill = node)) +
  geom_sankey(
    flow.alpha = 0.9,      # 桑基条带的不透明度
    # space = 50,          # 桑基节点间的距离
    # smooth = 6,          # 桑基条带的弯曲度
    # width = 0.1, alpha = 1 # 桑基节点的宽度和不透明度
  ) +  
  geom_sankey_label(size = 3, color = 'white', fill = NA) +
  theme_sankey(base_size = 18) +
  scale_fill_manual(values = ncols) + 
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Cancer cells")
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
  
  png(sprintf('inferCNV.cluster%s.png', i), height=500, width=800)
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


table(cnvmeta$cluster)
cnvdat <- cbind(
  genepos[cnvgenes, ],
  out[cnvgenes, slice_sample(cnvmeta, n=150, by = cluster) %>% rownames()] 
)
pdat <- apply(cnvdat[,4:ncol(cnvdat)], 2, scale) %>% t()

cells <- rbind(
  cluster=cnvcol[cnvmeta[rownames(pdat), 'cluster']],
  sample=scol[cnvmeta[rownames(pdat), 'sample']],
  tissue=tissuecol[cnvmeta[rownames(pdat), 'tissue']],
  group=groupcol[cnvmeta[rownames(pdat), 'group']]
)
colnames(cells) <- rownames(pdat)

hcc <- hclust(parallelDist::parDist(pdat, threads = 8, 
                                    method = "euclidean"), 
              method = "ward.D2")

png('inferCNV.all.png', height=600, width=1000)
copykat::heatmap.3(pdat, dendrogram="row", 
                   distfun = function(x) parallelDist::parDist(x, threads = 8, method = "euclidean"), 
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   ColSideColors=chr1, RowSideColors=cells, Colv=NA, Rowv=TRUE, labRow="", labCol="",
                   notecol="black", col=my_palette, breaks=col_breaks, key=TRUE,
                   keysize=1, density.info="none", trace="none",
                   cexRow=0.1, cexCol=0.1, cex.main=1, cex.lab=0.1,
                   symm=F, symkey=F, symbreaks=T, cex=1, cex.main=4)
dev.off()

cid <- hcc$labels[hcc$order]
cellcols <- (cells %>% t %>% as.data.frame())
png('inferCNV.colbar.pdf', height=6, width=8.27)

dev.off()



png('inferCNV.all.pdf', height=6, width=8.27)
ComplexHeatmap::pheatmap(pdat[cid, ], 
                         cluster_rows = F, cluster_cols = F, 
                         show_rownames = F, show_colnames = F)
dev.off()

args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)
print(args)
project <- args[1]

options(error = function() traceback(2))

library(Seurat)
library(infercnv)
library(dplyr)
rds <- sprintf('rds/%s_seurat0.rds', project)
ref_group_names <- "Mac/Mono"
if(file.exists(rds)){
  cells <- readRDS('infercnv/cell.metadata.rds') 
  x <- readRDS(rds)
  cids <- intersect(Cells(x), rownames(cells))
  annotations_file <- sprintf("infercnv/%s.cell.metadata.txt", project)
  write.table(cells[cids, ], annotations_file, sep="\t", 
              col.names = F, row.names = F, quote = F)
  counts_matrix <- x@assays$RNA@counts[,cids]
}else{
  counts_matrix <- readRDS("infercnv/counts_matrix.rds")
  annotations_file <- "infercnv/cell.metadata.txt"
  project <- "integrated"
}

out_dir <- sprintf("infercnv/%s_out", project)

# create the infercnv object
infercnvObj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix), 
                                    annotations_file=annotations_file,
                                    delim="\t",
                                    gene_order_file="gene_ordering_file.txt",
                                    ref_group_names=ref_group_names)

# infercnvObj@expr.data <- as.matrix(infercnvObj@expr.data)
# infercnvObj@count.data <- as.matrix(infercnvObj@count.data)

# perform infercnv operations to reveal cnv signal
infercnvObj <- infercnv::run(infercnvObj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=TRUE,
                             analysis_mode = "subclusters", 
                             num_threads=10,
                             output_format='pdf'
)
saveRDS(infercnvObj, sprintf('infercnv/%s_obj.rds', project))




q()
########## 
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)
print(args)
project <- args[1]

options(error = function() traceback(2))


rds <- sprintf('rds/%s_seurat0.rds', project)
if(file.exists(rds)){
  library(Seurat)
  library(infercnv)
  
  cells <- readRDS('infercnv/cell.metadata.rds')
  x <- readRDS(rds)
  cids <- intersect(Cells(x), rownames(cells))
  annotations_file <- sprintf("infercnv/%s.cell.metadata.txt", project)
  write.table(cells[cids, ], annotations_file, sep="\t", 
              col.names = F, row.names = F, quote = F)
  counts_matrix <- x@assays$RNA@counts[,cids]
  # create the infercnv object
  infercnvObj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix), 
                                       annotations_file=annotations_file,
                                       delim="\t",
                                       gene_order_file="gene_ordering_file.txt",
                                       ref_group_names=NULL)
  
  # infercnvObj@expr.data <- as.matrix(infercnvObj@expr.data)
  # infercnvObj@count.data <- as.matrix(infercnvObj@count.data)
  out_dir <- sprintf("infercnv/%s", project)
  # perform infercnv operations to reveal cnv signal
  infercnvObj <- infercnv::run(infercnvObj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=out_dir, 
                               analysis_mode = "subclusters", 
                               cluster_by_groups=TRUE, 
                               plot_steps=FALSE,
                               denoise=TRUE,
                               HMM=TRUE,
                               num_threads=10,
                               output_format='pdf'
  )
  saveRDS(infercnvObj, sprintf('infercnv/%s_obj.rds', project))
}