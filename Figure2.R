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
