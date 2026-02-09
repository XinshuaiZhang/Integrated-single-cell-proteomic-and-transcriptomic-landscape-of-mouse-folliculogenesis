library(TimeTalk)
library(Seurat)
library(tidyverse)
library(monocle3)
library(dplyr)
library(tibble)
library(purrr)
tmp.path <- system.file("extdata/Mouse-2020-Shao-LR-pairs.txt",
                        package = "TimeTalk")
LRpairs.df <- read.delim(file = tmp.path,stringsAsFactors = F)
load('./timetalk_all.RData')

tmp.cds=cds.Oocyte
tmp.seu=scdata_pro
tmp.orig.ident = "CellType"
tmp.ident.1="Granulosa"
tmp.ident.2 = "Oocyte"
LRpairs.df = LRpairs.df
tmp.mra.res = tmp.mra.res
tmp.winsz = 0.1
tmp.lags = 1
numPts = 200
tmp.SCC.cutoff = 0.2
tmp.granger.cutoff = 1e-2
tmp.cores = 10

cat(paste0("prepare ligand and receptor genelist"),sep = "\n")
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol
tmp.data <- GetAssayData(
  object = tmp.seu, 
  slot = "data", 
  assay = "RNA"
)
gene_symbols <- rownames(tmp.data)
l.remove <- setdiff(Lgenelist, gene_symbols)
r.remove <- setdiff(Rgenelist, gene_symbols)
index.remove <- c(
  which(Lgenelist %in% l.remove), 
  which(Rgenelist %in% r.remove)
)
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]


tmp.data_ident.1 <- GetAssayData(
  object = subset(tmp.seu, subset = CellType == tmp.ident.1),
  slot = "data",
  assay = "RNA"
)
l.cellnum <- tmp.data_ident.1[Lgenelist, ] %>% 
  as.data.frame() %>% 
  apply(., 1, function(gene_list) {
    sum(gene_list > 0)
  }) %>% 
  data.frame() %>% 
  rename_all(~c('L_cellnum')) %>% 
  mutate(L = Lgenelist) %>% 
  dplyr::select(c('L','L_cellnum'))
tmp.data_ident.2 <- GetAssayData(
  object = subset(tmp.seu, subset = CellType == tmp.ident.2),
  slot = "data",
  assay = "RNA"
)
r.cellnum <- tmp.data_ident.2[Rgenelist, ] %>% 
  as.data.frame() %>% 
  apply(., 1, function(gene_list) {
    sum(gene_list > 0)
  }) %>% 
  data.frame() %>% 
  rename_all(~c('R_cellnum')) %>% 
  mutate(R = Rgenelist) %>% 
  dplyr::select(c('R','R_cellnum'))


cat(paste0(tmp.ident.1, "-", tmp.ident.2, " start:"), sep = "\n")
tmp.df <- data.frame(pseudotime = pseudotime(tmp.cds, reduction_method = "UMAP"), 
                     stringsAsFactors = F)
tmp.df <- cbind(tmp.seu[[]], tmp.df)
if (!c("CellType") %in% colnames(tmp.df)) {
  stop("Please add CellType annotation in seurat object!")
}
cat(paste0("scale pseudotime"), sep = "\n")
tmp.cell.meta.1 <- tmp.df %>% rownames_to_column("cell_id") %>% 
  dplyr::filter(CellType == tmp.ident.1) %>% arrange(pseudotime)
x <- tmp.cell.meta.1$pseudotime
tmp.cell.meta.1$pseudotime <- MinMaxScale(x)
tmp.cell.meta.2 <- tmp.df %>% rownames_to_column("cell_id") %>% 
  dplyr::filter(CellType == tmp.ident.2) %>% arrange(pseudotime)
x <- tmp.cell.meta.2$pseudotime
tmp.cell.meta.2$pseudotime <- MinMaxScale(x)
tmp.mat.1 <- tmp.data[, tmp.cell.meta.1$cell_id]
tmp.mat.2 <- tmp.data[, tmp.cell.meta.2$cell_id]
tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
inter.tmp.mat.1 <- cellAlign::interWeights(
  expDataBatch = tmp.mat.1, 
  trajCond = tmp.mat.pseudotime.1, 
  winSz = tmp.winsz, 
  numPts = numPts)
inter.tmp.mat.2 <- cellAlign::interWeights(
  expDataBatch = tmp.mat.2,
  trajCond = tmp.mat.pseudotime.2,
  winSz = tmp.winsz,
  numPts = numPts
)
inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
time <- inter.tmp.mat.1$traj
inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
plan("multisession", workers = tmp.cores)
tmp.IS.list <- future.apply::future_lapply(seq_along(Lgenelist), FUN = function(ii) {
  cat(ii, sep = "\n")
  x <- inter.tmp.mat.1[Lgenelist[ii], ]
  y <- inter.tmp.mat.2[Rgenelist[ii], ]
  tmp.IS <- sqrt(x * y)
  return(tmp.IS)
})
names(tmp.IS.list) <- paste0(Lgenelist, "-", Rgenelist)
tmp.IS.df <- Reduce(rbind, tmp.IS.list) %>% 
  as.data.frame() %>% 
  mutate(
    LRpairs = paste0(Lgenelist, "-", Rgenelist)
  ) %>% 
  rownames_to_column(var = "temp_rowname") %>%
  dplyr::select(-temp_rowname) %>% 
  column_to_rownames('LRpairs')
colnames(tmp.IS.df) <- time %>% as.character()
data_si_use <- tmp.IS.df %>% 
  rownames_to_column('LRpairs') %>% 
  filter(LRpairs %in% (tmp.res.df_gc$LR %>% unique())) %>%
  column_to_rownames('LRpairs')
pheatmap::pheatmap(
  mat = data_si_use,
  scale = 'none',
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE
)
save.image('./cci_IS_G-to-O.RData')  
saveRDS(data_si_use,'./data_si_use_G-to-O.rds')
library(Mfuzz)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
data_si_use <- readRDS('./data_si_use_G-to-O.rds')

DEGs_exp_averp <- data_si_use %>% as.matrix()
boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- Mfuzz::filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 6
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)


acore.list <- acore(dat,cl=cl,min.acore=0)
tmp_memship <- Reduce(rbind,acore.list) %>% 
  rename_all(~c('Gene_symbol','MEM.SHIP'))
tmp_memship$Gene_symbol %>% duplicated() %>% table()

gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars <-
  colnames(tmp)[c(-length(colnames(tmp)), -length(colnames(tmp)) + 1)]
data_plot <-
  melt(
    tmp,
    id.vars = c('Gene_symbol', 'cluster'),
    variable.name = 'time',
    measure.vars = measure_vars,
    value.name = 'values'
  ) %>% 
  left_join(
    tmp_memship,by = 'Gene_symbol'
  )

plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  colo <- colorRampPalette(rev(c("#d62728","grey92")))(30)
  colorseq <- seq(
    min(data_plot$MEM.SHIP),
    max(data_plot$MEM.SHIP), 
    length = length(colo)
  )
  for (jj in 1:(length(colorseq) - 1)) {
    tmpcol <- (data_plot$MEM.SHIP >= colorseq[jj] & data_plot$MEM.SHIP <= 
                 colorseq[jj + 1])
    data_plot[which(tmpcol),'group'] <- jj
  }
  data_plot$group <- factor(data_plot$group,levels = sort(unique(data_plot$group)))
  
  p<-ggplot(data_plot, aes(x=time, y=(values), group=Gene_symbol)) + 
    geom_line(aes(color=as.factor(group)),linewidth=1,alpha=1) +
    scale_color_manual(values = colo[levels(data_plot$group) %>% as.numeric()])+
    scale_linetype_manual(values = 'dashed') + 
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:3){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,3,1)[-seq(1,3,3)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
p_mfuzz
data_plot <- data_plot %>% 
  mutate(
    cluster_raw = cluster,
    cluster = case_when(
      cluster_raw %in% c(5) ~ 1,
      cluster_raw %in% c(1,3,4) ~ 2,
      cluster_raw %in% c(2,6) ~ 3
    )
  )
p_mfuzz <- plot_mfuzz(data_plot)
p_mfuzz
gene_list <- data_plot %>% 
  select(c('Gene_symbol','cluster')) %>% 
  distinct() %>% 
  tidyr::separate(Gene_symbol,sep = '-',into = c('Lgene','Rgene'))
for(i in (gene_list$cluster %>% unique())){
  library(clusterProfiler)
  library(stringr)
  print(i)
  genes<- c(
    gene_list$Lgene[gene_list$cluster == i],
    gene_list$Rgene[gene_list$cluster == i]
  ) %>% unique()
  gene<-str_to_title(genes)
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  ego <- enrichKEGG(
    organism ="mmu",
    gene = gene$ENTREZID,
    minGSSize = 1,
    pvalueCutoff = 1,
    qvalueCutoff = 1) 
  dotplot(ego,showCategory=30,title=paste("Enrichment KEGG Cluster",i),font.size = 20)
  file_name <- paste(
    "./Granulosa_KEGG_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")
  data_go <- ego@result
  file_name <- paste(
    "./Granulosa_KEGG_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}

library(tibble)
data_plot %>% colnames()
data_save <- data_plot %>% 
  select(c('Gene_symbol','cluster','MEM.SHIP')) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster,levels = c(1,2,3))) %>% 
  arrange(cluster,desc(MEM.SHIP))
write.csv(data_save,'./Granulosa_cci_mfuzz_group.csv')
LRpairs_order <- data_plot %>% 
  select(c('Gene_symbol','cluster','MEM.SHIP')) %>% 
  distinct() %>% 
  mutate(cluster = factor(cluster,levels = c(1,2,3))) %>% 
  arrange(cluster,desc(MEM.SHIP)) %>% pull(Gene_symbol)
data_si_use <- data_si_use %>% 
  rownames_to_column('Gene_symbol') %>% 
  mutate(
    Gene_symbol = factor(Gene_symbol,levels = LRpairs_order)
  ) %>% 
  arrange(Gene_symbol) %>% 
  column_to_rownames('Gene_symbol')
pheatmap::pheatmap(
  mat = data_si_use,
  scale = 'none',
  show_rownames = FALSE,show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE
)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
color_function <- colorRamp2(
  breaks = c(min(data_si_use %>% colnames() %>% as.numeric()), 
             median(data_si_use %>% colnames() %>% as.numeric()), 
             max(data_si_use %>% colnames() %>% as.numeric())),
  colors = c("grey93", "#7EAED3FF", "#2A5783FF")
)

ha_top = HeatmapAnnotation(
  Pseudotime = anno_simple(
    data_si_use %>% colnames() %>% as.numeric(),
    col = color_function
  )
)
mark_gene <- c(
  'Jam3-Itgb1','Psen2-Notch2','Hsp90aa1-Egfr',
  'Amh-Amhr2','Inhba-Tgfbr1'
)
gene_pos <- which(rownames(data_si_use) %in% mark_gene)
right_annotation <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                         labels = mark_gene))


p <- Heatmap(
  data_si_use %>% as.matrix(),
  name = "mat",
  col = colorRampPalette(c("#716CAC", "white", "#c2473b"))(100),
  row_gap = unit(3,units = 'pt'),
  top_annotation = ha_top,
  right_annotation = right_annotation,
  column_title = NULL,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_labels = ifelse(
    rownames(data_si_use) %in% c(mark_gene), 
    rownames(data_si_use), ""),
  row_names_gp = gpar(fontsize = 12),
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(
    title = 'Interaction Score',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
    title_position = "leftcenter-rot",
    legend_direction = 'vertical',
    legend_position = "topright",
    legend_height = unit(30,units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
p
pdf('./heatmap_Granulosa.pdf')
p
dev.off()

tmp.res <- read.csv('./res_cci_time_talk_G-O.csv',row.names = 1) %>% 
  filter(category == 'PASS')
tmp.res$LR %>% table() %>% length()
tmp.res$TF %>% table() %>% length()
colnames(tmp.res)
tmp.res
TF_select <- tmp.res$TF %>% 
  table() %>% sort() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  pull('.') %>% .[1:30] %>% as.character()
LR_select <- tmp.res$LR %>% 
  table() %>% sort() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  pull('.') %>% .[1:30] %>% as.character()
data_plot <- tmp.res %>% 
  filter(LR %in% LR_select,TF %in% TF_select) %>% 
  select(c('L','R','TF','Interaction_Score','PCC')) %>% 
  arrange(L,R,TF) %>% 
  mutate(
    L = as.factor(L),
    R = as.factor(R),
    TF = as.factor(TF)
  )
data_plot$L %>% unique()
data_plot$R %>% unique()
data_plot$TF %>% unique()
library(networkD3)

nodes <- data.frame(
  name = unique(c(as.character(data_plot$L), as.character(data_plot$R), as.character(data_plot$TF)))
)

node_index <- function(node_name) {
  match(node_name, nodes$name) - 1
}

links <- data.frame(
  source = c(node_index(as.character(data_plot$L)), node_index(as.character(data_plot$R))),
  target = c(node_index(as.character(data_plot$R)), node_index(as.character(data_plot$TF))),
  value = c(data_plot$Interaction_Score,-log10(abs(data_plot$PCC)))
)
energy <- list(nodes = nodes, links = links)
energy$nodes$Node_Group <- c(
  rep('L',as.character(data_plot$L) %>% unique() %>% length()),
  rep('R',as.character(data_plot$R) %>% unique() %>% length()),
  rep('TF',as.character(data_plot$TF) %>% unique() %>% length())
)
energy$links$Link_Group <- c(
  rep('LR',nrow(data_plot)),
  rep('LR-TF',nrow(data_plot))
)
sankeyNetwork(
  Links = energy$links,
  Nodes = energy$nodes,
  NodeGroup = 'Node_Group',
  LinkGroup = 'Link_Group',
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "TWh",iterations = 3000,
  fontSize = 20,
  nodePadding = 10,
  nodeWidth = 90
)

library(future.apply)
library(TimeTalk)
library(Seurat)
library(tidyverse)
library(monocle3)
library(dplyr)
library(tibble)
library(purrr)
setwd('./result_0815_Fig6-CCI/')
RunTimeTalk_self <- function(
    tmp.cds, tmp.seu, tmp.ident.1, tmp.ident.2, cellnum.filter = 10,LRpairs.df, 
    tmp.mra.res, tmp.orig.ident, tmp.lags = 1, tmp.winsz = 0.1, 
    numPts = 200, tmp.cores = 10, tmp.SCC.cutoff = 0.2, tmp.granger.cutoff = 0.01
){
  cat(paste0("prepare ligand and receptor genelist"),sep = "\n")
  LRpairs <- LRpairs.df$lr_pair
  Lgenelist <- LRpairs.df$ligand_gene_symbol
  Rgenelist <- LRpairs.df$receptor_gene_symbol
  tmp.data <- GetAssayData(
    object = tmp.seu, 
    slot = "data", 
    assay = "RNA"
  )
  gene_symbols <- rownames(tmp.data)
  l.remove <- setdiff(Lgenelist, gene_symbols)
  r.remove <- setdiff(Rgenelist, gene_symbols)
  index.remove <- c(
    which(Lgenelist %in% l.remove), 
    which(Rgenelist %in% r.remove)
  )
  LRpairs <- LRpairs[-index.remove]
  Lgenelist <- Lgenelist[-index.remove]
  Rgenelist <- Rgenelist[-index.remove]
  
  tmp.data_ident.1 <- GetAssayData(
    object = subset(tmp.seu, subset = CellType == tmp.ident.1),
    slot = "data",
    assay = "RNA"
  )
  l.remove <- tmp.data_ident.1[Lgenelist, ] %>% 
    as.data.frame() %>% 
    apply(., 1, function(gene_list) {
      sum(gene_list > 0)
    }) %>% 
    data.frame() %>% 
    rename_all(~c('exp_cellnum')) %>% 
    mutate(Lgenelist = Lgenelist) %>% 
    filter(exp_cellnum<cellnum.filter) %>% 
    pull(Lgenelist) %>% unique()
  tmp.data_ident.2 <- GetAssayData(
    object = subset(tmp.seu, subset = CellType == tmp.ident.2),
    slot = "data",
    assay = "RNA"
  )
  r.remove <- tmp.data_ident.2[Rgenelist, ] %>% 
    as.data.frame() %>% 
    apply(., 1, function(gene_list) {
      sum(gene_list > 0)
    }) %>% 
    data.frame() %>% 
    rename_all(~c('exp_cellnum')) %>% 
    mutate(Rgenelist = Rgenelist) %>% 
    filter(exp_cellnum<cellnum.filter) %>% 
    pull(Rgenelist) %>% unique()  
  index.remove <- c(
    which(Lgenelist %in% l.remove), 
    which(Rgenelist %in% r.remove)
  )
  LRpairs <- LRpairs[-index.remove]
  Lgenelist <- Lgenelist[-index.remove]
  Rgenelist <- Rgenelist[-index.remove]
  
  cat(paste0(tmp.ident.1, "-", tmp.ident.2, " start:"), sep = "\n")
  tmp.df <- data.frame(pseudotime = pseudotime(tmp.cds, reduction_method = "UMAP"), 
                       stringsAsFactors = F)
  tmp.df <- cbind(tmp.seu[[]], tmp.df)
  if (!c("CellType") %in% colnames(tmp.df)) {
    stop("Please add CellType annotation in seurat object!")
  }
  cat(paste0("scale pseudotime"), sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>% rownames_to_column("cell_id") %>% 
    dplyr::filter(CellType == tmp.ident.1) %>% arrange(pseudotime)
  x <- tmp.cell.meta.1$pseudotime
  tmp.cell.meta.1$pseudotime <- MinMaxScale(x)
  tmp.cell.meta.2 <- tmp.df %>% rownames_to_column("cell_id") %>% 
    dplyr::filter(CellType == tmp.ident.2) %>% arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime
  tmp.cell.meta.2$pseudotime <- MinMaxScale(x)
  tmp.mat.1 <- tmp.data[, tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[, tmp.cell.meta.2$cell_id]
  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
  inter.tmp.mat.1 <- cellAlign::interWeights(expDataBatch = tmp.mat.1, 
                                             trajCond = tmp.mat.pseudotime.1, winSz = tmp.winsz, 
                                             numPts = numPts)
  inter.tmp.mat.2 <- cellAlign::interWeights(expDataBatch = tmp.mat.2, 
                                             trajCond = tmp.mat.pseudotime.2, winSz = tmp.winsz, 
                                             numPts = numPts)
  inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
  plan("multisession", workers = tmp.cores)
  tmp.res.list <- future.apply::future_lapply(seq_along(Lgenelist), FUN = function(ii) {
    cat(ii, sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii], ]
    y <- inter.tmp.mat.2[Rgenelist[ii], ]
    tmp.res <- cor(x, y, method = "pearson")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist, "-", Rgenelist)
  tmp.res.list.PCC <- tmp.res.list
  tmp.res.list <- future_lapply(seq_along(Lgenelist), FUN = function(ii) {
    cat(ii, sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii], ]
    y <- inter.tmp.mat.2[Rgenelist[ii], ]
    tmp.res <- cor(x, y, method = "spearman")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist, "-", Rgenelist)
  tmp.res.list.SCC <- tmp.res.list
  plan("sequential")
  tmp.res.cor <- data.frame(
    PCC = unlist(tmp.res.list), 
    SCC = unlist(tmp.res.list),
    stringsAsFactors = F) %>% 
    rownames_to_column("LRpairs") %>% 
    mutate(PCC = ifelse(is.na(PCC), 0, PCC)) %>% 
    mutate(SCC = ifelse(is.na(SCC), 0, SCC)) %>% 
    arrange(-SCC) %>% mutate(Rank = row_number())
  tmp.LR.list <- tmp.res.cor %>% 
    filter(abs(SCC) > tmp.SCC.cutoff) %>% 
    pull(LRpairs)
  tmp.TF.gene <- tmp.mra.res %>% 
    filter(group == tmp.ident.2) %>% 
    pull(Regulon)
  plan("multisession", workers = tmp.cores)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR) {
    cat(tmp.LR, sep = "\n")
    tmp.L.gene <- unlist(lapply(
      strsplit(tmp.LR, split = "-"),
      FUN = function(ii) {
        ii[1]
      }
    ))
    tmp.R.gene <- unlist(lapply(
      strsplit(tmp.LR, split = "-"),
      FUN = function(ii) {
        ii[2]
      }
    ))
    x <- inter.tmp.mat.1[tmp.L.gene, ]
    y <- inter.tmp.mat.2[tmp.R.gene, ]
    tmp.res.list <- lapply(tmp.TF.gene, function(tmp.TF.gene.use) {
      cat(paste0(tmp.ident.2, "_TF:", tmp.TF.gene.use), 
          sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use,]
      tmp.IS <- sqrt(x * y)
      tmp.res.LRtoTF.pvalue <- tryCatch(expr = {
        tmp.res <- lmtest::grangertest(tmp.IS, tmp.TF.level, 
                                       order = tmp.lags)
        tmp.res$`Pr(>F)`[2]
      }, error = function(e) {
        1
      })
      tmp.res.TFtoLR.pvalue <- tryCatch(expr = {
        tmp.res <- lmtest::grangertest(tmp.TF.level, tmp.IS, 
                                       order = tmp.lags)
        tmp.res$`Pr(>F)`[2]
      }, error = function(e) {
        1
      })
      tmp.PCC <- cor(tmp.IS, tmp.TF.level, method = "pearson")
      tmp.PCC <- ifelse(is.na(tmp.PCC), 0, tmp.PCC)
      tmp.SCC <- cor(tmp.IS, tmp.TF.level, method = "spearman")
      tmp.SCC <- ifelse(is.na(tmp.SCC), 0, tmp.SCC)
      tmp.res.df <- data.frame(
        LR = tmp.LR, L = tmp.L.gene,R = tmp.R.gene,Interaction_Score = mean(tmp.IS),
        TF = tmp.TF.gene.use, PCC = tmp.PCC, SCC = tmp.SCC, 
        LRtoTF = tmp.res.LRtoTF.pvalue, TFtoLR = tmp.res.TFtoLR.pvalue, 
        cell.L = tmp.ident.1, cell.R = tmp.ident.2, stringsAsFactors = F)
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind, tmp.res.list)
    tmp.res.df <- tmp.res.df %>% 
      mutate(category = ifelse(
        LRtoTF < tmp.granger.cutoff | TFtoLR < tmp.granger.cutoff, 
        "PASS", "SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  tmp.ttt.res.df <- Reduce(rbind, tmp.ttt.res)
  tmp.res.df <- tmp.ttt.res.df %>% 
    left_join(
      tmp.res.cor %>% rename_all(~c('LR','PCC_LR','SCC_LR','Rank')),
      by = 'LR'
    )
  return(tmp.res.df)
}

data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  dplyr::rename(protein = Accessions)
scdata_pro <- readRDS('./scdata_use_filter_outliers.rds')

meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')
data_singlestage_raw <- FetchData(
  scdata_pro,layer = 'counts',
  vars = rownames(scdata_pro)
) %>% 
  t() %>% as.data.frame()
tmp <- data_singlestage_raw[rownames(data_singlestage_raw) == 'Q9D2S9',] %>% as.numeric()
data_singlestage <-
  data_singlestage_raw %>%
  rownames_to_column('Accessions') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  dplyr::select(-c('Accessions')) %>%
  as.data.frame() %>%
  group_by(Genes) %>%
  summarise(across(.cols = everything(), .fns = ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup() %>%
  column_to_rownames(var = 'Genes')
data_singlestage[rownames(data_singlestage) == 'Lsm11',]
scdata_pro <- CreateSeuratObject(
  counts = data_singlestage %>% as.matrix(),
  data = data_singlestage %>% as.matrix(),
  meta.data = scdata_pro@meta.data
)
scdata_pro$CellType <- scdata_pro$group_cell


tmp.path <- system.file("extdata/Mouse-2020-Shao-LR-pairs.txt",
                        package = "TimeTalk")
LRpairs.df <- read.delim(file = tmp.path,stringsAsFactors = F)
tmp.mra.res <- readRDS('./timetalk_out/B_blastoid_RTN_mra_result.rds')
cds.Oocyte <- readRDS('./phylovelo_output/cds_monocle3.rds')
tmp.cds=cds.Oocyte
tmp.seu=scdata_pro
tmp.orig.ident = "CellType"
tmp.ident.1="Oocyte"
tmp.ident.2 = "Granulosa"
LRpairs.df = LRpairs.df
tmp.mra.res = tmp.mra.res
tmp.winsz = 0.1
tmp.lags = 1
numPts = 200
tmp.SCC.cutoff = 0.2
tmp.granger.cutoff = 1e-2
tmp.cores = 10
tmp.res.df <- RunTimeTalk_self(
  tmp.cds = tmp.cds,tmp.seu = tmp.seu,
  tmp.ident.1 = tmp.ident.1,
  tmp.ident.2 = tmp.ident.2,
  cellnum.filter = 5,LRpairs.df = LRpairs.df,
  tmp.mra.res = tmp.mra.res,
  tmp.orig.ident = tmp.orig.ident,tmp.lags = tmp.lags,
  tmp.winsz = tmp.winsz,numPts = numPts,tmp.cores = tmp.cores,
  tmp.SCC.cutoff = tmp.SCC.cutoff,tmp.granger.cutoff = tmp.granger.cutoff
)
write.csv(tmp.res.df,'./res_cci_time_talk_O-G.csv')
save.image('./timetalk_Oocyte.RData')
tmp.res.df_gc <- RunTimeTalk_self(
  tmp.cds = tmp.cds,tmp.seu = tmp.seu,
  tmp.ident.1 = tmp.ident.2,
  tmp.ident.2 = tmp.ident.1,
  cellnum.filter = 5,LRpairs.df = LRpairs.df,
  tmp.mra.res = tmp.mra.res,
  tmp.orig.ident = tmp.orig.ident,tmp.lags = tmp.lags,
  tmp.winsz = tmp.winsz,numPts = numPts,tmp.cores = tmp.cores,
  tmp.SCC.cutoff = tmp.SCC.cutoff,tmp.granger.cutoff = tmp.granger.cutoff
)
write.csv(tmp.res.df_gc,'./res_cci_time_talk_G-O.csv')
save.image('./timetalk_all.RData')

library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(Seurat)
library(stringr)
data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  dplyr::rename(protein = Accessions)
scdata_pro <- readRDS('./scdata_use_filter_outliers.rds')
rownames(scdata_pro) %>% head()
rownames(scdata_pro) %>% length()
scdata_pro@meta.data %>% colnames()
res_phylo_pseudotime <- read.csv('./res_phylo_pseudotime.csv') %>% 
  rename(sample = names(.)[1])
gene_list <- c(
  'Lta4h','Dpysl4','Rab27b'
) %>% 
  str_to_title() %>% unique()
protein_list <- data_gene_protein %>% 
  mutate(
    Genes = factor(Genes,levels = Genes %>% unique())
  ) %>%
  arrange(Genes)
data_tmp <- scdata_pro@assays$RNA$data %>% as.data.frame() %>% 
  rownames_to_column('Accessions') %>% 
  filter(
    Accessions %in% protein_list$Accessions
  ) %>% 
  left_join(protein_list,by = 'Accessions') %>% 
  select(
    Accessions,Genes,
    names(.)[!(names(.) %in% c('Accessions','Genes'))]
  )
data_plot <- data_plot <- data_tmp %>%
  select(-Accessions) %>%
  group_by(Genes) %>%
  mutate(dup_i = row_number()) %>%
  ungroup() %>%
  mutate(Genes = if_else(dup_i == 1L, Genes, paste0(Genes, "_", dup_i))) %>%
  select(-dup_i) %>%
  column_to_rownames("Genes") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('sample') %>%
  left_join(
    scdata_pro@meta.data %>% select(group_cell,group_stage,group_stage_cell) %>% rownames_to_column('sample'),by = 'sample'
  ) %>% 
  pivot_longer(cols = names(.)[2:(ncol(.)-3)],names_to = 'Genes',values_to = 'copy_number') %>% 
  mutate(
    group_stage = factor(group_stage,levels = c(
      'Secondary','Early antral','Antral','Preovulatory'
    )),
    group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary: Oocyte","Early antral: Oocyte",
      "Antral: Oocyte","Preovulatory: Oocyte", 
      "Secondary: Granulosa","Early antral: Granulosa",
      "Antral: Granulosa","Preovulatory: Granulosa"
    ))
  ) %>% 
  mutate(Genes = factor(Genes,levels = protein_list$Genes %>% levels())) %>% 
  arrange(Genes,group_cell,group_stage) %>% 
  filter(
    group_cell == 'Oocyte'
  ) %>% 
  group_by(Genes,group_cell) %>% 
  mutate(
    rela_copy_number = scale(copy_number,center = FALSE,scale = TRUE)
  ) %>% 
  filter(!is.na(rela_copy_number)) %>% 
  left_join(res_phylo_pseudotime %>% select(-c('group_stage_cell')),by = 'sample') %>% 
  group_by(group_stage_cell) %>% 
  arrange(pseudotime) %>% 
  mutate(sample = factor(sample,levels = unique(sample)))
colnames(data_plot)
stat_df <- data_plot %>%
  group_by(Genes) %>%
  summarize(
    n = sum(is.finite(pseudotime) & is.finite(rela_copy_number)),
    rho = suppressWarnings(cor(pseudotime, rela_copy_number, method = "spearman", use = "pairwise.complete.obs")),
    p_value = {
      x <- pseudotime
      y <- rela_copy_number
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3) NA_real_ else suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    star = case_when(
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label = sprintf("Spearman \u03c1=%.2f\nFDR=%.2g%s\nn=%d", rho, FDR, star, n)
  ) %>% 
  arrange(desc(rho)) %>% 
  mutate(order_num = 1:n())
stat_df <- data_plot %>%
  group_by(Genes) %>%
  summarize(
    n = sum(is.finite(pseudotime) & is.finite(rela_copy_number)),
    rho = suppressWarnings(cor(pseudotime, rela_copy_number, method = "spearman", use = "pairwise.complete.obs")),
    p_value = {
      x <- pseudotime
      y <- rela_copy_number
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3) NA_real_ else suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    
    FDR_top5_cutoff = quantile(FDR, probs = 0.05, na.rm = TRUE),
    FDR_top5 = !is.na(FDR) & FDR <= FDR_top5_cutoff,
    
    star = case_when(
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label = sprintf("Spearman \u03c1=%.2f\nFDR=%.2g%s\nn=%d", rho, FDR, star, n)
  ) %>%
  arrange(desc(rho)) %>%
  mutate(order_num = 1:n())

getwd()
write.csv(stat_df,'./stat_df.csv')
stat_df %>% 
  filter(Genes %in% c('Ooep','Zp2','Gdf9','Dnmt3l'))
