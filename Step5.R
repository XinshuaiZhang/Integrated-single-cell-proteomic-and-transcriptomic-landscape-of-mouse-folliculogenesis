library(Seurat)
library(Biobase)
library(dplyr)
setwd('./data_tf/')
scdata <- readRDS('./scdata_tpm.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata_sce <- as.SingleCellExperiment(scdata)
exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
meta <- scdata@meta.data
metadata <- data.frame(
  labelDescription = colnames(meta),
  row.names = colnames(meta)
)
phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                 phenoData=phenoData,
                                 annotation="scdata")
myExpressionSet

regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

print(regul)




signature <- rowTtest(myExpressionSet, "group_cell", 'Granulosa', "Oocyte")

head(signature[[2]])[1:6, ]

signature <-
  (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
     sign(signature$statistic))[,1]
head(signature)[1:6]
signature_sub <- signature[!is.na(signature)]




nullmodel <-
  ttestNull(
    myExpressionSet,
    "group_cell",
    'Granulosa',
    "Oocyte",
    per = 1000,
    repos = TRUE,
    verbose = FALSE
  )

head(nullmodel)[, 1:6]


regul

mrs <-
  msviper(
    ges = signature_sub,
    regulon =  regul,
    nullmodel =  nullmodel,
    verbose = T
  )

tmp <- summary(mrs,mrs = mrs$signature %>% length())
view(tmp)

plot(
  mrs,
  cex = 1,
  mrs = 5,
  density = 1,
  include = c("expression", "activity"),
  gama = 2
)
par(mai = c(0,0,0.5,1))
plot(
  mrs,
  cex = 2,
  mrs = 20,
  density = 1,
  include = c("expression", "activity"),
  gama = 2
)





mrs <- ledge(mrs)

summary(mrs)






vpres <- viper(myExpressionSet, regul, verbose = FALSE)

dim(vpres)

tmp <- rowTtest(vpres,  "group_cell", 'Granulosa', "Oocyte")

data.frame(
  Gene = rownames(tmp$p.value),
  t = round(tmp$statistic, 2),
  `p-value` = signif(tmp$p.value,
                     3)
)[order(tmp$p.value)[1:10],]





vpsig <- viperSignature(myExpressionSet, "group_cell", "Granulosa", verbose = FALSE)


vpres <- viper(vpsig, regul, verbose = FALSE) 






pos <- pData(vpres)[["group_stage_cell"]]  %in% c(pData(vpres)[["group_stage_cell"]] %>% unique())

d1 <- exprs(vpres)[, pos]

colnames(d1) <- pData(vpres)[["group_stage_cell"]][pos]

dd <- dist(t(d1), method = "euclidean")
tmp <- as.matrix(dd)
heatmap(
  as.matrix(dd), 
  Rowv = as.dendrogram(hclust(dd, method = "average")), 
  symm = T)

dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd), method = "average")),
        symm = T)

BiocManager::install("dorothea")
library(dorothea)
library(ggplot2)
library(dplyr)



dorothea_regulon_mouse <-
  get(data("dorothea_mm", package = "dorothea"))
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
scdata <- readRDS('./scdata_tpm.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Oocyte')
scdata@assays$RNA$data <- scdata@assays$RNA$counts
scdata_use <- CreateSeuratObject(
  counts = scdata@assays$RNA$counts,
  meta.data = scdata@meta.data
  )
scdata_use <- NormalizeData(scdata_use)
scdata_use <- run_viper(
  scdata_use,
  regulon,
  options = list(
    method = "scale",
    minsize = 4,
    eset.filter = FALSE,
    cores = 1,
    verbose = TRUE
  ),tidy = FALSE
)

Idents(scdata_use) <- scdata_use$group_stage
viper_scores_df <- GetAssayData(
    scdata_use,
    slot = "data",
    assay = "dorothea"
  ) %>%
  data.frame(check.names = F) %>%
  t()

CellsClusters <- data.frame(
  cell = names(Idents(scdata_use)),
  cell_type = Idents(scdata_use) %>% as.character(),
  check.names = F
)

viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>%
  rownames_to_column("cell") %>%
  gather(tf, activity,-cell) %>%
  inner_join(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  summarise(var = var(avg))  %>%
  ungroup() %>%
  arrange(desc(var)) %>%
  select(tf) %>% 
  dplyr::slice(1:n())

summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) %>% 
  t() %>% as.data.frame() %>% 
  select(c("Secondary","Early antral","Antral","Preovulatory"))
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white", "red"))(palette_length)

my_breaks <- c(
  seq(
    min(summarized_viper_scores_df),
    0,
    length.out = ceiling(palette_length / 2) + 1
  ),
  seq(
    max(summarized_viper_scores_df) / palette_length,
    max(summarized_viper_scores_df),
    length.out = floor(palette_length / 2)
  )
)

viper_hmap <- pheatmap::pheatmap(
  summarized_viper_scores_df,
  scale = 'none',
  fontsize = 14,
  fontsize_row = 10,
  color = my_color,
  breaks = my_breaks,
  main = "DoRothEA (ABC)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA
)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
ha_top = HeatmapAnnotation(
  foo_top = anno_block(
    gp = gpar(fill = c('#C34F73','#0088D4','#4B8600','#5F559B')),
    labels = summarized_viper_scores_df %>% colnames(),
    labels_gp = gpar(col = "white", fontsize = 18)
  )
)


rows_to_label <- c(
  "Zfpm1","Nfkb2","Hnf1b","Rela","Tob1","Mapk14","Arid3a",
  "Sall2","Crtc1","Bhlhe40","Pias1","Pparg","Zbtb32")
rows_to_label <- mrs_select
p <- Heatmap(
  summarized_viper_scores_df %>% as.matrix(),
  name = "mat",
  col = colorRampPalette(c("blue","white", "#A90C38FF"))(100),
  cluster_rows = TRUE,
  row_split = 3,
  column_split = c(1,2,3,4),
  row_gap = unit(5,units = 'pt'),
  top_annotation = ha_top,
  column_title = NULL,
  cluster_columns = F,
  show_row_names = TRUE,
  row_labels = ifelse(
    rownames(summarized_viper_scores_df) %in% c(rows_to_label),
    rownames(summarized_viper_scores_df), ""),
  row_names_gp = gpar(fontsize = 12),
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 30),
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(
    title = 'TF activity',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
    title_position = "leftcenter-rot",
    legend_direction = 'vertical',
    legend_position = "topright",
    legend_height = unit(30,units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
p
row_order_list <-
  lapply(1:length(row_order(p)), function(i) {
    rownames(summarized_viper_scores_df)[row_order(p)[[i]]]
  })
names(row_order_list) <- paste('Cluster ',1:3,sep = '')

df_hiecluster_res <- data.frame(
  Cluster = rep(
    names(row_order_list), 
    times = sapply(row_order_list, length)
    ),
  Gene = unlist(row_order_list)
)
write.csv(df_hiecluster_res,'./df_hiecluster_dorothea.csv')
tf_select <- row_order_list[[1]]
VlnPlot(scdata_use,features = tf_select)
data_plot <-
  FetchData(scdata_use, vars = c(tf_select, 'group_stage')) %>%
  pivot_longer(cols = names(.)[-ncol(.)],
               names_to = 'TF',
               values_to = 'tpm')
colnames(data_plot)
ggplot(data_plot, aes(x = group_stage, y = tpm)) +
  geom_violin(
    trim = F,
    size = 0.2,
    show.legend = FALSE,
    width = 0.9
  ) + labs(y = NULL, x = NULL) +
  geom_boxplot(aes(fill = group_stage),width = 0.1,show.legend = TRUE) +
  facet_wrap(~TF,nrow = 5,scales = 'free',switch = 'y') +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(
    plot.margin = margin(l = 5,r = 5,b = 10,t = 5),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_line(size = 0.2, color = "black"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length = unit(.5, "mm"),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.text = element_text(size = 18)
  )

library(viper)
library(Seurat)
library(Biobase)
library(dplyr)
setwd('./outputFolder_Oocyte_1108/')
load('./alltf_Oocyte.RData')

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Oocyte')



scdata$group_use <- ifelse(
  scdata$group_stage == 'Preovulatory',
  yes = 'Preovulatory',no = 'others'
)
scdata_sce <- as.SingleCellExperiment(scdata)
exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
meta <- scdata@meta.data
metadata <- data.frame(
  labelDescription = colnames(meta),
  row.names = colnames(meta)
)
phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                 phenoData=phenoData,
                                 annotation="scdata")
myExpressionSet

regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

print(regul)


scdata$group_use %>% unique()


signature <-
  rowTtest(myExpressionSet, "group_use", 'Preovulatory', "others")

head(signature[[2]])[1:6,]

signature <-
  (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
     sign(signature$statistic))[,1]
head(signature)[1:6]
signature_sub <- signature[!is.na(signature)]




nullmodel <-
  ttestNull(
    x = myExpressionSet,
    pheno = "group_use", 
    group1 = 'Preovulatory',
    group2 = "others",
    per = 1000,
    repos = TRUE,
    verbose = FALSE
  )

head(nullmodel)[, 1:6]


regul

mrs <-
  msviper(
    ges = signature_sub,
    regulon =  regul,
    nullmodel =  nullmodel,
    verbose = T
  )

tmp <- summary(mrs,mrs = mrs$signature %>% length())
view(tmp)

plot(
  mrs,
  cex = 2,
  mrs = 20,
  density = 1,
  include = c("expression", "activity"),
  gama = 2
)
mrs_select <- c('Gnrh1','Hmgb2','Sub1')
tmp <- summary(mrs,mrs = mrs$signature %>% length())
tmp %>% 
  filter(Regulon %in% mrs_select)
tmp$Regulon[grepl('Hmg',tmp$Regulon)]
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)
pdf(
  './TF_trans_Oocyte.pdf',
  width = 8,height = 2
)
par(mai = c(0.5,1,0.5,1))
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)
dev.off()





mrs_ledge <- ledge(mrs)

summary(mrs_ledge)
df_res <- summary(mrs_ledge,mrs = mrs_ledge$signature %>% length()) %>% 
  select(-c('Ledge')) %>% 
  filter(!is.na(Regulon))
df_res$Ledge <- lapply(df_res$Regulon,FUN = function(regulon){
  mrs_ledge$ledge[[regulon]] %>% paste(collapse = ';')
}) %>% unlist()
view(df_res)
write.csv(df_res,'./df_res_Oocyte_Preo.csv')






vpres <- viper(myExpressionSet, regul, verbose = FALSE)

dim(vpres)

tmp <- rowTtest(vpres,  "group_use", 'Preovulatory', "others")

data.frame(
  Gene = rownames(tmp$p.value),
  t = round(tmp$statistic, 2),
  `p-value` = signif(tmp$p.value,
                     3)
)[order(tmp$p.value)[1:10],]





vpsig <- viperSignature(myExpressionSet, "group_use", "Preovulatory", verbose = FALSE)


vpres <- viper(vpsig, regul, verbose = FALSE) 






pos <- pData(vpres)[["group_stage_cell"]]  %in% c(pData(vpres)[["group_stage_cell"]] %>% unique())

d1 <- exprs(vpres)[, pos]

colnames(d1) <- pData(vpres)[["group_stage_cell"]][pos]

dd <- dist(t(d1), method = "euclidean")
tmp <- as.matrix(dd)
heatmap(
  as.matrix(dd), 
  Rowv = as.dendrogram(hclust(dd, method = "average")), 
  symm = T)

dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd), method = "average")),
        symm = T)


save.image('./alltf_Oocyte.RData')

library(viper)
library(Seurat)
library(Biobase)
library(dplyr)
setwd('./outputFolder_Granulosa_1108/')
load('./alltf_GC.RData')

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Granulosa')



scdata$group_use <- ifelse(
  scdata$group_stage == 'Preovulatory',
  yes = 'Preovulatory',no = 'others'
)
scdata_sce <- as.SingleCellExperiment(scdata)
exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
meta <- scdata@meta.data
metadata <- data.frame(
  labelDescription = colnames(meta),
  row.names = colnames(meta)
)
phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                 phenoData=phenoData,
                                 annotation="scdata")
myExpressionSet

regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

print(regul)


scdata$group_use %>% unique()


signature <-
  rowTtest(myExpressionSet, "group_use", 'Preovulatory', "others")

head(signature[[2]])[1:6,]

signature <-
  (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
     sign(signature$statistic))[,1]
head(signature)[1:6]
signature_sub <- signature[!is.na(signature)]




nullmodel <-
  ttestNull(
    x = myExpressionSet,
    pheno = "group_use", 
    group1 = 'Preovulatory',
    group2 = "others",
    per = 1000,
    repos = TRUE,
    verbose = FALSE
  )

head(nullmodel)[, 1:6]


regul

mrs <-
  msviper(
    ges = signature_sub,
    regulon =  regul,
    nullmodel =  nullmodel,
    verbose = T
  )

tmp <- summary(mrs,mrs = mrs$signature %>% length())
view(tmp)

viper:::plot.msviper
plot(
  mrs,
  cex = 2,
  mrs = 20,
  density = 1,
  include = c("expression", "activity"),
  gama = 2
)
mrs_select <- c('Foxl2','Foxo1','Creb3l2')
tmp <- summary(mrs,mrs = mrs$signature %>% length())
tmp %>% 
  filter(Regulon %in% mrs_select)
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)
mrs_select <- c('Cebpb','Foxo1','Stat3')
tmp <- summary(mrs,mrs = mrs$signature %>% length())
tmp %>% 
  filter(Regulon %in% mrs_select)
tmp$Regulon[grepl('Hmg',tmp$Regulon)]
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)
pdf(
  './TF_trans_Granulosa.pdf',
  width = 8,height = 2
)
par(mai = c(0.5,1,0.5,1))
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)
dev.off()





mrs_ledge <- ledge(mrs)

summary(mrs_ledge)
df_res <- summary(mrs_ledge,mrs = mrs_ledge$signature %>% length()) %>% 
  select(-c('Ledge')) %>% 
  filter(!is.na(Regulon))
df_res$Ledge <- lapply(df_res$Regulon,FUN = function(regulon){
  mrs_ledge$ledge[[regulon]] %>% paste(collapse = ';')
}) %>% unlist()
view(df_res)
write.csv(df_res,'./df_res_Oocyte_Preo.csv')






vpres <- viper(myExpressionSet, regul, verbose = FALSE)

dim(vpres)

tmp <- rowTtest(vpres,  "group_use", 'Preovulatory', "others")

data.frame(
  Gene = rownames(tmp$p.value),
  t = round(tmp$statistic, 2),
  `p-value` = signif(tmp$p.value,
                     3)
)[order(tmp$p.value)[1:10],]





vpsig <- viperSignature(myExpressionSet, "group_use", "Preovulatory", verbose = FALSE)


vpres <- viper(vpsig, regul, verbose = FALSE) 






pos <- pData(vpres)[["group_stage_cell"]]  %in% c(pData(vpres)[["group_stage_cell"]] %>% unique())

d1 <- exprs(vpres)[, pos]

colnames(d1) <- pData(vpres)[["group_stage_cell"]][pos]

dd <- dist(t(d1), method = "euclidean")
tmp <- as.matrix(dd)
heatmap(
  as.matrix(dd), 
  Rowv = as.dendrogram(hclust(dd, method = "average")), 
  symm = T)

dd <- viperSimilarity(d1)
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd), method = "average")),
        symm = T)



save.image('./alltf_GC.RData')

library(Seurat)
library(Biobase)
library(dplyr)
library(viper)
setwd('./data_tf/')

scdata <- readRDS('./scdata_tpm.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Oocyte')

group_stage_use <- 'Secondary'
print(group_stage_use)
scdata$group_use <- ifelse(
  scdata$group_stage == group_stage_use,
  yes = group_stage_use,no = 'others'
)
scdata_sce <- as.SingleCellExperiment(scdata)
exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
meta <- scdata@meta.data
metadata <- data.frame(
  labelDescription = colnames(meta),
  row.names = colnames(meta)
)
phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                 phenoData=phenoData,
                                 annotation="scdata")

regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

print(regul)

scdata$group_use %>% unique()


signature <-
  rowTtest(myExpressionSet, "group_use", group_stage_use, "others")


signature <-
  (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
     sign(signature$statistic))[,1]
signature_sub <- signature[!is.na(signature)]




nullmodel <-
  ttestNull(
    x = myExpressionSet,
    pheno = "group_use", 
    group1 = group_stage_use,
    group2 = "others",
    per = 1000,
    repos = TRUE,
    verbose = FALSE
  )




mrs <-
  msviper(
    ges = signature_sub,
    regulon =  regul,
    nullmodel =  nullmodel,
    verbose = T
  )
res_mrs <- mrs

tmp <- summary(mrs,mrs = mrs$signature %>% length())
View(tmp)


par(mai = c(0,0,0.5,1))
tmp <- data.frame(
  nes = mrs$es$nes,
  p.value = mrs$es$p.value,
  activity = (maobject$es$nes/max(abs(maobject$es$nes)))
)
tmp$activity %>% unique()
mrs_select <- tmp %>% 
  tibble::rownames_to_column('mrs') %>% 
  mutate(
    group_activity = ifelse(activity>0,'1','-1') %>% as.numeric()
  ) %>% 
  group_by(group_activity) %>%
  arrange(p.value) %>% 
  filter((activity == 1 & row_number() <= 10) | (activity == -1 & row_number() <= 3)) %>%
  arrange(desc(nes)) %>% 
  pull(mrs)

mrs_select <- c("Zfp263","Zbtb7a","Mxi1","Mbd1","Mnt")
plot(
  x = mrs,
  cex = 1,
  mrs = mrs_select,
  density = 0,
  include = c("expression", "activity"),
  gama = 0
)





mrs <- ledge(mrs)


df_res <- summary(mrs,mrs = mrs$signature %>% length()) %>% 
  select(-c('Ledge'))
df_res$Ledge <- lapply(df_res$Regulon,FUN = function(regulon){
  mrs$ledge[[regulon]] %>% paste(collapse = ';')
}) %>% unlist()
view(df_res)
file_name <- paste(
  './res',
  '/',
  group_stage_use %>% str_replace_all(' ', '_'),
  '/df_res.csv',
  sep = ''
)
write.csv(df_res,file = file_name)

library(viper)
library(Seurat)
library(Biobase)
library(dplyr)
library(stringr)
library(tibble)
setwd('./outputFolder_Granulosa_1108/')



scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Granulosa')
path_dir <- './res_Granulosa_1110/'
dir.create(path_dir,showWarnings = FALSE)


for(group_stage_use in  (scdata$group_stage %>% unique())){
  print(group_stage_use)
  dir.create(paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    sep = ''
  ))
  scdata$group_use <- ifelse(
    scdata$group_stage == group_stage_use,
    yes = group_stage_use,no = 'others'
  )
  scdata_sce <- as.SingleCellExperiment(scdata)
  exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
  meta <- scdata@meta.data
  metadata <- data.frame(
    labelDescription = colnames(meta),
    row.names = colnames(meta)
  )
  phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
  myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                   phenoData=phenoData,
                                   annotation="scdata")

  regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

  print(regul)

  scdata$group_use %>% unique()


  signature <-
    rowTtest(myExpressionSet, "group_use", group_stage_use, "others")


  signature <-
    (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
       sign(signature$statistic))[,1]
  signature_sub <- signature[!is.na(signature)]




  nullmodel <-
    ttestNull(
      x = myExpressionSet,
      pheno = "group_use", 
      group1 = group_stage_use,
      group2 = "others",
      per = 1000,
      repos = TRUE,
      verbose = FALSE
    )


  

  mrs <-
    msviper(
      ges = signature_sub,
      regulon =  regul,
      nullmodel =  nullmodel,
      verbose = T
    )

  tmp <- summary(mrs,mrs = mrs$signature %>% length())
  view(tmp)

  file_name <- paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    '/heatmap_top30.pdf',
    sep = ''
  )
  pdf(file = file_name,width = 16,height = 25)
  par(mai = c(0,0,0.5,1))
  plot(
    mrs,
    cex = 5,
    mrs = 100,
    density = 1,
    include = c("expression", "activity"),
    gama = 2
  )
  dev.off()





  mrs <- ledge(mrs)


  df_res <- summary(mrs,mrs = mrs$signature %>% length()) %>% 
    select(-c('Ledge'))
  df_res$Ledge <- lapply(df_res$Regulon,FUN = function(regulon){
    mrs$ledge[[regulon]] %>% paste(collapse = ';')
  }) %>% unlist()
  view(df_res)
  file_name <- paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    '/df_res.csv',
    sep = ''
  )
  write.csv(df_res,file = file_name)
}




library(ggplot2)
library(tibble)
library(pheatmap)
library(Seurat)
library(stringr)
library(purrr)
library(reshape2)
library(tidyr)
path_dir <- './res_Granulosa_1110/'
setwd(path_dir)
file_names <- c("Secondary","Early_antral","Antral","Preovulatory")
data_list <- lapply(file_names, function(file_name){
  print(file_name)
  data_tmp <-
    read.csv(
      paste('./', file_name, '/', 'df_res.csv', sep = ''), 
      row.names = 1) %>% 
    mutate(group = file_name) %>% 
    filter(!is.na(Regulon))
}) %>% 
  do.call(rbind,.) %>% 
  mutate(group = factor(group,levels = file_names))

colnames(data_list)

data_plot <- data_list %>% 
  select(Regulon,NES,group) %>% 
  pivot_wider(names_from = group,values_from = NES) %>% 
  column_to_rownames('Regulon')
pheatmap(
  mat = data_plot,
  scale = 'row',
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  cutree_rows = 6,
  fontsize_col = 20,
  gaps_col = c(1,2,3)
)


standardize <- function(x) {
  (x - mean(x)) / sd(x)
}
data_standardized <- data_list %>% 
  select(Regulon,NES,group) %>% 
  pivot_wider(names_from = group,values_from = NES) %>% 
  column_to_rownames('Regulon')
data_standardized <- apply(data_standardized, 1, standardize) %>% 
  t() %>% as.data.frame()
file_path <- paste(path_dir,'data_standardized.csv',sep = '')
write.csv(data_standardized,file_path)

library(Seurat)
library(Biobase)
library(dplyr)
library(viper)
setwd('./outputFolder_Oocyte_1108/')

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Oocyte')
path_dir <- './res_Oocyte_1110/'
dir.create(path_dir,showWarnings = FALSE)

for(group_stage_use in  (scdata$group_stage %>% unique())){
  print(group_stage_use)
  dir.create(paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    sep = ''
  ))
  scdata$group_use <- ifelse(
    scdata$group_stage == group_stage_use,
    yes = group_stage_use,no = 'others'
  )
  scdata_sce <- as.SingleCellExperiment(scdata)
  exprMatrix <- scdata@assays$RNA$counts %>% as.matrix()
  meta <- scdata@meta.data
  metadata <- data.frame(
    labelDescription = colnames(meta),
    row.names = colnames(meta)
  )
  phenoData <- new("AnnotatedDataFrame",data=meta,varMetadata=metadata)
  myExpressionSet <- ExpressionSet(assayData=exprMatrix,
                                   phenoData=phenoData,
                                   annotation="scdata")

  regul <- aracne2regulon('./network.txt', myExpressionSet, verbose = FALSE)

  print(regul)

  scdata$group_use %>% unique()


  signature <-
    rowTtest(myExpressionSet, "group_use", group_stage_use, "others")


  signature <-
    (qnorm(signature$p.value / 2, lower.tail = FALSE) * 
       sign(signature$statistic))[,1]
  signature_sub <- signature[!is.na(signature)]




  nullmodel <-
    ttestNull(
      x = myExpressionSet,
      pheno = "group_use", 
      group1 = group_stage_use,
      group2 = "others",
      per = 1000,
      repos = TRUE,
      verbose = FALSE
    )


  

  mrs <-
    msviper(
      ges = signature_sub,
      regulon =  regul,
      nullmodel =  nullmodel,
      verbose = T
    )

  tmp <- summary(mrs,mrs = mrs$signature %>% length())
  view(tmp)

  file_name <- paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    '/heatmap_top30.pdf',
    sep = ''
  )
  pdf(file = file_name,width = 16,height = 25)
  par(mai = c(0,0,0.5,1))
  plot(
    mrs,
    cex = 5,
    mrs = 100,
    density = 1,
    include = c("expression", "activity"),
    gama = 2
  )
  dev.off()





  mrs <- ledge(mrs)


  df_res <- summary(mrs,mrs = mrs$signature %>% length()) %>% 
    select(-c('Ledge'))
  df_res$Ledge <- lapply(df_res$Regulon,FUN = function(regulon){
    mrs$ledge[[regulon]] %>% paste(collapse = ';')
  }) %>% unlist()
  view(df_res)
  file_name <- paste(
    path_dir,
    group_stage_use %>% str_replace_all(' ', '_'),
    '/df_res.csv',
    sep = ''
  )
  write.csv(df_res,file = file_name)
}


library(ggplot2)
library(tibble)
library(pheatmap)
library(Seurat)
library(stringr)
library(purrr)
library(reshape2)
library(tidyr)
path_dir <- './res_Oocyte_1110/'
setwd(path_dir)
file_names <- c("Secondary","Early_antral","Antral","Preovulatory")
data_list <- lapply(file_names, function(file_name){
  print(file_name)
  data_tmp <-
    read.csv(
      paste('./', file_name, '/', 'df_res.csv', sep = ''), 
      row.names = 1) %>% 
    mutate(group = file_name) %>% 
    filter(!is.na(Regulon))
}) %>% 
  do.call(rbind,.) %>% 
  mutate(group = factor(group,levels = file_names))

colnames(data_list)

data_plot <- data_list %>% 
  select(Regulon,NES,group) %>% 
  pivot_wider(names_from = group,values_from = NES) %>% 
  column_to_rownames('Regulon')
pheatmap(
  mat = data_plot,
  scale = 'row',
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  cutree_rows = 6,
  fontsize_col = 20,
  gaps_col = c(1,2,3)
)


standardize <- function(x) {
  (x - mean(x)) / sd(x)
}
data_standardized <- data_list %>% 
  select(Regulon,NES,group) %>% 
  pivot_wider(names_from = group,values_from = NES) %>% 
  column_to_rownames('Regulon')
data_standardized <- apply(data_standardized, 1, standardize) %>% 
  t() %>% as.data.frame()
file_path <- paste(path_dir,'data_standardized.csv',sep = '')
write.csv(data_standardized,file_path)


maobject <- mrs
mrs <- 10
color = c("cornflowerblue", "salmon") 
pval = NULL
bins = 500
cex = 0
density = 0
smooth = 0
sep = 0.2
hybrid = TRUE
include = c("expression","activity")
gama = 2

plot_self <- function(
    x, mrs = 10, color = c("cornflowerblue", "salmon"), 
    pval = NULL, bins = 500, cex = 0, density = 0, smooth = 0, 
    sep = 0.2, hybrid = TRUE, include = c("expression", "activity"), 
    gama = 2, ...){
  maobject <- x
  rm(x)
  marg <- par("mai")
  rlist <- maobject$signature
  if (ncol(rlist) > 0) 
    rlist <- rowMeans(rlist)
  if (length(mrs) == 1 & is.numeric(mrs[1])) {
    mrs <- names(maobject$es$nes)[order(maobject$es$p.value)[1:round(mrs)]]
    mrs <- mrs[order(maobject$es$nes[match(mrs, names(maobject$es$nes))], 
                     decreasing = TRUE)]
  }else{
    mrs <- mrs
  }
  groups <- maobject$regulon[match(mrs, names(maobject$regulon))]
  if (is.null(pval)) 
    pval <- maobject$es$p.value[match(mrs, names(maobject$es$nes))]
  if (is.data.frame(pval)) 
    pval <- as.matrix(pval)
  if (is.matrix(rlist)) 
    rlist <- rlist[, 1]
  if (min(rlist) < 0){
    rlist <- sort(rlist)
  }else{
    rlist <- sort(-rlist)
  }
  groups <- groups[length(groups):1]
  color1 <- color
  layout(
    matrix(1:2, 1, 2), 
    widths = c(10 - length(include),length(include))
  )
  color <- rgb2hsv(col2rgb(color))
  satval <- color[3, ]
  color <- color[1, ]
  preg <- as.numeric(any(sapply(groups, function(x) any(x$tfmode < 0)))) + 1
  textsize <- 1
  xlimit <- c(0, length(rlist) * (1 + 0.2 * max(nchar(names(groups)))/8))
  if (!is.null(pval)) {
    if (is.matrix(pval)) {
      pval <- pval[nrow(pval):1, ]
      pval <- pval[, ncol(pval):1]
      xlimit <- c(-length(rlist) * (sep * ncol(pval) + 
                                      0.02), length(rlist) * (1 + 0.2 * max(nchar(names(groups)))/9))
    }else {
      pval <- pval[length(pval):1]
      xlimit <- c(-length(rlist) * 0.12, length(rlist) * 
                    (1 + 0.2 * max(nchar(names(groups)))/9))
      textsize = 0.8
    }
  }
  if (length(include) > 0 & include[1] != "") 
    par(mai = c(marg[1:3], 0.05))
  plot(1, type = "n", ylim = c(0, length(groups)), xlim = xlimit, 
       axes = FALSE, ylab = "", xlab = "", yaxs = "i")
  if (cex > 0){
    textsize <- (length(groups) <= 20) * cex + 
      (length(groups) > 20) * (20/length(groups)) * cex
    
  } 
  switch(preg, {
    for (i in 1:length(groups)) {
      densi <- rep(0, length(rlist))
      x <- which(names(rlist) %in% names(groups[[i]]$tfmode))
      if (length(x) > 0) {
        densi[x] <- 1
        denStep <- round(length(densi)/bins)
        x1 <- x[x < denStep]
        x2 <- x[x >= denStep & x <= (length(rlist) - 
                                       denStep)]
        x3 <- x[x > (length(rlist) - denStep)]
        densiRes <- sapply(x2, function(i, densi, denStep) sum(densi[(i - 
                                                                        denStep):(i + denStep)]), densi = densi, denStep = denStep)
        densiRes <- densiRes/max(densiRes)
        temp <- rlist[x]
        if (satval[1] == 0) temp <- hsv((temp < 0) * 
                                          color[1] + (temp > 0) * color[2], satval, 
                                        1 - densiRes) else temp <- hsv((sign(temp) < 
                                                                          0) * color[1] + (sign(temp) > 0) * color[2], 
                                                                       densiRes, satval)
        for (ii in order(densiRes)) lines(c(x[ii], x[ii]), 
                                          c(i - 1, i), col = temp[ii])
        if (density > 0) {
          denStep <- round(length(densi)/density)
          xpos <- seq(denStep, length(rlist) - denStep, 
                      length = density)
          densiRes <- sapply(xpos, function(i, densi, 
                                            denStep) {
            sum(densi[(i - denStep):(i + denStep)])
          }, densi = densi, denStep = denStep)
          densiRes <- densiRes/max(densiRes)
          if (smooth > 0) densiRes <- smooth.spline(xpos, 
                                                    densiRes, spar = smooth)$y
          lines(xpos, i + densiRes - 1)
        }
      }
    }
    text(rep(length(rlist) * 1.02, length(groups)), 1:length(groups) - 
           0.5, names(groups), adj = 0, cex = textsize)
    if (!is.null(pval)) {
      if (is.matrix(pval)) {
        for (i in 1:ncol(pval)) {
          text(rep(-length(rlist) * (sep * (i - 1) + 
                                       0.02), length(groups)), 1:length(groups) - 
                 0.5, pval[, i], adj = 1, cex = 0.85 * textsize)
          text(-length(rlist) * (sep * (i - 1) + 0.02), 
               length(groups) + 0.5, colnames(pval)[i], 
               adj = 1, cex = 1)
        }
      } else {
        text(rep(-length(rlist) * 0.02, length(groups)), 
             1:length(groups) - 0.5, signif(pval, 3), adj = 1, 
             cex = 0.85 * textsize)
        text(0, length(groups) + 0.5, ifelse(max(pval) > 
                                               1, "oddsR", "p-value"), adj = 1, cex = 1.2)
      }
    }
    text(length(rlist) * 1.02, length(groups) + 0.5, "Set", 
         adj = 0, cex = 1.2)
  }, {
    for (i in 1:length(groups)) {
      for (ii in 1:2) {
        tset <- groups[[i]]$tfmode
        tset <- tset[(tset < 0 & ii == 1) | (!(tset < 0) & ii == 2)]
        tset1 <- names(tset)[abs(tset) > 0.5]
        tset2 <- names(tset)[abs(tset) < 0.5]
        if (length(tset) > 1) {
          densi <- rep(0, length(rlist))
          x <- match(names(tset), names(rlist))
          tw1 <- rep(1, length(x))
          if (hybrid) {
            x <- match(tset1, names(rlist))
            if (ii == 1) {
              x <- c(x, match(tset2, names(sort(-abs(rlist) * 
                                                  sign(maobject$es$nes[names(maobject$es$nes) == 
                                                                         names(groups)[i]])))))
            } else {
              x <- c(x, match(tset2, names(sort(abs(rlist) * 
                                                  sign(maobject$es$nes[names(maobject$es$nes) == 
                                                                         names(groups)[i]])))))
            }
            tw1 <- groups[[i]]$likelihood[match(c(tset1, 
                                                  tset2), names(groups[[i]]$tfmode))]
            tw1 <- tw1/max(tw1)
          }
          densi[x] <- 1
          denStep <- round(length(densi)/bins)
          x1 <- x[x < denStep]
          x2 <- x[x >= denStep & x <= (length(rlist) - 
                                         denStep)]
          x3 <- x[x > (length(rlist) - denStep)]
          densiRes <- sapply(x2, function(i, densi, 
                                          denStep) {
            sum(densi[(i - denStep):(i + denStep)])
          }, densi = densi, denStep = denStep)
          densiRes <- densiRes * (tw1[x >= denStep & 
                                        x <= (length(rlist) - denStep)])
          densiRes <- densiRes/max(densiRes)
          temp <- rlist[x]
          temp <- hsv(color[ii], densiRes, satval)
          for (iii in order(densiRes)) {
            lines(c(x[iii], x[iii]), c(i - 1 + (ii - 
                                                  1)/2, i - 1 + ii/2), col = temp[iii])
          }
          if (density > 0) {
            denStep <- round(length(densi)/density)
            xpos <- seq(denStep, length(rlist) - denStep, 
                        length = density)
            densiRes <- sapply(xpos, function(i, densi, 
                                              denStep) sum(densi[(i - denStep):(i + 
                                                                                  denStep)]), densi = densi, denStep = denStep)
            densiRes <- densiRes/max(densiRes)
            if (smooth > 0) densiRes <- smooth.spline(xpos, 
                                                      densiRes, spar = smooth)$y
            lines(xpos, i - 1 + densiRes/2 + (ii - 1)/2)
          }
        }
      }
    }
    text(rep(length(rlist) * 1.02, length(groups)), 1:length(groups) - 
           0.5, names(groups), adj = 0, cex = textsize)
    if (!is.null(pval)) {
      if (is.matrix(pval)) {
        for (i in 1:ncol(pval)) {
          text(rep(-length(rlist) * (sep * (i - 1) + 
                                       0.02), length(groups)), 1:length(groups) - 
                 0.5, pval[, i], adj = 1, cex = 0.85 * textsize)
          text(-length(rlist) * (sep * (i - 1) + 0.02), 
               length(groups) + 0.5, colnames(pval)[i], 
               adj = 1, cex = 1)
        }
      } else {
        text(rep(-length(rlist) * 0.02, length(groups)), 
             1:length(groups) - 0.5, signif(pval, 3), adj = 1, 
             cex = 0.85 * textsize)
        axis(3, -0.05 * length(rlist), ifelse(max(pval) > 
                                                1, "oddsR", "p-value"), adj = 1, cex = 1.2, 
             tick = FALSE, line = -0.5)
      }
    }
    axis(3, length(rlist) * 1.05, "Set", adj = 0, cex = 1.2, 
         line = -0.5, tick = FALSE)
  })
  abline(h = 0:length(groups))
  lines(c(0, 0), c(0, length(groups)))
  lines(c(length(rlist), length(rlist)), c(0, length(groups)))
  ss <- maobject$signature
  if (!is.null(dim(ss))) 
    ss <- ss[, 1]
  x <- NULL
  xn <- NULL
  xpos <- NULL
  if ("activity" %in% include) {
    x <- maobject$es$nes[match(mrs, names(maobject$es$nes))]/max(abs(maobject$es$nes))
    xn <- "Act"
  }
  if ("expression" %in% include) {
    x <- cbind(x, ss[match(mrs, names(ss))]/max(abs(ss)))
    xn <- c(xn, "Exp")
    xpos <- rank(-abs(ss))[match(mrs, names(ss))]
  }
  if (!is.null(x)) {
    if (is.null(dim(x))) 
      dim(x) <- c(length(x), 1)
    rownames(x) <- colnames(x) <- NULL
    marg[2] <- 0.05
    par(mai = marg)
    scmax <- max(abs(x), na.rm = TRUE)
    x <- abs(x/scmax)^gama * sign(x)
    x <- filterRowMatrix(x, nrow(x):1)
    x1 <- x
    x1[is.na(x1)] <- 0
    coli <- hsv(ifelse(x1 < 0, color[1], color[2]), abs(x1), 
                1)
    coli[is.na(x)] <- hsv(0, 0, 0.5)
    image(1:ncol(x), 1:nrow(x), t(matrix(1:(ncol(x) * nrow(x)), 
                                         nrow(x), ncol(x))), col = coli, ylab = "", xlab = "", 
          axes = FALSE, yaxs = "i")
    box()
    grid(ncol(x), nrow(x), col = "black", lty = 1)
    axis(3, 1:length(xn), xn, las = 2, line = -0.5, tick = FALSE)
    axis(4, length(xpos):1, xpos, las = 1, tick = FALSE, 
         line = -0.4, cex.axis = 0.85 * textsize)
  }
  par(mai = marg)
}



scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte')
data_tran <- scdata_tran_sub@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column('gene')
dir.create('./ARACNe_input',showWarnings = FALSE)
write.table(
  data_tran,
  './Oocyte.txt',
  row.names = FALSE,quote = FALSE,sep = '\t')
scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Granulosa')
data_tran <- scdata_tran_sub@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column('gene')
write.table(
  data_tran,
  './Granulosa.txt',
  row.names = FALSE,quote = FALSE,sep = '\t')



load('./mouse_regulators.rda')
trrust_tf <- read.csv(
  './tfs.txt',
  col.names = 'TF') %>% 
  mutate(Group = 'trrust')
trrust_tf$TF %>% unique() %>% length()
data_PISCES <- mouse_regulators %>% stack() %>% rename_all(~c('TF','Group'))
data_PISCES$TF
trrust_tf$TF[!(trrust_tf$TF %in% data_PISCES$TF)]
data_tf_all <- rbind(
  data_PISCES,trrust_tf
) %>% 
  distinct()




write.table(
  data_tf_all$TF,
  './tfs_all.txt',
  row.names = FALSE,quote = FALSE,col.names = FALSE)


scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
save_path <- './ARACNe_input_1128/'
dir.create(save_path,showWarnings = FALSE)
for(group_stage_cell_use in (scdata_tran$group_stage_cell %>% unique())){
  scdata_tran_sub <- subset(scdata_tran,subset = group_stage_cell == group_stage_cell_use)
  data_tran <- scdata_tran_sub@assays$RNA$counts %>% as.data.frame() %>% 
    rownames_to_column('gene')
  file_name <- paste(
    save_path,group_stage_cell_use,'.txt',sep = ''
  )
  write.table(
    x = data_tran,
    file = file_name,
    row.names = FALSE,quote = FALSE,sep = '\t')
}


data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
res_tf <- data.frame(
  tf_name = c(
    data_raw_1$Symbol,data_raw_2$tf
  ) %>% unique()
)




write.table(
  res_tf$tf_name,
  './tfs_all.txt',
  row.names = FALSE,quote = FALSE,col.names = FALSE)
