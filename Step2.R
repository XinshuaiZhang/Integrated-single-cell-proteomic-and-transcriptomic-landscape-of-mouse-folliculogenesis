my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
set.seed(1234)
scdata <- readRDS('./scdata_use_filter_outliers.rds')
scdata <- Seurat::NormalizeData(scdata)
scdata <- ScaleData(object = scdata,features = rownames(scdata))
scdata <- FindVariableFeatures(scdata,nfeatures = 2000)
scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata),npcs = 50)
ElbowPlot(scdata, ndims = 40)
ndim = 30
{
  scdata <- FindNeighbors(scdata, dims = 1:ndim,k.param = 5)
  scdata <- FindClusters(scdata, resolution = 0.5)
  table(Idents(scdata))
}


scdata <- RunTSNE(scdata, dims.use = 1:ndim,dims = 2,perplexity  = 5)
Idents(scdata) <- scdata$group_cell
PCAPlot(scdata, pt.size = 1, label = T)$data %>% 
  rownames_to_column('sample') %>% 
  ggplot(.,mapping = aes(
    x = PC_1,
    y = PC_2,
    color = ident,
    label = sample
  )) +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  guides(color = guide_legend(
    title = 'Cluster',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )
TSNEPlot(scdata, pt.size = 1, label = T)$data %>%
  rownames_to_column('sample') %>%
  ggplot(.,mapping = aes(
    x = tSNE_1,
    y = tSNE_2,
    color = ident,
    label = sample
  )) +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  guides(color = guide_legend(
    title = 'Cluster',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  labs(title = 't-SNE') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )

scdata <- RunUMAP(object = scdata, dims = 1:30,n.neighbors = 28,seed.use = 1234)#28
scdata <- RunUMAP(object = scdata, dims = 1:30,n.neighbors = 28,seed.use = 1234)#28
DimPlot(object = scdata, pt.size=1, reduction = 'umap',group.by = 'group_stage_cell')
data_plot <- DimPlot(object = scdata, pt.size=1, reduction = 'umap',group.by = 'group_stage_cell')$data %>% 
  rownames_to_column('sample') %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    "Secondary_Granulosa","Early antral_Granulosa","Antral_Granulosa","Preovulatory_Granulosa",
    "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte","Preovulatory_Oocyte"
  ) %>% str_replace('_',': ')))
color_group_stage_cell <- c(
  '#B9DDF1FF','#7EAED3FF','#5081AEFF','#2A5783FF',
  "#FFBEB2FF","#F8826BFF","#E33E43FF","#AE123AFF"
)
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage_cell,
         label = sample
       )) +
  geom_point(size = 4) +
  scale_color_manual(values = color_group_stage_cell) +#my36colors[c(1,3)]
  guides(color = guide_legend(
    title = 'Cluster',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  labs(x = 'Umap 1',y = 'Umap 2') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 20),
    axis.text = element_blank(),
    legend.text = element_text(size = 18)
  )
ggsave('./umap_protein.png',width = 9,height = 6)
ggsave('./umap_protein.pdf',width = 9,height = 6)





data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
scdata_pro <- readRDS('./scdata_use_filter_outliers.rds')

meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')
data_singlestage_raw <- FetchData(
  scdata_pro,layer = 'counts',
  vars = rownames(scdata_pro)
) %>% 
  t() %>% as.data.frame()
data_singlestage <-
  data_singlestage_raw %>%
  rownames_to_column('Accessions') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  dplyr::select(-c('Accessions')) %>%
  as.data.frame() %>%
  group_by(Genes) %>%
  summarise(across(.cols = everything(), .fns = ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup() %>% 
  column_to_rownames('Genes') %>% 
  {log10(.+1)}
scdata <- CreateSeuratObject(
  counts = data_singlestage,
  data = data_singlestage,
  meta.data = meta_data_pro,
  min.cells = 2, 
  min.features = 2000)
gNames <- c("bmp15","gdf9","kit","dppa3","zp3","ddx4","gja4",
            "zp1","zp2","zar1","dazl","Rnf24",
            "Hsf1","Ncoa2","Slc39a10",
            "cnmd","Amh","Nr5a2","gja1","Hmgcs2","Foxl2",
            "Runx1","Amhr2") %>% str_to_title()
gNames <- c(
  "Bmp15","Gja4","Zp1","Zp2","Zar1",
  "Gja1","Amh","Hmgcs2","Foxl2",'Amhr2'
)
gNames <- c(
  "Bmp15","Gdf9","Kit","Dppa3",
  "Ncoa2","Cnmd","Runx1"
)
DotPlot(object = scdata,features = gNames,scale = FALSE,group.by = 'group_cell') +
  scale_color_gradient(
    name = 'Mean Exp',
    low = 'grey90',
    high = '#3131F2'
  ) +
  scale_size_continuous(name = 'Fraction') +
  guides(
    size = guide_legend(ncol = 1,title.position = 'left'),
    color = guide_colourbar(nrow = 1,title.position = 'left')
  ) +
  labs(x = 'Features', y= 'Cell Type') +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.direction = 'vertical',
    legend.title = element_text(angle = 90,hjust = 0.5,vjust = 0),
    axis.title.x = element_blank()
  )
ggsave(
  './dotplot_protein.png',
  width = 10,height = 4
  )
ggsave(
  './dotplot_protein.pdf',
  width = 10,height = 4
)



data_pro_clean <- read.csv('./data_singlestage.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
scdata_pro <- readRDS('./scdata_use_filter_outliers.rds')

meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')
data_singlestage_raw <- FetchData(
  scdata_pro,layer = 'counts',
  vars = rownames(scdata_pro)
) %>% 
  t() %>% as.data.frame()
data_singlestage <-
  data_singlestage_raw %>%
  rownames_to_column('Accessions') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  dplyr::select(-c('Accessions')) %>%
  as.data.frame() %>%
  group_by(Genes) %>%
  summarise(across(.cols = everything(), .fns = ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup()
sample_level <- scdata$sample %>% levels()
sample_level_group <- c(
  sample_level[grepl('O',sample_level)],
  sample_level[grepl('GC',sample_level)]
)
gNames <- c("kitl2","Amh","Nr5a2","lnha",#"cnmd",
            "gja1","Fst","Lgr5","Hmgcs2","Foxl2",
            "Runx1","Kitl","Amhr2") %>% str_to_title()
gNames[!(gNames %in% data_singlestage$Genes)]
gNames <- gNames[gNames %in% data_singlestage$Genes]
ProNames <- data_gene_protein$Accessions[data_gene_protein$Genes %in% gNames]
names(gNames) <- ProNames
Idents(scdata_pro) <- scdata_pro$group_cell
VlnPlot(scdata_pro, features = ProNames,cols = my36colors,ncol = 5)
DotPlot(scdata_pro, features = ProNames,group.by = 'group_cell')

data_plot <- FetchData(
  object = scdata_pro,layer = 'counts',
  clean = 'none',vars = c(ProNames,'group_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  as.data.frame() %>%
  mutate(gene = gNames[protein])


ggplot(data = data_plot,mapping = aes(x = group_cell,y = expression,fill = group)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'none',
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

data_plot$group_tmp <- ifelse(
  data_plot$sample %in% c('AO1','AO6','AO8','EO4'),
  yes = 'red',
  no = 'blue'
)
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group_tmp)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 10),
    strip.text = element_text(size = 20)
  )

gNames <- c("bmp15","gdf9","kit","dppa3","zp3","ddx4","gja4",
            "zp1","zp2","h1oo","zar1","Dnmt1o","Mater",
            "Hsf1","Ncoa2","Rnf24","Slc39a10","Tgfbr",
            "dazl","sycp3") %>% str_to_title()
gNames[!(gNames %in% data_singlestage$Genes)]
gNames <- gNames[gNames %in% data_singlestage$Genes]
ProNames <- data_singlestage$Accessions[data_singlestage$Genes %in% gNames]
names(gNames) <- ProNames
Idents(scdata) <- scdata$group_cell
VlnPlot(scdata, features = ProNames,cols = my36colors,ncol = 5)


data_plot <- FetchData(
  object = scdata,layer = 'data',
  clean = 'none',vars = c(ProNames,'group_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  as.data.frame() %>%
  mutate(gene = gNames[protein])


ggplot(data = data_plot,mapping = aes(x = group_cell,y = expression,fill = group)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'none',
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

data_plot$group_tmp <- ifelse(
  data_plot$sample %in% c('AO1','AO6','AO8','EO4'),
  yes = 'red',
  no = 'blue'
)
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group_tmp)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 10),
    strip.text = element_text(size = 20)
  )
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 19),
    strip.text = element_text(size = 20)
  )


FeaturePlot(scdata,features = c('Lamb1','Dag1'),pt.size = 2)
sample_level_group <- c(
  sample_level[seq(1,length(sample_level),2)],
  sample_level[seq(2,length(sample_level),2)]
)
lavel_group_tmp <- c(
  "Granulosa_Secondary", "Granulosa_Early antral","Granulosa_Antral","Granulosa_Preovulatory",
  "Oocyte_Secondary", "Oocyte_Early antral", "Oocyte_Antral","Oocyte_Preovulatory" 
)
gNames <- c("Tgfb2","Tgfb3")
gNames <- c( "Dll3","Jag1","Jag2","Notch2")
gNames <- c('BMP7','BMP15','GDF9','INHBB','INHA') %>% str_to_title()
gNames <- c("Bmpr1a","Bmpr1b","Bmpr2","Acvr1a","Acvr1b","Acvr2a","Acvr2b","Tgfbr1","Tgfbr3")
data_plot <- FetchData(object = scdata,layer = 'count',clean = 'none',vars = c(gNames,'group_cell','group_stage')) %>% 
  rownames_to_column('sample') %>% 
  melt(id.var = c('sample','group_cell','group_stage'),variable.name = 'gene',value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  mutate(group_tmp = paste(group_cell,group_stage,sep = '_')) %>% 
  mutate(group_tmp = factor(group_tmp,levels = lavel_group_tmp)) %>% 
  as.data.frame()
ggplot(data = data_plot,mapping = aes(x = group_tmp,y = expression,fill = group_tmp)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.5),fill = NA) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )




scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata <- Seurat::NormalizeData(scdata)
scdata <- ScaleData(object = scdata,features = rownames(scdata))
scdata <- FindVariableFeatures(scdata,nfeatures = 2000)
scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata),npcs = 50)
ElbowPlot(scdata, ndims = 50)
ndim = 40
{
  scdata <- FindNeighbors(scdata, dims = 1:ndim,k.param = 10)
  scdata <- FindClusters(scdata, resolution = 0.5)
  table(Idents(scdata))
}


scdata <- RunUMAP(object = scdata, dims = 1:30,n.neighbors = 35,seed.use = 1234)
scdata <- RunUMAP(object = scdata, dims = 1:30,n.neighbors = 35,seed.use = 1234)
DimPlot(object = scdata, pt.size=1, reduction = 'umap',group.by = 'group_stage_cell')
data_plot <- DimPlot(object = scdata, pt.size=1, reduction = 'umap',group.by = 'group_stage_cell')$data %>% 
  rownames_to_column('sample') %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    "Secondary_Granulosa","Early antral_Granulosa","Antral_Granulosa","Preovulatory_Granulosa",
    "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte","Preovulatory_Oocyte"
  )))
color_group_stage_cell <- c(
  '#B9DDF1FF','#7EAED3FF','#5081AEFF','#2A5783FF',
  "#FFBEB2FF","#F8826BFF","#E33E43FF","#AE123AFF"
)
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage_cell,
         label = sample
       )) +
  geom_point(size = 4) +
  scale_color_manual(values = color_group_stage_cell) +#my36colors[c(1,3)]
  guides(color = guide_legend(
    title = 'Cluster',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  labs(x = 'Umap 1',y = 'Umap 2') +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 20),
    axis.text = element_blank(),
    legend.text = element_text(size = 18)
  )
ggsave('./umap_transcript.png',width = 9,height = 6)
ggsave('./umap_transcript.pdf',width = 9,height = 6)


gNames <- c("bmp15","gdf9","kit","dppa3","zp3","ddx4","gja4",
            "zp1","zp2","zar1","dazl","Rnf24",
            "Hsf1","Ncoa2","Slc39a10",
            "cnmd","Amh","Nr5a2","gja1","Hmgcs2","Foxl2",
            "Runx1","Amhr2") %>% str_to_title()

gNames <- c(
  "Bmp15","Gja4","Zp1","Zp2",#"Zar1",
  "Gja1","Fst","Hmgcs2"#,"Foxl2",'Amhr2'
)
DotPlot(object = scdata,features = gNames,scale = FALSE,group.by = 'group_cell') +
  scale_color_gradient(
    name = 'Mean Exp',
    low = 'grey90',
    high = '#3131F2'
  ) +
  scale_size_continuous(name = 'Fraction') +
  guides(
    size = guide_legend(ncol = 1,title.position = 'left'),
    color = guide_colourbar(nrow = 1,title.position = 'left')
  ) +
  labs(x = 'Features', y= 'Cell Type') +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.direction = 'vertical',
    legend.title = element_text(angle = 90,hjust = 0.5,vjust = 0),
    axis.title.x = element_blank()
  )
ggsave(
  './dotplot_transcript.png',
  width = 10,height = 4
)
ggsave(
  './dotplot_transcript.pdf',
  width = 10,height = 4
)



sample_level <- scdata$sample %>% levels()
sample_level_group <- c(
  sample_level[grepl('O',sample_level)],
  sample_level[grepl('GC',sample_level)]
)
gNames <- c("kitl2","Amh","Nr5a2","lnha",#"cnmd",
            "gja1","Fst","Lgr5","Hmgcs2","Foxl2",
            "Runx1","Kitl","Amhr2") %>% str_to_title()
gNames[!(gNames %in% data_copy_num_sorted$Genes)]
gNames <- gNames[gNames %in% data_copy_num_sorted$Genes]
ProNames <- data_copy_num_sorted$ProteinGroups[data_copy_num_sorted$Genes %in% gNames]
names(gNames) <- ProNames
Idents(scdata) <- scdata$group_cell
VlnPlot(scdata, features = ProNames,cols = my36colors,ncol = 5)


data_plot <- FetchData(
  object = scdata,layer = 'data',
  clean = 'none',vars = c(ProNames,'group_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  as.data.frame() %>%
  mutate(gene = gNames[protein])


ggplot(data = data_plot,mapping = aes(x = group_cell,y = expression,fill = group)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'none',
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

data_plot$group_tmp <- ifelse(
  data_plot$sample %in% c('AO1','AO6','AO8','EO4'),
  yes = 'red',
  no = 'blue'
)
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group_tmp)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 10),
    strip.text = element_text(size = 20)
  )

gNames <- c("bmp15","gdf9","kit","dppa3","zp3","ddx4","gja4",
            "zp1","zp2","h1oo","zar1","Dnmt1o","Mater",
            "Hsf1","Ncoa2","Rnf24","Slc39a10","Tgfbr",
            "dazl","sycp3") %>% str_to_title()
gNames[!(gNames %in% data_copy_num_sorted$Genes)]
gNames <- gNames[gNames %in% data_copy_num_sorted$Genes]
ProNames <- data_copy_num_sorted$Accessions[data_copy_num_sorted$Genes %in% gNames]
names(gNames) <- ProNames
Idents(scdata) <- scdata$group_cell
VlnPlot(scdata, features = ProNames,cols = my36colors,ncol = 5)


data_plot <- FetchData(
  object = scdata,layer = 'data',
  clean = 'none',vars = c(ProNames,'group_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  as.data.frame() %>%
  mutate(gene = gNames[protein])


ggplot(data = data_plot,mapping = aes(x = group_cell,y = expression,fill = group)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'none',
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

data_plot$group_tmp <- ifelse(
  data_plot$sample %in% c('AO1','AO6','AO8','EO4'),
  yes = 'red',
  no = 'blue'
)
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group_tmp)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 10),
    strip.text = element_text(size = 20)
  )
ggplot(data = data_plot,mapping = aes(x = sample,y = expression,fill = group)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = my36colors[c(1,3)]) +
  facet_wrap(~gene,ncol = 2,scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 19),
    strip.text = element_text(size = 20)
  )


FeaturePlot(scdata,features = c('Lamb1','Dag1'),pt.size = 2)
sample_level_group <- c(
  sample_level[seq(1,length(sample_level),2)],
  sample_level[seq(2,length(sample_level),2)]
)
lavel_group_tmp <- c(
  "Granulosa_Secondary", "Granulosa_Early antral","Granulosa_Antral","Granulosa_Preovulatory",
  "Oocyte_Secondary", "Oocyte_Early antral", "Oocyte_Antral","Oocyte_Preovulatory" 
)
gNames <- c("Tgfb2","Tgfb3")
gNames <- c( "Dll3","Jag1","Jag2","Notch2")
gNames <- c('BMP7','BMP15','GDF9','INHBB','INHA') %>% str_to_title()
gNames <- c("Bmpr1a","Bmpr1b","Bmpr2","Acvr1a","Acvr1b","Acvr2a","Acvr2b","Tgfbr1","Tgfbr3")
data_plot <- FetchData(object = scdata,layer = 'count',clean = 'none',vars = c(gNames,'group_cell','group_stage')) %>% 
  rownames_to_column('sample') %>% 
  melt(id.var = c('sample','group_cell','group_stage'),variable.name = 'gene',value.name = 'expression') %>% 
  mutate(group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa')) %>%
  mutate(sample = factor(sample,levels = sample_level_group)) %>% 
  mutate(group_tmp = paste(group_cell,group_stage,sep = '_')) %>% 
  mutate(group_tmp = factor(group_tmp,levels = lavel_group_tmp)) %>% 
  as.data.frame()
ggplot(data = data_plot,mapping = aes(x = group_tmp,y = expression,fill = group_tmp)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.5),fill = NA) +
  facet_wrap(~gene,ncol = 3,scales = 'free_y',strip.position = 'left') +
  scale_fill_manual(values = my36colors) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = '',y = '') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 25),
    legend.position = 'top',
    legend.title = element_blank(),
    strip.text = element_text(size = 20),
    strip.placement = 'outside',
    strip.background = element_blank()
  )


scdata <- readRDS('./scdata_use.rds')
tmp <- FetchData(
  scdata,layer = 'counts',
  vars = rownames(scdata)
) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  left_join(protein_list,by = 'protein') %>% 
  dplyr::select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
scdata <- CreateSeuratObject(
  tmp, 
  meta.data = scdata@meta.data,
  min.cells = 2, 
  min.features = 2000)
scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)
all.genes <- rownames(scdata)
scdata <- ScaleData(scdata, features = all.genes)
gNames <- c("bmp15","gdf9","kit","dppa3","zp3","ddx4","gja4",
            "zp1","zp2","zar1","dazl","Rnf24",
            "Hsf1","Ncoa2","Slc39a10",
            "cnmd","Amh","Nr5a2","gja1","Hmgcs2","Foxl2",
            "Runx1","Amhr2") %>% str_to_title()
gNames <- c(
  "Bmp15","Gja4","Zp1","Zp2","Zar1",
  "Gja1","Amh","Hmgcs2","Foxl2",'Amhr2'
)
DotPlot(object = scdata,features = gNames,scale = FALSE,group.by = 'group_cell') +
  scale_color_gradient(
    name = 'Mean Exp',
    low = 'grey90',
    high = '#3131F2'
  ) +
  scale_size_continuous(name = 'Fraction') +
  guides(
    color = guide_colourbar(nrow = 1,title.position = 'top'),
    size = guide_legend(ncol = 1,title.position = 'top')
  ) +
  labs(x = 'Features', y= 'Cell Type') +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
    legend.direction = 'vertical',
    axis.title.x = element_blank()
  )


setwd('./protein')
library(stringr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(paletteer)
library(patchwork)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)
library(org.Mm.eg.db)
library(pathview)
library(paletteer)
data_raw <- read.csv('./20240315_MBR/report.tsv',sep = '\t')
data_raw<-data_raw[data_raw$Proteotypic == 1 & data_raw$Protein.Q.Value <= 0.01,]
colnames(data_raw)
data_protein<-data_raw[,c('File.Name',"Run" ,"Protein.Group" ,"Protein.Names","Genes","PG.Quantity")]
data_protein<-data_protein[!duplicated(data_protein),]
tmp <- data_raw %>% 
  filter(Genes %in% c("Tex19.2","Ezhip","H1.8","Brme1","Pgk2")) %>% 
  group_by(Genes) %>% 
  filter(PG.Quantity == max(PG.Quantity)) %>% 
  arrange(desc(Genes))
tmp$Genes %>% unique()
write.csv(tmp,'./file.csv')

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Seurat)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
setwd('./')





data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  dplyr::rename(protein = Accessions)
saveRDS(protein_list,'./protein_list.rds')

scdata_pro <- readRDS('./scdata_use.rds')
data_pro <- scdata_pro@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  filter(protein %in% protein_list$protein) %>% 
  left_join(protein_list,by = 'protein') %>% 
  dplyr::select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
tmp <- data_pro %>% 
  filter(rowSums(.>0)>1)
write.csv(data_pro,'./protein-copy_number-diann_mbr.csv')
meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')


scdata_tran <- readRDS('./scdata_tpm.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
write.csv(data_tran,'./transcription-tpm.csv')
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))

histone_data<-read.csv('./histone_gene_use.csv',sep=",",header = T,check.names = F) %>% 
  pull(names(.)[2])



gene_protein_coding <-
  read.csv('./gene_protein_coding.csv', row.names = 1)
head(gene_protein_coding)




data_plot_tran <- data_tran %>% 
  rownames_to_column('Genes') %>% 
  reshape2::melt(id.vars = 'Genes',variable.name = 'sample',value.name = 'abundance') %>% 
  left_join(meta_data_tran %>% dplyr::select(sample,group_omics,group_stage_cell),by = 'sample') %>% 
  group_by(Genes,group_stage_cell) %>% 
  summarise(
    abundance = mean(abundance)
  ) %>% 
  rename(abundance_tran = abundance)

data_plot_pro <- data_pro %>% 
  rownames_to_column('Genes') %>% 
  reshape2::melt(id.vars = 'Genes',variable.name = 'sample',value.name = 'abundance') %>% 
  left_join(meta_data_pro %>% dplyr::select(sample,group_omics,group_stage_cell),by = 'sample') %>% 
  group_by(Genes,group_stage_cell) %>% 
  summarise(
    abundance = mean(abundance)
  ) %>% 
  rename(abundance_pro = abundance)

pwthway_select <- read.csv('./fig2E_OC_GC.csv')

data_plot <- full_join(data_plot_tran, data_plot_pro, by = c("Genes", "group_stage_cell")) %>% 
  na.omit(.) %>% 
  mutate(
    abundance_tran = log10(abundance_tran+1),
    abundance_pro = log10(abundance_pro+1)
  ) %>% 
  filter(
    abundance_tran>0,abundance_pro>0,
    grepl('Oocyte',group_stage_cell)
  ) %>% 
  ungroup()

model <- lm(abundance_pro ~ abundance_tran, data = data_plot)


coefficients <- coef(model)


equation <- paste0("y = ", round(coefficients[2], 2), "x + ", round(coefficients[1], 2))

rho <- cor(data_plot$abundance_pro, data_plot$abundance_tran, method = "spearman")
rho_label <- paste('rho = ',round(rho,2),sep = '')



density_data <- density(data_plot$abundance_tran)

valid_indices <- which(density_data$x < 1.5 & density_data$x>0.1)

min_density_idx <- valid_indices[which.min(density_data$y[valid_indices])]
min_density_x <- density_data$x[min_density_idx]
min_density_x
min_density_y <- density_data$y[min_density_idx]

library(ggplot2)
library(dplyr)
library(viridis) 
library(ggpointdensity)
library(aplot)

model <- lm(abundance_pro ~ abundance_tran, data = data_plot)
intercept <- coef(model)[1]
slope <- coef(model)[2]
shift_distance <- 1

p1 <- ggplot(data = data_plot, mapping = aes(x = abundance_tran, y = abundance_pro)) +
  geom_pointdensity(adjust = 0.1) +
  geom_smooth(method = 'lm',formula = 'y~x',color = 'grey80') +
  geom_abline(intercept = intercept + shift_distance, slope = slope, color = 'red', linetype = "longdash",linewidth = 1.5) + 
  geom_abline(intercept = intercept - shift_distance, slope = slope, color = 'red', linetype = "longdash",linewidth = 1.5) +
  geom_vline(
    xintercept = min_density_x, 
    linetype = "dashed",linewidth = 1.5,
    color = "blue") +
  annotate(geom = 'text',x = 0.25,y=9.2,label = equation,hjust = 0,size = 8) +
  annotate(geom = 'text',x = 0.25,y=8.8,label = rho_label,hjust = 0,size = 8) +
  scale_color_viridis_c(breaks = c(10,280),labels = c('low','high')) +
  guides(
    color = guide_colorbar(
      title = 'Density',
      title.position = 'left',
      title.vjust = 1,
      label.position = 'bottom')
  ) +
  scale_y_continuous(
    breaks = c(2.5,5,7.5),
    expand = c(0,0)) +
  scale_x_continuous(
    breaks = c(0.5,1.5,2.5,3.5),
    expand = c(0,0)) +
  labs(x = 'Log10(Transcript)',y = 'Log10(Protein)') +
  theme_classic() +
  theme(
    legend.title = element_text(face = 'bold'),
    legend.text = element_text(size = 20),
    legend.direction = 'horizontal',
    legend.position = c(0.8,0.1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
p1

p2 <- ggplot(data_plot, aes(abundance_tran)) +
  geom_density(fill = "#AE123AFF", alpha = 0.5) +
  annotate(
    geom = 'segment',
    x = min_density_x, 
    xend = min_density_x,
    y = 0, yend = min_density_y, 
    linetype = "dashed",linewidth = 1.5,
    color = "blue") +
  scale_y_continuous(expand = c(0, 0),breaks = c(0.2,0.4,0.6)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )
p2
p3<-ggplot(data_plot,aes(abundance_pro))+
  geom_density(fill="#2A5783FF",alpha=0.5)+
  scale_y_continuous(expand = c(0,0),breaks = c(0.2,0.3,0.4))+
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )+
  coord_flip()
p_res <- p1%>%
  insert_top(p2,height = 0.2)%>%
  insert_right(p3,0.2)
p_res
ggsave(
  plot = p_res,
  filename = './fig2E_Oc.png',
  width = 15,
  height = 12
)
ggsave(
  plot = p_res,
  filename = './fig2E_Oc.pdf',
  width = 15,
  height = 12
)

data_plot_lowcorrelation <- data_plot %>% 
  mutate(group = ifelse(
    test = abundance_tran < min_density_x,
    yes = 'low correlation',
    no = 'others'
  ))
data_plot_lowcorrelation$group %>% table()
write.csv(data_plot_lowcorrelation,'./fig2f_Oocyte_lowcorrelation.csv')


library(clusterProfiler)
gene_selecy <- data_plot %>% 
  filter(abundance_tran<=min_density_x) %>% 
  pull(Genes) %>% unique()
gene_symbol <-
  bitr(
    gene_selecy,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = "org.Mm.eg.db"
  )



enrich.go <- enrichGO(
  gene = gene_symbol$ENTREZID,

  OrgDb = 'org.Mm.eg.db',
  keyType = 'ENTREZID',
  ont = 'BP',

  pAdjustMethod = 'fdr',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = FALSE
)
colnames(enrich.go@result)
enrich.go <- enrich.go
write.csv(enrich.go@result,'./fig2E_go_Oo.csv')
data_plot <- read.csv(
  './fig2E_go_Oo.csv',
  row.names = 1)
pathway_select_oo <- pwthway_select$OC %>% unique() %>% .[nchar(.)>0]
data_plot <- data_plot %>% 
  arrange(desc(Count)) %>% 
  filter(Description %in% pathway_select_oo) %>% 
  mutate(Description = str_wrap(Description,width = 50))  %>% 
  mutate(
    Description = factor(Description,levels = rev(Description))
  ) %>% 
  filter(pvalue<0.05) %>% 
  slice(1:15)
data_plot$Description
plot_go <- ggplot(data_plot, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot,
    aes(x = 1,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 0
  ) +
  scale_fill_gradient(
    low = 'white',
    high = '#4d627c',
    breaks = c(
      min(-log10(data_plot$pvalue)) %>% ceiling(), 
      max(-log10(data_plot$pvalue)) %>% floor()/2 + 
        min(-log10(data_plot$pvalue)) %>% ceiling()/2,
      max(-log10(data_plot$pvalue)) %>% floor()
    )
  ) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16),
    legend.text = element_text(size = 16)
  )
plot_go
ggsave(
  filename = './fig2E_go_Oo.png',
  bg = 'white',
  plot = plot_go,
  width = 10,
  height = 12,
  dpi = 300
)
data_plot <- full_join(data_plot_tran, data_plot_pro, by = c("Genes", "group_stage_cell")) %>% 
  na.omit(.) %>% 
  mutate(
    abundance_tran = log10(abundance_tran+1),
    abundance_pro = log10(abundance_pro+1)
  ) %>% 
  filter(
    abundance_tran>0,abundance_pro>0,
    grepl('Granulosa',group_stage_cell)
  ) %>% 
  ungroup()

model <- lm(abundance_pro ~ abundance_tran, data = data_plot)


coefficients <- coef(model)


equation <- paste0("y = ", round(coefficients[2], 2), "x + ", round(coefficients[1], 2))

rho <- cor(data_plot$abundance_pro, data_plot$abundance_tran, method = "spearman")
rho_label <- paste('rho = ',round(rho,2),sep = '')



density_data <- density(data_plot$abundance_tran)

valid_indices <- which(density_data$x < 1.5 & density_data$x>0.1)

min_density_idx <- valid_indices[which.min(density_data$y[valid_indices])]
min_density_x <- density_data$x[min_density_idx]
min_density_x
min_density_y <- density_data$y[min_density_idx]

library(ggplot2)
library(dplyr)
library(viridis) 
library(ggpointdensity)
library(aplot)

p1<-ggplot(data = data_plot, mapping = aes(x = abundance_tran, y = abundance_pro)) +
  geom_pointdensity(adjust = 0.1) +
  geom_smooth(method = 'lm',formula = 'y~x',color = 'grey80') +
  geom_vline(
    xintercept = min_density_x, 
    linetype = "dashed",linewidth = 1.5,
    color = "blue") +
  annotate(geom = 'text',x = 0.25,y=7.8,label = equation,hjust = 0,size = 8) +
  annotate(geom = 'text',x = 0.25,y=7.6,label = rho_label,hjust = 0,size = 8) +
  scale_color_viridis_c(breaks = c(10,280),labels = c('low','high')) +
  guides(
    color = guide_colorbar(
      title = 'Density',
      title.position = 'left',
      title.vjust = 1,
      label.position = 'bottom')
  ) +
  scale_y_continuous(
    breaks = c(2.5,5,7.5),
    expand = c(0,0)) +
  scale_x_continuous(
    breaks = c(0.5,1.5,2.5,3.5),
    expand = c(0,0)) +
  labs(x = 'Log10(Transcript)',y = 'Log10(Protein)') +
  theme_classic() +
  theme(
    legend.title = element_text(face = 'bold'),
    legend.text = element_text(size = 20),
    legend.direction = 'horizontal',
    legend.position = c(0.8,0.1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
p1

p2 <- ggplot(data_plot, aes(abundance_tran)) +
  geom_density(fill = "#AE123AFF", alpha = 0.5) +
  annotate(
    geom = 'segment',
    x = min_density_x, 
    xend = min_density_x,
    y = 0, yend = min_density_y, 
    linetype = "dashed",linewidth = 1.5,
    color = "blue") +
  scale_y_continuous(expand = c(0, 0),breaks = c(0.2,0.4,0.6)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )
p2
p3<-ggplot(data_plot,aes(abundance_pro))+
  geom_density(fill="#2A5783FF",alpha=0.5)+
  scale_y_continuous(expand = c(0,0),breaks = c(0.2,0.3,0.4))+
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )+
  coord_flip()
p_res <- p1%>%
  insert_top(p2,height = 0.2)%>%
  insert_right(p3,0.2)
p_res
ggsave(
  plot = p_res,
  filename = './Fig2E_gc.png',
  width = 15,
  height = 12)
ggsave(
  plot = p_res,
  filename = './Fig2E_gc.pdf',
  width = 15,
  height = 12)

data_plot_lowcorrelation <- data_plot %>% 
  mutate(group = ifelse(
    test = abundance_tran < min_density_x,
    yes = 'low correlation',
    no = 'others'
  ))
data_plot_lowcorrelation$group %>% table()
write.csv(data_plot_lowcorrelation,'./fig2f_GC_lowcorrelation.csv')

library(clusterProfiler)
gene_selecy <- data_plot %>% 
  filter(abundance_tran<=min_density_x) %>% 
  pull(Genes) %>% unique()
gene_symbol <-
  bitr(
    gene_selecy,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = "org.Mm.eg.db"
  )



enrich.go <- enrichGO(
  gene = gene_symbol$ENTREZID,

  OrgDb = 'org.Mm.eg.db',
  keyType = 'ENTREZID',
  ont = 'BP',

  pAdjustMethod = 'fdr',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = FALSE
)
colnames(enrich.go@result)
enrich.go <- enrich.go
write.csv(enrich.go@result,'./fig2E_go_GC.csv')
data_plot <- read.csv(
  './fig2E_go_GC.csv',
  row.names = 1)
pathway_select_gc <- pwthway_select$GC %>% unique() %>% .[nchar(.)>0]
data_plot <- data_plot %>% 
  arrange(desc(Count)) %>% 
  filter(Description %in% pathway_select_oo) %>% 
  mutate(Description = str_wrap(Description,width = 50))  %>% 
  mutate(
    Description = factor(Description,levels = rev(Description))
  ) %>% 
  filter(pvalue<0.05) %>% 
  slice(1:15)
data_plot$Description
plot_go <- ggplot(data_plot, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot,
    aes(x = 1,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 0
  ) +
  scale_fill_gradient(
    low = 'white',
    high = '#865254',
    breaks = c(
      min(-log10(data_plot$pvalue)) %>% ceiling(), 
      max(-log10(data_plot$pvalue)) %>% floor()/2 + 
        min(-log10(data_plot$pvalue)) %>% ceiling()/2,
      max(-log10(data_plot$pvalue)) %>% floor()
    )
  ) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.y = element_blank(),
    axis.title.x = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16),
    legend.text = element_text(size = 16)
  )
plot_go
ggsave(
  filename = './fig2E_go_GC.png',
  bg = 'white',
  plot = plot_go,
  width = 10,
  height = 12,
  dpi = 300
)

data_plot_tran <- data_tran %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rename_all(~c('abundance')) %>% 
  mutate(percentage = abundance/sum(abundance)) %>% 
  arrange(desc(percentage)) %>% 
  mutate(
    cumulative_intensity = cumsum(percentage),
    group_omits = 'transcription',
    rank = 1:nrow(.)
  ) %>% 
  rownames_to_column('Genes')
data_plot_pro <- data_pro %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rename_all(~c('abundance')) %>% 
  mutate(percentage = abundance/sum(abundance)) %>% 
  arrange(desc(percentage)) %>% 
  mutate(
    cumulative_intensity = cumsum(percentage),
    group_omits = 'protein',
    rank = 1:nrow(.)
  ) %>% 
  rownames_to_column('Genes')
data_plot <- rbind(data_plot_tran, data_plot_pro) %>% 
  ungroup()
data_plot <- data_plot_pro
data_plor_label <- data_plot %>% 
  arrange(group_omits,cumulative_intensity) %>% 
  group_by(group_omits) %>% 
  dplyr::filter(Genes %in% c("Ezhip","H1.8"))
ggplot(data = data_plot, aes(x = rank, y = cumulative_intensity, color = group_omits)) +
  geom_point() +
  geom_point(
    data = data_plor_label,
    aes(x = rank, y = cumulative_intensity),
    color = 'grey10',
    size = 2
  ) +
  geom_text_repel(
    data = data_plor_label,
    aes(x = rank, y = cumulative_intensity, label = Genes),
    color = 'grey10',
    show.legend = FALSE,
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0.5,
    segment.size  = 1,
    force = 400,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.03, "npc")),
    max.overlaps = Inf
  ) +
  facet_wrap(~ group_omits, nrow = 1, scales = 'free_x') +
  scale_color_manual(values = my36colors) +
  labs(x = 'Abundance rank', y = 'cumulative indensity') +
  theme_classic() +
  theme(
    legend.title = element_text(face = 'bold'),
    legend.text = element_text(size = 20),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20, margin = margin(b = 6, unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  './protein_abundance.png',
  width = 5,height = 8
)
ggsave(
  './protein_abundance.pdf',
  width = 5,height = 8
)

colnames(meta_data_tran)
group_stage_cell <- meta_data_tran$group_stage_cell %>% unique()
data_plot_tran <- lapply(group_stage_cell, function(group_use){
  sample <- meta_data_tran$sample[meta_data_tran$group_stage_cell %in% group_use]
  data_tmp <- data_tran[,colnames(data_tran) %in% sample] %>% 
    rownames_to_column('Genes')
  protein_all <- apply(data_tmp[,-1],2, function(sample){
    data_tmp$Genes[sample>0] %>% .[!is.na(.)]
  }) %>% unlist() %>% 
    table() %>% as.data.frame() %>% 
    rename_all(~c('protein','num')) %>% 
    filter(num>2) %>% 
    pull(protein) %>% as.character()
  protein_all <- protein_all[protein_all %in% gene_protein_coding$gene_name]
  return(protein_all)
}) 
names(data_plot_tran) <- group_stage_cell

colnames(meta_data_pro)
group_stage_cell <- meta_data_pro$group_stage_cell %>% unique()
group_use <- group_stage_cell[1]
data_plot_pro <- lapply(group_stage_cell, function(group_use){
  sample <- meta_data_pro$sample[meta_data_pro$group_stage_cell %in% group_use]
  data_tmp <- data_pro[,colnames(data_pro) %in% sample] %>% 
    rownames_to_column('Genes')
  protein_all <- apply(data_tmp[,-1],2, function(sample){
    data_tmp$Genes[sample>0] %>% .[!is.na(.)]
  }) %>% unlist() %>% 
    table() %>% as.data.frame() %>% 
    rename_all(~c('protein','num')) %>% 
    filter(num>2) %>% 
    pull(protein) %>% as.character()
  return(protein_all)
}) 
names(data_plot_pro) <- group_stage_cell
data_plot <- lapply(group_stage_cell, function(stage_cell){
  all_gene <- union(data_plot_tran[[stage_cell]],data_plot_pro[[stage_cell]]) %>% length()
  share_gene <- intersect(data_plot_tran[[stage_cell]],data_plot_pro[[stage_cell]]) %>% length()
  gene_tran <- setdiff(data_plot_tran[[stage_cell]],data_plot_pro[[stage_cell]]) %>% length()
  gene_pro <- setdiff(data_plot_pro[[stage_cell]],data_plot_tran[[stage_cell]]) %>% length()
  return(c(share_gene/all_gene,gene_tran/all_gene,gene_pro/all_gene))
}) %>% do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rename_all(~c('share_gene','gene_tran','gene_pro')) %>% 
  mutate(stage_cell = group_stage_cell) %>% 
  reshape2::melt(id.var = 'stage_cell',variable.name = 'group_gene',value.name = 'percentage') %>% 
  mutate(group_gene = factor(group_gene,levels = c(
    'gene_tran','share_gene','gene_pro'
    )))
ggplot(data = data_plot,mapping = aes(x = stage_cell,y = percentage*100,fill = group_gene)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(
    values = c('#0D7291','#FC7136','#716CAC'),
    labels = c("Transcription Only","Shared Genes","Protein Only"),
    guide = guide_legend(ncol = 1, label.position = "bottom")
    ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Abundance rank', y = '% of genes') +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    plot.margin = margin(l = 55,t = 10),
    legend.title = element_blank(),
    legend.position = 'right',
    legend.text = element_text(size = 20,angle = -90),
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 18),
    axis.text.y = element_text(size = 18),
    axis.line = element_blank()
  )
ggsave(
  filename = './D_co_exp.png',
  width = 8,
  height = 12
)
ggsave(
  filename = './D_co_exp.pdf',
  width = 8,
  height = 12
)


tmp_pro <- data_plot_pro[seq(1,8,2)] %>% unlist() %>% unique()
tmp_tran <- data_plot_tran[seq(1,8,2)] %>% unlist() %>% unique()
tmp_res <- tmp_pro[!(tmp_pro %in% tmp_tran)]
write.csv(tmp_res,'./fig2c_protein_unique_GC.csv')
data_plot <- list(
  Transcriptome = data_plot_tran[seq(1,8,2)] %>% unlist() %>% unique(),
  proteome = data_plot_pro[seq(1,8,2)] %>% unlist() %>% unique()
)
p_venn <- plot(
  euler(data_plot,shape = "ellipse"),
  quantities = list(
    type = c("percent", "counts"),
    col = "black",
    font = 2,
    cex = 2
  ),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 3
  ),
  edges = list(col = my36colors, lwd = 5, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
p_venn
pdf(
  file = './Venn_pro-tran_GC.pdf',
  width = 28,height = 20
)
p_venn
dev.off()

data_plot <- list(
  Transcriptome = data_plot_tran[seq(1,8,2)] %>% unlist() %>% unique(),
  proteome = data_plot_pro[seq(1,8,2)] %>% unlist() %>% unique()
)
data_plot %>% names()
gene_only_tran <- setdiff(data_plot[['Transcriptome']],data_plot[['proteome']])
gene_only_pro <- setdiff(data_plot[['proteome']],data_plot[['Transcriptome']])
data_plot <- data_pro %>% 
  dplyr::select(names(.)[!grepl('O',names(.))]) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_pro) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =4) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/4,
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  geom_pointdensity(size =4,adjust = 0.1) +
  scale_size_manual(values = c(1,2,3)) +
  labs(
    x = 'Log10(CN_Protein)',
    y = 'Identification frequency'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_pro_GC.png',
  width = 9,height = 8)
ggsave(
  './e_special_pro_GC.pdf',
  width = 9,height = 8)

data_plot <- data_tran %>% 
  dplyr::select(names(.)[!grepl('O',names(.))]) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_tran) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =2) +
  geom_pointdensity(size =2,adjust = 0.01) +
  scale_color_gradient(low = "#2B5C8A", high = "#9E3D22") +
  labs(
    x = 'Log10(TPM_Transcript)',
    y = 'Identification frequency'
  ) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/4,
    color = 'grey',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'grey',size = 2,
    linetype = "dashed") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_trans_GC.png',
  width = 9,height = 8)
ggsave(
  './e_special_trans_GC.pdf',
  width = 9,height = 8)

tmp_pro <- data_plot_pro[seq(2,8,2)] %>% unlist() %>% unique()
tmp_tran <- data_plot_tran[seq(2,8,2)] %>% unlist() %>% unique()
tmp_res <- tmp_pro[!(tmp_pro %in% tmp_tran)]
write.csv(tmp_res,'./fig2c_protein_unique_OO.csv')
data_plot <- list(
  Transcriptome = data_plot_tran[seq(2,8,2)] %>% unlist() %>% unique(),
  proteome = data_plot_pro[seq(2,8,2)] %>% unlist() %>% unique()
)

p_venn <- plot(
  euler(data_plot,shape = "ellipse"),
  quantities = list(
    type = c("percent", "counts"),
    col = "black",
    font = 2,
    cex = 2
  ),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 3
  ),
  edges = list(col = my36colors, lwd = 5, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
p_venn
pdf(
  file = './Venn_pro-tran_Oo.pdf',
  width = 28,height = 20
)
p_venn
dev.off()


data_plot <- list(
  Transcriptome = data_plot_tran[seq(2,8,2)] %>% unlist() %>% unique(),
  proteome = data_plot_pro[seq(2,8,2)] %>% unlist() %>% unique()
)
data_plot %>% names()
gene_only_tran <- setdiff(data_plot[['Transcriptome']],data_plot[['proteome']])
gene_only_pro <- setdiff(data_plot[['proteome']],data_plot[['Transcriptome']])
data_plot <- data_pro %>% 
  dplyr::select(names(.)[grepl('O',names(.))]) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_pro) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =4) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/4,
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  geom_pointdensity(size =4,adjust = 0.1) +
  scale_size_manual(values = c(1,2,3)) +
  labs(
    x = 'Log10(CN_Protein)',
    y = 'Identification frequency'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_pro_Oo.png',
  width = 9,height = 8)
ggsave(
  './e_special_pro_Oo.pdf',
  width = 9,height = 8)

data_plot <- data_tran %>% 
  dplyr::select(names(.)[grepl('O',names(.))]) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_tran) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =2) +
  geom_pointdensity(size =2,adjust = 0.01) +
  scale_color_gradient(low = "#2B5C8A", high = "#9E3D22") +
  labs(
    x = 'Log10(TPM_Transcript)',
    y = 'Identification frequency'
  ) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/4,
    color = 'grey',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'grey',size = 2,
    linetype = "dashed") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_trans_Oo.png',
  width = 9,height = 8)
ggsave(
  './e_special_trans_Oo.pdf',
  width = 9,height = 8)


data_plot <- list(
  Transcriptome = data_plot_tran %>% unlist() %>% unique(),
  proteome = data_plot_pro %>% unlist() %>% unique()
)
library(eulerr)
p_venn <- plot(
  euler(data_plot,shape = "ellipse"),
  quantities = list(
    type = c("percent", "counts"),
    col = "black",
    font = 2,
    cex = 2
  ),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 3
  ),
  edges = list(col = my36colors, lwd = 5, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
p_venn
ggsave(
  filename = './C_venn.png',
  bg = 'white',
  plot = p_venn,
  width = 15,
  height = 12
)

data_plot <- list(
  Transcriptome = data_plot_tran %>% unlist() %>% unique(),
  proteome = data_plot_pro %>% unlist() %>% unique()
)
data_plot %>% names()
gene_only_tran <- setdiff(data_plot[['Transcriptome']],data_plot[['proteome']])
gene_only_pro <- setdiff(data_plot[['proteome']],data_plot[['Transcriptome']])
data_plot <- data_pro %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_pro) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =4) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/2,
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'lightblue',size = 2,
    linetype = "dashed") +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  geom_pointdensity(size =4,adjust = 0.1) +
  scale_size_manual(values = c(1,2,3)) +
  labs(
    x = 'Log10(CN_Protein)',
    y = 'Identification frequency'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_pro_All.png',
  width = 9,height = 8)
ggsave(
  './e_special_pro_All.pdf',
  width = 9,height = 8)

data_plot <- data_tran %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% gene_only_tran) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  group_by(Genes) %>% 
  summarise(
    mean_exp = mean(copy_number) %>% log10(),
    frequency = sum(copy_number>0)
  ) %>% 
  arrange(desc(frequency),mean_exp)
ggplot(data = data_plot,aes(x = mean_exp,y = frequency)) +
  geom_point(size =2) +
  geom_pointdensity(size =2,adjust = 0.01) +
  scale_color_gradient(low = "#2B5C8A", high = "#9E3D22") +
  labs(
    x = 'Log10(TPM_Transcript)',
    y = 'Identification frequency'
  ) +
  geom_hline(
    yintercept = (colnames(scdata_pro) %>% length())/2,
    color = 'grey',size = 2,
    linetype = "dashed") +
  geom_vline(
    xintercept = data_plot$mean_exp %>% median(),
    color = 'grey',size = 2,
    linetype = "dashed") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  './e_special_trans_All.png',
  width = 9,height = 8)
ggsave(
  './e_special_trans_All.pdf',
  width = 9,height = 8)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = data_plot,
  col =my3color,
  alpha = 0.30,
  print.mode=c("raw","percent"),
  filename =NULL
)

pdf("venn.pdf")
grid.draw(venn.plot)
dev.off()


data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
protein_list <- data_pro_clean %>% 
  select(Genes,Accessions) %>% 
  rename(protein = Accessions)
data_pro <- scdata_pro@assays$RNA$data %>% as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  filter(protein %in% protein_list$protein) %>% 
  left_join(protein_list,by = 'protein') %>% 
  select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
colnames(meta_data_tran)
group_stage_cell <- meta_data_tran$group_stage_cell %>% unique() %>% 
  .[grepl('Oocyte',.)] %>% as.character()
group_stage_cell
group_use <- group_stage_cell[2]
data_plot_tran <- lapply(group_stage_cell, function(group_use){
  print(group_use)
  sample <- meta_data_tran$sample[meta_data_tran$group_stage_cell %in% group_use]
  data_tmp <- data_tran[,colnames(data_tran) %in% sample] %>% 
    rownames_to_column('Genes') %>% 
    rowwise() %>% 
    summarise(
      Genes = Genes,
      row_sum = sum(c_across(-Genes)),
      cell_num = sum(c_across(-Genes) > 0)
    ) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(
      group_gene = ifelse(
        test = cell_num >= (length(sample)/2),
        yes = 'high conf',
        no = 'low conf'
        )) %>% 
    rename_all(~c(
      'Genes',
      group_use %>% as.character(),
      paste('cell num of (',group_use %>% as.character(),')',sep = ''),
      paste('group of gene (',group_use %>% as.character(),')',sep = '')
      ))
  return(data_tmp)
}) %>% Reduce(function(x, y) merge(x, y, by = "Genes", all = TRUE), .) %>% 
  pivot_longer(cols = group_stage_cell,names_to = 'stage_cell',values_to = 'abundance')

colnames(meta_data_pro)
group_stage_cell <- meta_data_pro$group_stage_cell %>% unique() %>% 
  .[grepl('Oocyte',.)] %>% as.character()
group_use <- group_stage_cell[1]
data_plot_pro <- lapply(group_stage_cell, function(group_use){
  print(group_use)
  sample <- meta_data_pro$sample[meta_data_pro$group_stage_cell %in% group_use]
  data_tmp <- data_pro[,colnames(data_tran) %in% sample] %>% 
    rownames_to_column('Genes') %>% 
    rowwise() %>% 
    summarise(
      Genes = Genes,
      row_sum = sum(c_across(-Genes)),
      cell_num = sum(c_across(-Genes) > 0)
    ) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    mutate(
      group_gene = ifelse(
        test = cell_num >= (length(sample)/2),
        yes = 'high conf',
        no = 'low conf'
      )) %>% 
    rename_all(~c(
      'Genes',
      group_use %>% as.character(),
      paste('cell num of (',group_use %>% as.character(),')',sep = ''),
      paste('group of gene (',group_use %>% as.character(),')',sep = '')
    ))
  return(data_tmp)
}) %>% 
  Reduce(function(x, y) merge(x, y, by = "Genes", all = TRUE), .) %>% 
  pivot_longer(cols = group_stage_cell,names_to = 'stage_cell',values_to = 'abundance')


data_plot <-
  inner_join(
    data_plot_tran %>% as.data.frame(), 
    data_plot_pro %>% as.data.frame(), 
    by = c('Genes', 'stage_cell')) %>% 
  rename_with(~ gsub("\\.x$", "_tran", .), everything()) %>%
  rename_with(~ gsub("\\.y$", "_pro", .), everything()) %>% 
  mutate(stage_cell = factor(stage_cell,levels = stage_cell %>% unique())) %>% 
  arrange(stage_cell) %>% 
  group_by(Genes) %>% 
  summarise(
    spearman_rho = tryCatch(
      {cor.test(abundance_tran, abundance_pro, method = "pearson")$estimate},
      warning = function(w) NA,
      error = function(e) NA
    ),
    p_value = tryCatch(
      {cor.test(abundance_tran, abundance_pro, method = "pearson")$p.value},
      warning = function(w) NA,
      error = function(e) NA
    ),
    across(starts_with("cell num of"), first),
    across(starts_with("group of gene"), first),
    .groups = 'drop'
  ) %>% 
  filter(!is.na(spearman_rho)) %>% 
  mutate(
    group = case_when(
      spearman_rho < 0 ~ 'negative',
      spearman_rho > 0 ~ 'positive'
    )
  )
table(
  abs(
    data_plot$spearman_rho %>% as.numeric()
  ) + data_plot$p_value %>% as.numeric()
)
data_plot_pro[data_plot_pro$Genes == 'Nfyb','abundance']
data_plot_tran[data_plot_tran$Genes == 'Nfyb','abundance']


colnames(data_plot)
data_plot$group_p <- ifelse(
  test = data_plot$p_value<0.05,yes = 'significant',no = data_plot$group
)
tmp <- data_plot %>% 
  filter(group == 'significant') %>% 
  arrange(desc(spearman_rho))
data_plot$spearman_rho %>% unique() %>% length()
data_plot$group %>% table()
legend_label <- data_plot %>% 
  summarise(spearman_rho = median(spearman_rho))
breaks <- seq(
  from = data_plot$spearman_rho %>% min(), 
  to = data_plot$spearman_rho %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  dplyr::reframe(
    Counts = cut(spearman_rho, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(spearman_rho, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names(),
    breaks = breaks[-1] %>% round(.,2)
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique()),
    group_cor = ifelse(test = breaks>0,yes = 'positive',no = 'negative')
  ) %>% 
  filter(Counts>0) %>% 
  ungroup()
data_plot_tmp_sig <- data_plot %>% 
  filter(p_value<0.05) %>% 
  dplyr::reframe(
    Counts = cut(spearman_rho, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(spearman_rho, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names(),
    breaks = breaks[-1] %>% round(.,2)
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique()),
    group_cor = ifelse(test = breaks>0,yes = 'positive',no = 'negative')
  ) %>% 
  filter(Counts>0) %>% 
  ungroup()
ggplot(data = data_plot_tmp, aes(x = group, y = Counts,fill = group_cor)) +
  geom_bar(
    stat = 'identity',
    position = 'identity',
    width = 1,
    alpha = 0.9
  ) +
  geom_bar(
    data = data_plot_tmp_sig,
    aes(x = group, y = Counts),
    stat = 'identity',
    position = 'identity',
    width = 1,
    fill = 'grey30'
  ) +
  annotate(
    geom = 'text',
    x = 80,
    y = 60,
    label = paste('median = ', round(legend_label,2))
  ) +
  labs(x = 'pearson correlation', y = '# of genes') +
  scale_x_discrete(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = my36colors) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20, margin = margin(b = 6, unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )


colnames(data_plot)
data_plot_use <- data_plot %>% 
  mutate(spearman_rho = spearman_rho %>% as.numeric()) %>% 
  arrange(desc(spearman_rho)) %>% 
  mutate(rank_value = row_number()) %>% 
  mutate(
    group = case_when(
      spearman_rho < 0 ~ 'negative',
      spearman_rho > 0 ~ 'positive'
    ),
    group_use = case_when(
      rank_value > (max(rank_value)-100) ~ 'high positive',
      (spearman_rho>(-0.01)) & (spearman_rho<0.01) ~ 'low corelation',
      rank_value < 101 ~ 'high negative',
    )
  )
data_plot_use$group_use %>% table()
write.csv(data_plot_use,'./Fig2_i.csv')
head(data_plot_use)
ggplot(
  data = data_plot_use, 
  aes(x = rank_value, y = spearman_rho,fill = spearman_rho)
) +
  geom_bar(stat = 'identity',width = 1.1) +
  scale_fill_gradient2(low = "#1f77b4",mid = '#FEFFD9',high = "#d62728") +
  labs(y = 'correlation coefficient',x = 'RNA-Protein pairs') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )



func_gene_cor <- function(cell_select,gene_select){
  protein_select <- protein_list$protein[protein_list$Genes ==gene_select]
  data_exp_pro <- FetchData(
    object = scdata_pro,
    slot = 'data',
    vars = c(protein_select,"group_cell","group_stage")) %>% 
    filter(group_cell == cell_select) %>% 
    rename_all(~c('abundance',"group_cell","group_stage")) %>% 
    group_by(group_stage) %>% 
    summarise(
      abundance = mean(abundance)
    ) %>% 
    mutate(abundance = (abundance-min(abundance))/(max(abundance)-min(abundance))) %>%
    rename(abundance_pro = abundance) %>% 
    rownames_to_column(var = 'sample')
  data_exp_tran <- FetchData(
    object = scdata_tran,
    slot = 'data_tpm',
    vars = c(gene_select,"group_cell","group_stage")) %>% 
    filter(group_cell == cell_select) %>% 
    rename_all(~c('abundance',"group_cell","group_stage")) %>% 
    group_by(group_stage) %>% 
    summarise(
      abundance = mean(abundance)
    ) %>% 
    mutate(abundance = (abundance-min(abundance))/(max(abundance)-min(abundance))) %>%
    rename(abundance_tran = abundance) %>%
    rownames_to_column(var = 'sample')
  data_exp <- data_exp_pro %>% 
    left_join(data_exp_tran,by = 'sample')
  p <- ggplot(
    data = data_exp,aes(x = abundance_tran,y = abundance_pro)
  ) +
    geom_point() +
    theme_bw()
  return(p)
}
gene_select <- read.csv('./strong related gene 0825.csv')
filtered_data <- read.csv('./Fig2_i.csv') %>%
  filter(if_all(contains("group.of.gene"), ~ . == "high conf"))
tmp <- filtered_data %>% 
  filter(Genes %in% gene_select$gene)
table(gene_select$gene %in% filtered_data$Genes)
filtered_data$group_use %>% table()
gene_select$gene


cell_select <- 'Oocyte'
gene_select <- 'Tulp3'
func_gene_cor(cell_select = cell_select,gene_select = gene_select)


data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
protein_list <- data_pro_clean %>% 
  select(Genes,Accessions) %>% 
  rename(protein = Accessions)
data_pro <- scdata_pro@assays$RNA$data %>% as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  filter(protein %in% protein_list$protein) %>% 
  left_join(protein_list,by = 'protein') %>% 
  select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
colnames(meta_data_tran)
group_stage_cell <- meta_data_tran$group_stage_cell %>% unique() %>% 
  .[grepl('Oocyte',.)] %>% as.character()
group_stage_cell
group_use <- group_stage_cell[1]
data_plot_tran <- lapply(group_stage_cell, function(group_use){
  sample <- meta_data_tran$sample[meta_data_tran$group_stage_cell %in% group_use]
  data_tmp <- data_tran[,colnames(data_tran) %in% sample] %>% 
    rowSums() %>% 
    as.data.frame() %>% 
    rownames_to_column('Genes') %>% 
    rename_all(~c('Genes',group_use %>% as.character()))
  return(data_tmp)
}) %>% Reduce(function(x, y) merge(x, y, by = "Genes", all = TRUE), .) %>% 
  reshape2::melt(id.var = 'Genes',variable.name = 'stage_cell',value.name = 'abundance')

colnames(meta_data_pro)
group_stage_cell <- meta_data_pro$group_stage_cell %>% unique() %>% 
  .[grepl('Oocyte',.)] %>% as.character()
group_use <- group_stage_cell[1]
data_plot_pro <- lapply(group_stage_cell, function(group_use){
  sample <- meta_data_pro$sample[meta_data_pro$group_stage_cell %in% group_use]
  data_tmp <- data_pro[,colnames(data_pro) %in% sample] %>% 
    rowSums() %>% 
    as.data.frame() %>% 
    rownames_to_column('Genes') %>% 
    rename_all(~c('Genes',group_use %>% as.character()))
  return(data_tmp)
}) %>% Reduce(function(x, y) merge(x, y, by = "Genes", all = TRUE), .) %>% 
  reshape2::melt(id.var = 'Genes',variable.name = 'stage_cell',value.name = 'abundance')

data_plot <-
  inner_join(
    data_plot_tran, 
    data_plot_pro, 
    by = c('Genes', 'stage_cell')) %>% 
  rename_all(~c('Genes', 'stage_cell','abundance_tran','abundance_pro')) %>% 
  mutate(stage_cell = factor(stage_cell,levels = stage_cell %>% unique())) %>% 
  arrange(stage_cell) %>% 
  separate(col = 'stage_cell',into = c('stage','celltype'),sep = ':')
dat_plot_sub <- data_plot %>% 
  filter(stage == 'Secondary') %>%
  filter((abundance_tran + abundance_pro) != 0) %>% 
  mutate(
    abundance_tran = log10(abundance_tran+1),
    abundance_pro = log10(abundance_pro+1)
  )
ggplot(
  data = dat_plot_sub,
  aes(x = abundance_tran,y = abundance_pro)) +
  geom_point() +
  theme_bw()



library(tidyverse)
library(org.Mm.eg.db)


data_scale <- read.csv('./data_scale.csv',row.names = 1)
colnames(data_scale)
sample_meta_data <- read.csv('./meta_data_qc.csv',row.names = 1)
data_scale <- data_scale %>% 
  dplyr::select(all_of(c(
    names(.)[1:3],
    sample_meta_data$sample
  )))
data_gene_protein <- data_scale[,c("Genes","Accessions")]


data_cellular <- read.csv('./cellular component.csv')

GOID <- c("GO:0005739")
data_cellular_gene <- lapply(data_cellular$GO, function(GOID){
  GOID <- GOID %>% str_remove_all(' ')
  GOgeneID <- get(GOID, org.Mm.egGO2ALLEGS) %>% 
    mget(org.Mm.egSYMBOL) %>% unlist() %>% as.character()
  return(GOgeneID)
})
names(data_cellular_gene) <- data_cellular$Cellular.component
data_cellular_protein <- lapply(names(data_cellular_gene), function(cellular){
  gene_all <- data_cellular_gene[[cellular]]
  protein_use <- data_gene_protein$Accessions[
    match(gene_all,data_gene_protein$Genes)
  ] %>% .[!is.na(.)]
})
names(data_cellular_protein) <- names(data_cellular_gene)


data_plot <- lapply(names(data_cellular_protein), function(cellular){
  protein_use <- data_cellular_protein[[cellular]]
  data_exp <- data_scale %>% 
    filter(Accessions %in% protein_use) %>% 
    dplyr::select(-c(Genes,MolecularWeight)) %>% 
    column_to_rownames('Accessions') %>% 
    t() %>% as.data.frame() %>% 
    rowSums() %>% as.data.frame() %>% 
    rename_all(~c('exp')) %>% 
    rownames_to_column('sample') %>% 
    left_join(
      sample_meta_data %>% dplyr::select(sample,group_cell,group_stage_cell),
      by = 'sample'
    ) %>% 
    column_to_rownames('sample') %>% 
    filter(group_cell == 'Oocyte') %>% 
    group_by(group_stage_cell) %>% 
    summarise(exp_all = sum(exp)) %>% 
    mutate(cellular = cellular)
    return(data_exp)
}) %>% do.call(rbind,.) %>% 
  mutate(
    cellular = factor(cellular,levels = c(
      'white',"Nucleus","Spindle", "Mitochondria", 
      "Cytoplasm", "Golgi apparatus", "Cytosol",
      "Membrane","Ribosome", "Nucleoplasm",
      "Endoplasmic Reticulum","Cytoskeleton"
    ))
  ) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    'all',"Secondary: Oocyte","Early antral: Oocyte",
    "Antral: Oocyte", "Preovulatory: Oocyte"
    ))) %>% 
  group_by(group_stage_cell) %>%
  mutate(pencent = exp_all/sum(exp_all)*100) %>% 
  arrange(desc(cellular)) %>% 
  mutate(
    pencent_label = pencent %>% 
      round(.,digits = 1) %>% paste('%',sep = '')
    ) %>% 
  group_by(group_stage_cell) %>%
  mutate(
    label_position = cumsum(pencent)-pencent/2
  )

data_plot <- data_plot %>% as.data.frame()
data_plot[nrow(data_plot)+1,] <- c('all',1,'white',NA,NA,NA)
data_plot <- data_plot %>% 
  mutate(
    exp_all = exp_all %>% as.numeric(),
    pencent = pencent %>% as.numeric(),
    label_position = label_position %>% as.numeric()
  )
colnames(data_plot)
ggplot(data_plot, aes(
  x = as.numeric(group_stage_cell),
  y = pencent,
  fill = cellular
)) +
  geom_bar(
    aes(alpha = group_stage_cell),
    stat = 'identity',

    position = 'stack',

    width = 0.8,

    color = 'white'
  ) +
  geom_text(aes(
    x = as.numeric(group_stage_cell) + 0.2,
    y = label_position,
    label = pencent_label
  )) +
  coord_polar(theta = 'y') +
  scale_alpha_manual(values = c(0.6,0.75,0.9,1,0)) +
  guides(alpha = 'none') +
  scale_fill_manual(
    values = c(paletteer_d("ggthemes::Tableau_20")),
    breaks = c(
      "Nucleus",
      "Spindle",
      "Mitochondria",
      "Cytoplasm",
      "Golgi apparatus",
      "Cytosol",
      "Membrane",
      "Ribosome",
      "Nucleoplasm",
      "Endoplasmic Reticulum",
      "Cytoskeleton"
    ),
    na.value = 'white',
    guide = guide_legend(title = 'Cellular')
  ) +
  labs(title = 'Oocyte') +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20),
    legend.position = c(1.05,0.5),
    plot.margin = margin(t = 0,r = 150,b = 0,l = 0),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = -12)
  )

data_plot <- lapply(names(data_cellular_protein), function(cellular){
  protein_use <- data_cellular_protein[[cellular]]
  data_exp <- data_scale %>% 
    filter(Accessions %in% protein_use) %>% 
    dplyr::select(-c(Genes,MolecularWeight)) %>% 
    column_to_rownames('Accessions') %>% 
    t() %>% as.data.frame() %>% 
    rowSums() %>% as.data.frame() %>% 
    rename_all(~c('exp')) %>% 
    rownames_to_column('sample') %>% 
    left_join(
      sample_meta_data %>% dplyr::select(sample,group_cell,group_stage_cell),
      by = 'sample'
    ) %>% 
    column_to_rownames('sample') %>% 
    filter(group_cell == 'Granulosa') %>% 
    group_by(group_stage_cell) %>% 
    summarise(exp_all = sum(exp)) %>% 
    mutate(cellular = cellular)
  return(data_exp)
}) %>% do.call(rbind,.) %>% 
  mutate(
    cellular = factor(cellular,levels = c(
      'white',"Nucleus","Spindle", "Mitochondria", 
      "Cytoplasm", "Golgi apparatus", "Cytosol",
      "Membrane","Ribosome", "Nucleoplasm",
      "Endoplasmic Reticulum","Cytoskeleton"
    ))
  ) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    'all',"Secondary: Granulosa", "Early antral: Granulosa",
    "Antral: Granulosa", "Preovulatory: Granulosa"
  ))) %>% 
  group_by(group_stage_cell) %>%
  mutate(pencent = exp_all/sum(exp_all)*100) %>% 
  arrange(desc(cellular)) %>% 
  mutate(
    pencent_label = pencent %>% 
      round(.,digits = 1) %>% paste('%',sep = '')
  ) %>% 
  group_by(group_stage_cell) %>%
  mutate(
    label_position = cumsum(pencent)-pencent/2
  )

data_plot <- data_plot %>% as.data.frame()
data_plot[nrow(data_plot)+1,] <- c('all',1,'white',NA,NA,NA)
data_plot <- data_plot %>% 
  mutate(
    exp_all = exp_all %>% as.numeric(),
    pencent = pencent %>% as.numeric(),
    label_position = label_position %>% as.numeric()
  )
colnames(data_plot)
ggplot(data_plot, aes(
  x = as.numeric(group_stage_cell),
  y = pencent,
  fill = cellular
)) +
  geom_bar(
    aes(alpha = group_stage_cell),
    stat = 'identity',

    position = 'stack',

    width = 0.8,

    color = 'white'
  ) +
  geom_text(aes(
    x = as.numeric(group_stage_cell) + 0.2,
    y = label_position,
    label = pencent_label
  )) +
  coord_polar(theta = 'y') +
  scale_alpha_manual(values = c(0.6,0.75,0.9,1,0)) +
  guides(alpha = 'none') +
  scale_fill_manual(
    values = c(paletteer_d("ggthemes::Tableau_20")),
    breaks = c(
      "Nucleus",
      "Spindle",
      "Mitochondria",
      "Cytoplasm",
      "Golgi apparatus",
      "Cytosol",
      "Membrane",
      "Ribosome",
      "Nucleoplasm",
      "Endoplasmic Reticulum",
      "Cytoskeleton"
    ),
    na.value = 'white',
    guide = guide_legend(title = 'Cellular')
  ) +
  labs(title = 'Granulosa cell') +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 20),
    legend.position = c(1.05,0.5),
    plot.margin = margin(t = 0,r = 150,b = 0,l = 0),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = -12)
  )

library(stringr)
library(msigdbr)
library(dplyr)
library(Seurat)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)


gene_select_Oo <- read.csv(
  './fig2c_protein_unique_OO.csv',
  row.names = 1) %>% pull(var = 1) %>% unique()
gene_select_GC <- read.csv(
  './fig2c_protein_unique_GC.csv',
  row.names = 1) %>% pull(var = 1) %>% unique()

my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
setwd('./')




data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  rename(protein = Accessions)

scdata_pro <- readRDS('./scdata_use.rds')
data_pro <- scdata_pro@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  filter(protein %in% protein_list$protein) %>% 
  left_join(protein_list,by = 'protein') %>% 
  dplyr::select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')

data_plot_pro <- data_pro %>% 
  rownames_to_column('Genes') %>%  
  reshape2::melt(
    id.vars = 'Genes',
    measure.vars = colnames(.)[2:ncol(.)],
    variable.name = 'sample',
    value.name = 'copy_num'
  ) %>% 
  filter(copy_num>0) %>% 
  left_join(meta_data_pro,by = 'sample') %>% 
  mutate(
    log10_cp = log10(copy_num),
    group_omics = 'protein'
  ) %>% 
  select(
    Genes,sample,copy_num,group_cell,
    group_stage,group_sample,log10_cp, group_omics
  )
data_plot_pro_Oo <- data_plot_pro %>% 
  filter(group_cell == 'Oocyte')

data_plot <- rbind(
  data_plot_pro %>% 
    filter(group_cell == 'Oocyte') %>% 
    mutate(group = 'All Protein groups'),
  data_plot_pro %>% 
    filter(group_cell == 'Oocyte',Genes %in% gene_select_Oo) %>% 
    mutate(group = 'Detected only in proteins')
)

ggplot() +   
  geom_histogram(
    data = data_plot_pro_Oo,
    mapping = aes(x = log10_cp,y = ..density..),
    bins = 100, alpha = 0.3,color = 'white',
    position = 'identity'
  ) +
  geom_density(
    data = data_plot, 
    mapping = aes(x = log10_cp,color = group),
    alpha = 0.8,linewidth = 2,key_glyph = "path"
  ) +
  geom_vline(
    xintercept = median(data_plot_pro_Oo$log10_cp),
    color = 'black',linetype = 'longdash',
    linewidth = 2
  ) +
  annotate(
    geom = 'text',
    x = (median(data_plot_pro_Oo$log10_cp)-0.6),
    y= 0.01,label = paste('Median:',median(data_plot_pro_Oo$log10_cp) %>% round(.,digits = 2)),
    size = 6,
    color = 'black'
  ) + 
  guides(color = guide_legend(override.aes = list(linetype = 1,linewidth = 2))) +
  scale_color_manual(values = c(
    '#00468BFF','#ED0000FF','#42B540FF'
  )) +
  scale_fill_manual(values = c(
    '#00468BFF','#ED0000FF','#42B540FF'
  )) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Log10(abundance)", y = "Density") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = 'top',
    legend.direction = 'vertical',
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig2f_proteins_Detected only in proteins_OC.pdf',
  width = 6,height = 6
)
data_plot_pro_Granulosa <- data_plot_pro %>% 
  filter(group_cell == 'Granulosa')

data_plot <- rbind(
  data_plot_pro %>% 
    filter(group_cell == 'Granulosa') %>% 
    mutate(group = 'All Protein groups'),
  data_plot_pro %>% 
    filter(group_cell == 'Granulosa',Genes %in% gene_select_GC) %>% 
    mutate(group = 'Detected only in proteins')
)

ggplot() +   
  geom_histogram(
    data = data_plot_pro_Granulosa,
    mapping = aes(x = log10_cp,y = ..density..),
    bins = 100, alpha = 0.3,color = 'white',
    position = 'identity',show.legend = FALSE
  ) +
  geom_density(
    data = data_plot, 
    mapping = aes(x = log10_cp,color = group),
    alpha = 0.3,linewidth = 2,key_glyph = "path"
  ) +
  geom_vline(
    xintercept = median(data_plot_pro_Granulosa$log10_cp),
    color = 'black',linetype = 'longdash',linewidth = 2
  ) +
  annotate(
    geom = 'text',
    x = (median(data_plot_pro_Granulosa$log10_cp)-0.6),
    y= 0.01,label = paste('Median:',median(data_plot_pro_Granulosa$log10_cp) %>% round(.,digits = 2)),
    size = 6,color = 'black'
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1,linewidth = 2))) +
  scale_color_manual(values = c(
    '#00468BFF','#ED0000FF','#42B540FF'
  )) +
  scale_fill_manual(values = c(
    '#00468BFF','#ED0000FF','#42B540FF'
  )) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Log10(abundance)", y = "Density") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = 'top',
    legend.direction = 'vertical',
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig2f_proteins_Detected only in proteins_GC.pdf',
  width = 6,height = 6
)
ggsave(
  filename = './fig2f_proteins_Detected only in proteins_GC.png',
  width = 6,height = 6
)



data_Maternal <- read.csv('./maternal_gene.txt',header = FALSE)
source_species_dataset <- 'mmusculus_gene_ensembl'

library(biomaRt)


human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")


gene_list <- data_Maternal$V1


homology <- getLDS(
  attributes = c('hgnc_symbol'),
  filters = 'hgnc_symbol',
  values = gene_list,
  mart = human_mart,
  attributesL = c('external_gene_name', 'ensembl_gene_id'),
  martL = mouse_mart,
  uniqueRows = TRUE
)


head(homology)




tmp_res <- read.csv('./fig2c_protein_unique_OO.csv',row.names = 1)

data_plot <- list(
  `Maternal genes` = homology$Gene.name %>% unique(),
  proteome = tmp_res$x %>% unique()
)

p_venn <- plot(
  euler(data_plot,shape = "ellipse"),
  quantities = list(
    type = c("percent", "counts"),
    col = "black",
    font = 2,
    cex = 2
  ),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 3
  ),
  edges = list(col = my36colors, lwd = 5, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
p_venn

library(clusterProfiler)
gene_select <- tmp_res$x %>% unique()
gene_symbol <-
  bitr(
    gene_select,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = "org.Mm.eg.db"
  )



enrich.go <- enrichGO(
  gene = gene_symbol$SYMBOL,

  OrgDb = 'org.Mm.eg.db',
  keyType = 'SYMBOL',
  ont = 'BP',

  pAdjustMethod = 'fdr',
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = FALSE
)
colnames(enrich.go@result)
result %>% dim()
result <- enrich.go@result %>% 
  filter(p.adjust<0.05)
write.csv(result,'./result_Oocte_special_protein.csv')
result$Description
library(aPEAR)
p <- enrichmentNetwork(
  enrich.go@result %>% filter(p.adjust<0.05),
  drawEllipses = TRUE, fontSize = 5,minClusterSize = 2,
  plotOnly = FALSE
)
p$plot
p$clusters$Cluster %>% unique()


data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  dplyr::rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  dplyr::rename(protein = Accessions)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
colnames(data_pro_clean)
scdata_pro <- readRDS('./scdata_use.rds')
meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')
data_singlestage_raw <- FetchData(
  scdata_pro,layer = 'counts',
  vars = rownames(scdata_pro)
) %>% 
  t() %>% as.data.frame()
data_singlestage <-  
  data_singlestage_raw %>% t() %>% as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  left_join(meta_data_pro %>% dplyr::select(c('sample','group_stage_cell')),by = 'sample') %>% 
  group_by(group_stage_cell) %>% 
  mutate(across(.cols = where(is.numeric),.fns = replace_outliers)) %>% 
  ungroup() %>% 
  dplyr::select(-group_stage_cell) %>% 
  column_to_rownames('sample') %>% 
  t() %>% as.data.frame()

scdata_pro <- CreateSeuratObject(
  counts = data_singlestage %>% as.matrix(),
  data = data_singlestage %>% as.matrix(),
  meta.data = scdata_pro@meta.data
)
saveRDS(scdata_pro,'./scdata_use_filteroutliers.rds')
data_pro <- scdata_pro@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column('protein') %>% 
  filter(protein %in% protein_list$protein) %>% 
  left_join(protein_list,by = 'protein') %>% 
  dplyr::select(-c('protein')) %>% 
  group_by(Genes) %>% 
  summarise(across(.cols = everything(), .fns = mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  column_to_rownames(var = 'Genes')
scdata_tran <- readRDS('./scdata_tpm.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_omics = 'transcription')
data_singlestage_raw <- FetchData(
  scdata_tran,layer = 'data_tpm',
  vars = rownames(scdata_tran)
) %>% 
  t() %>% as.data.frame()
data_singlestage <-  
  data_singlestage_raw %>% t() %>% as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  left_join(meta_data_tran %>% dplyr::select(c('sample','group_stage_cell')),by = 'sample') %>% 
  group_by(group_stage_cell) %>% 
  mutate(across(.cols = where(is.numeric),.fns = replace_outliers)) %>% 
  ungroup() %>% 
  dplyr::select(-group_stage_cell) %>% 
  column_to_rownames('sample') %>% 
  t() %>% as.data.frame()
scdata_tran <- CreateSeuratObject(
  counts = data_singlestage %>% as.matrix(),
  data = data_singlestage %>% as.matrix(),
  meta.data = scdata_tran@meta.data
)
saveRDS(scdata_tran,'./scdata_tpm_filteroutliters.rds')

meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
histone_data<-read.csv('./histone_gene_use.csv',sep=",",header = T,check.names = F) %>% 
  pull(names(.)[2])


gene_protein_coding <-
  read.csv('./gene_protein_coding.csv', row.names = 1)
head(gene_protein_coding)


pvalue_calculate <- function(adata = scdata_pro,group_stage = NA,sample_num = 10){
  library(coin)
  data_singlestage <- FetchData(
    adata,layer = 'counts',
    vars = rownames(adata)
  ) %>% 
    t() %>% as.data.frame() 
  gene_select <- data_singlestage %>% 
    rownames_to_column('Protein') %>% 
    rowwise() %>%
    filter(sum(c_across(-Protein) > 0) >= sample_num) %>%
    ungroup() %>%
    pull(Protein)
  data_singlestage <- data_singlestage[rownames(data_singlestage) %in% gene_select,]
  p_Permutation <- lapply(gene_select, function(protein){
    print(protein)
    tmp <- data.frame(
      value = data_singlestage[rownames(data_singlestage) == protein,] %>% as.numeric(),
      group = adata$group_stage
    )
    p_value <- kruskal.test(value ~ group, data = tmp)$p.value
    return(p_value)
  })
  p_values <- p_Permutation %>%
    as.data.frame() %>% t() %>% 
    as.data.frame() %>% 
    rename_all(~c('pvalue_Permutation')) %>% 
    mutate(Genes = gene_select) %>%
    rownames_to_column('row_name') %>% 
    select(-row_name)
  return(p_values)
}
scdata_tran$group_cell %>% unique()
scdata_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = data_gene_protein$Genes)
protein_use <- data_gene_protein$Accessions[match(rownames(scdata_sub),data_gene_protein$Genes)]
rownames(scdata_sub)[1:10]
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
data_tran_p <-  pvalue_calculate(adata = scdata_sub) %>%
  mutate(omics = 'Transcription') %>%
  column_to_rownames('Genes') %>% 
  rownames_to_column('Genes')
scdata_pro$group_cell %>% unique()
meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')

scdata_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = protein_use)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()

data_pro_p <-  pvalue_calculate(adata = scdata_sub) %>%
  mutate(omics = 'Protein') %>%
  rename('Accessions' = 'Genes') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  select(-c('Accessions')) %>% 
  column_to_rownames('Genes') %>% 
  rownames_to_column('Genes')
data_p <- data_pro_p %>% 
  full_join(data_tran_p,by = 'Genes',suffix = c("_pro", "_tran")) %>% 
  mutate(
    group = case_when(
      (pvalue_Permutation_pro < 0.05 & pvalue_Permutation_tran < 0.05) ~ 'Sig',
      (pvalue_Permutation_pro < 0.05 & pvalue_Permutation_tran > 0.05) ~ 'Only pro Sig',
      (pvalue_Permutation_pro > 0.05 & pvalue_Permutation_tran < 0.05) ~ 'Only tran Sig',
      (pvalue_Permutation_pro > 0.05 & pvalue_Permutation_tran > 0.05) ~ 'NS',
    )
  )
data_p$group %>% table()
write.csv(data_p,'./data_p_kruskal_test_Oocyte.csv')
data_p <- read.csv('./data_p_kruskal_test_Oocyte.csv',row.names = 1)
data_p <- data_p %>% filter(!is.na(group))
gene_use <- c(
  'Wtap','Npm1','Nelfa','Zcchc8','Bod1','Saraf','Cks2','Spry4','ska1','Gdf9','Dnmt1','Dnmt3a'
)
data_plot_select_point <- data_p %>% 
  dplyr::filter(Genes %in% gene_use)
p1 <- ggplot(
  data = data_p,
  aes(
    x = -log10(pvalue_Permutation_pro),
    y = -log10(pvalue_Permutation_tran),
    color = group,
  )) +
  geom_point() +
  geom_text_repel(
    data = data_plot_select_point,
    mapping = aes(
      x = -log10(pvalue_Permutation_pro),
      y = -log10(pvalue_Permutation_tran),
      label = Genes
    ),
    color = 'black',
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0.5,
    segment.size  = 1,
    force = 400,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.03, "npc")),
    max.overlaps = Inf
  ) +
  geom_point(
    data = data_plot_select_point,mapping = aes(
      x = -log10(pvalue_Permutation_pro),
      y = -log10(pvalue_Permutation_tran)
    ),
    color = 'black',size = 3
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  scale_size_manual(values = c(1,2,3)) +
  scale_color_manual(values = paletteer_d("ggsci::nrc_npg")) +
  labs(
    x = '-Log10(Pvalue of Protein)',
    y = '-Log10(Pvalue of Transcript)'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'top',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
data_plot_summary <- data_p$group %>% table() %>% 
  as.data.frame() %>% 
  dplyr::rename(group = names(.)[1]) %>% 
  mutate(
    prop = Freq/sum(Freq),
    group = factor(group,levels = group %>% unique() %>% rev())
  ) %>% 
  arrange(rev(group)) %>% 
  mutate(
    label = paste(Freq,'(',round(prop,digits = 2)*100,'%)',sep = ''),
    y_loc = prop/2 + c(0,cumsum(prop)[-length(cumsum(prop))])
  )
data_plot_summary
p2 <- ggplot(data = data_plot_summary,aes(x = 1,y = prop,fill = group)) +
  geom_bar(stat = 'identity',color = 'white',show.legend = FALSE) +
  geom_text(aes(x = 1,y = y_loc,label = label),color = 'white') +
  scale_fill_manual(values = paletteer_d("ggsci::nrc_npg")[1:4] %>% rev()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.border = element_rect(fill = NA,color = 'white'),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'right',
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
cowplot::plot_grid(p1,p2,align = 'h',nrow = 1,rel_widths = c(1,0.1))
ggsave('./pvalue_Oocyte.png',width = 12,height = 12)
ggsave('./pvalue_Oocyte.pdf',width = 12,height = 12)
scdata_tran$group_cell %>% unique()
scdata_sub <- subset(scdata_tran,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = data_gene_protein$Genes)
protein_use <- data_gene_protein$Accessions[match(rownames(scdata_sub),data_gene_protein$Genes)]
rownames(scdata_sub)[1:10]
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
data_tran_p <-  pvalue_calculate(adata = scdata_sub) %>%
  mutate(omics = 'Transcription') %>%
  column_to_rownames('Genes') %>% 
  rownames_to_column('Genes')
scdata_pro$group_cell %>% unique()
meta_data_pro <- scdata_pro@meta.data %>% 
  mutate(group_omics = 'protein')

scdata_sub <- subset(scdata_pro,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = protein_use)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()

data_pro_p <-  pvalue_calculate(adata = scdata_sub) %>%
  mutate(omics = 'Protein') %>%
  rename('Accessions' = 'Genes') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  select(-c('Accessions')) %>% 
  column_to_rownames('Genes') %>% 
  rownames_to_column('Genes')
data_p <- data_pro_p %>% 
  full_join(data_tran_p,by = 'Genes',suffix = c("_pro", "_tran")) %>% 
  mutate(
    group = case_when(
      (pvalue_Permutation_pro < 0.05 & pvalue_Permutation_tran < 0.05) ~ 'Sig',
      (pvalue_Permutation_pro < 0.05 & pvalue_Permutation_tran > 0.05) ~ 'Only pro Sig',
      (pvalue_Permutation_pro > 0.05 & pvalue_Permutation_tran < 0.05) ~ 'Only tran Sig',
      (pvalue_Permutation_pro > 0.05 & pvalue_Permutation_tran > 0.05) ~ 'NS',
    )
  )
data_p$group %>% table()
write.csv(data_p,'./data_p_kruskal_test_Granulosa.csv')
data_p <- read.csv('./data_p_kruskal_test_Granulosa.csv')
data_p <- data_p %>% filter(!is.na(group))
data_p$group %>% table()
gene_use <- c(
  "Ucp2","Mrpl23","Utp11","Hnrnpa0","Map2k3","Hmga1",
  "N6amt1","Dpysl2","Chmp1a","Eif5","Sema7a","Cdc73",
  "Eif4g3","Inhba","Khdc3"
)
data_plot_select_point <- data_p %>%
  dplyr::filter(Genes %in% gene_use)
p1 <- ggplot(
  data = data_p,
  aes(
    x = -log10(pvalue_Permutation_pro),
    y = -log10(pvalue_Permutation_tran),
    color = group,
  )) +
  geom_point() +
  geom_text_repel(
    data = data_plot_select_point,
    mapping = aes(
      x = -log10(pvalue_Permutation_pro),
      y = -log10(pvalue_Permutation_tran),
      label = Genes
    ),
    color = 'black',
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0,
    segment.size  = 1,
    force = 30,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = Inf
  ) +
  geom_point(
    data = data_plot_select_point,mapping = aes(
      x = -log10(pvalue_Permutation_pro),
      y = -log10(pvalue_Permutation_tran)
    ),
    color = 'black'
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3))
  ) +
  scale_size_manual(values = c(1,2,3)) +
  scale_color_manual(values = paletteer_d("ggsci::nrc_npg")) +
  labs(
    x = '-Log10(Pvalue of Protein)',
    y = '-Log10(Pvalue of Transcript)'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'top',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
p1
data_plot_summary <- data_p$group %>% table() %>% 
  as.data.frame() %>% 
  dplyr::rename(group = names(.)[1]) %>% 
  mutate(
    prop = Freq/sum(Freq),
    group = factor(group,levels = group %>% unique() %>% rev())
  ) %>% 
  arrange(rev(group)) %>% 
  mutate(
    label = paste(Freq,'(',round(prop,digits = 2)*100,'%)',sep = ''),
    y_loc = prop/2 + c(0,cumsum(prop)[-length(cumsum(prop))])
  )
data_plot_summary
p2 <- ggplot(data = data_plot_summary,aes(x = 1,y = prop,fill = group)) +
  geom_bar(stat = 'identity',color = 'white',show.legend = FALSE) +
  geom_text(aes(x = 1,y = y_loc,label = label),color = 'white') +
  scale_fill_manual(values = paletteer_d("ggsci::nrc_npg")[1:4] %>% rev()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    panel.border = element_rect(fill = NA,color = 'white'),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'right',
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
cowplot::plot_grid(p1,p2,align = 'h',nrow = 1,rel_widths = c(1,0.1))

ggsave('./pvalue_Granulosa.png',width = 12,height = 12)
ggsave('./pvalue_Granulosa.pdf',width = 12,height = 12)

library(tidyr)
library(tibble)
library(ggplot2)
library(dplyr)
library(Seurat)
setwd('./')

scdata_pro <- readRDS('./scdata_use.rds')
data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  rename(Accessions = ProteinGroups)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
colnames(data_pro_clean)

scdata_tran <- readRDS('./scdata_tpm.rds')
color_group_stage_cell <- c('#C34F73','#0088D4','#4B8600','#5F559B',
                            '#E64B35','#4DBBD5','#00A087','#3C5488')
names(color_group_stage_cell) <- c(
  "Secondary: Granulosa", "Early antral: Granulosa",
  "Antral: Granulosa", "Preovulatory: Granulosa",
  "Secondary: Oocyte","Early antral: Oocyte",
  "Antral: Oocyte", "Preovulatory: Oocyte"
)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 


data_plot <- scdata_pro@assays$RNA$counts %>% 
  as.data.frame() %>% 
  apply(., 1, function(x){sum(x>0)}) %>% 
  as.data.frame() %>% 
  rename_all(~c('count')) %>% 
  table() %>% as.data.frame() %>% 
  rename_all(~c('sample_num','count'))
ggplot(data_plot,aes(x = sample_num,y = count)) +
  geom_bar(stat = 'identity')+theme_classic()
gene_select <- scdata_pro@assays$RNA$counts %>% 
  rownames(.)

group_stage <- 'Secondary'
cv_calculate <- function(adata = scdata_pro,group_stage = 'Secondary',sample_num = 2){
  data_singlestage <- FetchData(
    adata,layer = 'counts',
    vars = rownames(adata),
    cells = adata$sample[adata$group_stage==group_stage]
    ) %>% 
    t()
  gene_select <- data_singlestage %>% 
    apply(., 1, function(x){sum(x>0)}) %>%
    as.data.frame() %>%
    rename_all(~c('count')) %>%
    filter(count>sample_num) %>%
    rownames(.)
  data_singlestage <- data_singlestage[rownames(data_singlestage) %in% gene_select,]
  col_means <- apply(data_singlestage, 1, function(x){mean(x)}) %>% unlist()
  col_sds <- apply(data_singlestage, 1, function(x){sd(x)}) %>% unlist()

  cv_values <- (col_sds / col_means) %>%
    as.data.frame() %>% 
    rename_all(~c('CV_score')) %>% 
    mutate(stage = group_stage) %>%
    rownames_to_column('row_name')
  return(cv_values)
}
cv_calculate_all <- function(adata = scdata_pro,group_stage = NA,sample_num = 10){
  data_singlestage <- FetchData(
    adata,layer = 'counts',
    vars = rownames(adata)
  ) %>% 
    t() %>% as.data.frame() 
  gene_select <- data_singlestage %>% 
    rownames_to_column('Protein') %>% 
    rowwise() %>%
    filter(sum(c_across(-Protein) > 0) >= sample_num) %>%
    ungroup() %>%
    pull(Protein)
  data_singlestage <- data_singlestage[rownames(data_singlestage) %in% gene_select,]
  col_means <- apply(data_singlestage, 1, function(x){mean(x)}) %>% unlist()
  col_sds <- apply(data_singlestage, 1, function(x){sd(x)}) %>% unlist()

  cv_values <- (col_sds / col_means) %>%
    as.data.frame() %>% 
    rename_all(~c('CV_score')) %>% 
    mutate(stage = group_stage) %>%
    rownames_to_column('row_name')
  return(cv_values)
}



scdata_tran$group_cell %>% unique()
scdata_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = data_gene_protein$Genes)
protein_use <- data_gene_protein$Accessions[match(rownames(scdata_sub),data_gene_protein$Genes)]
rownames(scdata_sub)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
res_cv <- lapply(group_stage_list,FUN = function(x){
  cv_calculate(adata = scdata_sub,group_stage = x,sample_num = 2)
})
data_tran <- do.call(rbind,res_cv) %>%
  mutate(omics = 'Transcription') %>%
  rename('Genes' = 'row_name')

scdata_pro$group_cell %>% unique()
scdata_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = protein_use)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
res_cv <- lapply(group_stage_list,FUN = function(x){
  cv_calculate(adata = scdata_sub,group_stage = x,sample_num = 2)
})
data_pro <- do.call(rbind,res_cv) %>%
  mutate(omics = 'Protein') %>%
  rename('Accessions' = 'row_name') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  dplyr::select(-c('Accessions'))

data_plot <- rbind(data_pro,data_tran) %>% 
  group_by(Genes,stage) %>% 
  arrange(omics, .by_group = TRUE) %>%
  summarise(diff_cv = diff(CV_score), .groups = 'drop') %>% 
  pivot_wider(names_from = stage, values_from = diff_cv) %>%
  column_to_rownames('Genes') %>%
  rowwise() %>%
  na.omit() %>% 
  dplyr::select(c("Secondary","Early antral","Antral","Preovulatory"))
pheatmap::pheatmap(
  mat = data_plot,
  show_rownames = FALSE,
  scale = 'row',
  cluster_rows = TRUE,
  cluster_cols = FALSE
)


scdata_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = data_gene_protein$Genes)
gene_unselect_tran <- scdata_sub@assays$RNA$data_tpm %>% rowSums() %>% 
  as.data.frame() %>% rename_all(~c('Exp')) %>% 
  arrange(desc(Exp)) %>% rownames_to_column('Genes') %>% 
  slice(1:as.integer(n()*0.2),as.integer(n()*0.8):n()) %>% 
  pull(Genes)
scdata_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte') %>% 
  subset(.,features = protein_use)
gene_unselect_pro <- scdata_sub@assays$RNA$counts %>% rowSums() %>% 
  as.data.frame() %>% rename_all(~c('Exp')) %>% 
  arrange(desc(Exp)) %>% rownames_to_column('Genes') %>% 
  slice(1:as.integer(n()*0.2),as.integer(n()*0.8):n()) %>% 
  pull(Genes)
gene_unselect <- c(gene_unselect_tran,gene_unselect_pro) %>% unique()

data_plot <- data_pro %>% 
  left_join(data_tran,by = c('Genes','stage')) %>%
  na.omit() %>% 
  rename(
    cv_protein = 'CV_score.x',
    cv_tran = 'CV_score.y'
  ) %>%
  dplyr::select(c('Genes','cv_protein','cv_tran')) %>% 
  mutate(
    cv_protein_lg = log2(cv_protein+1),
    cv_tran_lg = log2(cv_tran+1),
    group = case_when(
      (cv_protein_lg > 1.3) & (cv_tran_lg > 0.55) ~ 'Double High Variability',
      (cv_protein_lg > 1.3) & (cv_tran_lg < 0.55) ~ 'High Transcription Variability',
      (cv_protein_lg < 1.3) & (cv_tran_lg > 0.55) ~ 'High Protein Variability',
      (cv_protein_lg < 1.3) & (cv_tran_lg < 0.55) ~ 'Double Low Variability'
    )
  ) %>% 
  mutate(
    group = factor(group,levels = c(
      'Double High Variability','High Transcription Variability',
      'High Protein Variability','Double Low Variability'))
  )
write.csv(data_plot,'./fig2g_cv_Oocyte.csv')

ggplot(
  data = data_plot,
  aes(
    x = cv_protein_lg,
    y = cv_tran_lg,
    color = group
  )) +
  geom_point() +
  geom_hline(yintercept = 0.55,color = 'black',linetype = 'dashed',linewidth = 1.5) +
  geom_vline(xintercept = 1.3,color = 'black',linetype = 'dashed',linewidth = 1.5) +
  guides(
    color = guide_legend(
      title = 'Type of Coefficient of Variation',
      override.aes = list(size = 5))
    ) +
  scale_color_manual(values = paletteer_d("ggsci::default_nejm")) +
  labs(
    x = 'Protien Coefficient of Variation',
    y = 'Transcription Coefficient of Variation'
    ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'top',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
ggsave(
  filename = './p_cv_dotplot.png',
  width = 15,
  height = 12
  )


data_cv <- data_plot
for(group_use in data_cv$group %>% levels()){
  library(clusterProfiler)
  library(stringr)
  print(group_use)
  gene_use <- data_cv$Genes[data_cv$group == group_use] %>% 
    unique()
  gene_symbol <-
    bitr(
      gene_use,
      fromType = 'SYMBOL',
      toType = 'ENTREZID',
      OrgDb = "org.Mm.eg.db"
    )

  enrich.go <- enrichGO(
    gene = gene_symbol$ENTREZID,

    OrgDb = 'org.Mm.eg.db',
    keyType = 'ENTREZID',
    ont = 'BP',

    pAdjustMethod = 'fdr',
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = FALSE
  )
  file_name <- paste(
    './cv_go_',
    group_use,'.csv',
    sep = ''
    )
  write.csv(enrich.go@result,file_name)
  data_plot <- enrich.go@result %>% 
    arrange(desc(Count)) %>% 
    mutate(Description = str_wrap(Description,width = 50))  %>% 
    mutate(
      Description = factor(Description,levels = rev(Description))
    ) %>% 
    filter(pvalue<0.05) %>% 
    slice(1:15)
  data_plot$Description
  plot_go <- ggplot(data_plot, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
    geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
    scale_fill_gradient(low = 'white',high = '#78A5CAFF') +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = '', y = 'Count') +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        vjust = 0.5,
        size = 20,
        face = 'bold'
      ),
      axis.text = element_text(
        size = 16,
        face = 'bold',
        color = 'black'
      ),
      axis.text.y = element_text(
        vjust = 0.5,hjust = 1,face = 'bold',size = 16),
      axis.title.x = element_text(
        size = 16,
        face = 'bold',
        hjust = 0.5),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(face = 'bold',size = 16),
      legend.text = element_text(size = 16)
    )
  file_name <- paste(
    './cv_go_',
    group_use,'.png',
    sep = ''
  )
  ggsave(
    plot = plot_go,
    filename = file_name,
    bg = 'white',
    width = 10,
    height = 12
  )
}


color_use <- paletteer_d("ggsci::default_nejm")
pathway_all <- read.csv('./fig2g_cv.csv')
colnames(pathway_all)
data_cv <- data_plot

pathway_select <- pathway_all$double.high %>% 
  unique() %>% .[nchar(.)>0] %>% .[1:9]
group_use <- data_cv$group %>% unique() %>% levels() %>% .[1]
file_name <- paste(
  './cv_go_',
  group_use,'.csv',
  sep = ''
)
colnames(data_plot_use)
data_plot_use <- read.csv(file_name,row.names = 1) %>% 
  filter(
    Description %in% pathway_select
  ) %>% 
  arrange(Count) %>% 
  mutate(Description = str_wrap(Description,width = 30))  %>%
  mutate(
    Description = factor(Description,levels = rev(Description))
  )
ggplot(data_plot_use, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot_use,
    aes(x = 0,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 90
  ) +
  scale_fill_gradient(low = 'white',high = color_use[1]) +
  scale_x_continuous(expand = c(0,0)) +
  coord_flip() +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.title = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16),
    legend.text = element_text(size = 16)
  )
pathway_select <- pathway_all$high.protein %>% 
  unique() %>% .[nchar(.)>0] %>% .[1:9]
group_use <- data_cv$group %>% levels() %>% .[3]
file_name <- paste(
  './cv_go_',
  group_use,'.csv',
  sep = ''
)
colnames(data_plot_use)
data_plot_use <- read.csv(file_name,row.names = 1) %>% 
  filter(
    Description %in% pathway_select
  ) %>% 
  arrange(Count) %>% 
  mutate(Description = str_wrap(Description,width = 40))  %>%
  mutate(
    Description = factor(Description,levels = rev(Description))
  )
ggplot(data_plot_use, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot_use,
    aes(x = 0,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 90
  ) +
  scale_fill_gradient(
    low = 'white',
    high = color_use[2],
    breaks = c(
      min(-log10(data_plot_use$pvalue)) %>% ceiling(), 
      max(-log10(data_plot_use$pvalue)) %>% floor()/2 + 
        min(-log10(data_plot_use$pvalue)) %>% ceiling()/2,
      max(-log10(data_plot_use$pvalue)) %>% floor()
      )
    ) +
  scale_x_continuous(expand = c(0,0)) +
  coord_flip() +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.title = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16),
    legend.text = element_text(size = 16)
  )

pathway_select <- pathway_all$high.transcriptome %>% 
  unique() %>% .[nchar(.)>0] %>% .[1:9]
group_use <- data_cv$group %>% levels() %>% .[2]
file_name <- paste(
  './cv_go_',
  group_use,'.csv',
  sep = ''
)
colnames(data_plot_use)
data_plot_use <- read.csv(file_name,row.names = 1) %>% 
  filter(
    Description %in% pathway_select
  ) %>% 
  arrange(Count) %>% 
  mutate(Description = str_wrap(Description,width = 40))  %>%
  mutate(
    Description = factor(Description,levels = Description)
  )
ggplot(data_plot_use, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot_use,
    aes(x = 0,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 90
  ) +
  scale_fill_gradient(
    low = 'white',
    high = color_use[3],
    breaks = c(
      min(-log10(data_plot_use$pvalue)) %>% ceiling(), 
      max(-log10(data_plot_use$pvalue)) %>% floor()/2 + 
        min(-log10(data_plot_use$pvalue)) %>% ceiling()/2,
      max(-log10(data_plot_use$pvalue)) %>% floor()
    )
  ) +
  scale_x_continuous(expand = c(0,0),position = 'top') +
  coord_flip() +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.title = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16,angle = 90,hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.text.position = 'left',
    legend.title.position = 'right',
    legend.position = 'left'
  )

pathway_select <- pathway_all$double.low %>% 
  unique() %>% .[nchar(.)>0] %>% .[1:9]
group_use <- data_cv$group %>% levels() %>% .[4]
file_name <- paste(
  './cv_go_',
  group_use,'.csv',
  sep = ''
)
colnames(data_plot_use)
data_plot_use <- read.csv(file_name,row.names = 1) %>% 
  filter(
    Description %in% pathway_select
  ) %>% 
  arrange(Count) %>% 
  mutate(Description = str_wrap(Description,width = 40))  %>%
  mutate(
    Description = factor(Description,levels = Description)
  )
ggplot(data_plot_use, aes(x = Count, y = Description,fill = -log10(pvalue),)) +
  geom_bar(stat = "identity",color = 'black', show.legend = TRUE,linewidth = 0.3,width = 0.8) +
  geom_text(
    data = data_plot_use,
    aes(x = 0,y = Description,label = Description),
    hjust = 0,vjust = 0.5,size = 6,angle = 90
  ) +
  scale_fill_gradient(
    low = 'white',
    high = color_use[4],
    breaks = c(
      min(-log10(data_plot_use$pvalue)) %>% ceiling(), 
      max(-log10(data_plot_use$pvalue)) %>% floor()/2 + 
        min(-log10(data_plot_use$pvalue)) %>% ceiling()/2,
      max(-log10(data_plot_use$pvalue)) %>% floor()
    )
  ) +
  scale_x_continuous(expand = c(0,0),position = 'top') +
  coord_flip() +
  labs(x = '', y = 'Count') +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      size = 20,
      face = 'bold'
    ),
    axis.text = element_text(
      size = 16,
      face = 'bold',
      color = 'black'
    ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 16,
      face = 'bold',
      hjust = 0.5),
    axis.title = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(face = 'bold',size = 16,angle = 90,hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.text.position = 'left',
    legend.title.position = 'right',
    legend.position = 'left'
  )



scdata_tran$group_cell %>% unique()
scdata_sub <- subset(scdata_tran,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = data_gene_protein$Genes)
protein_use <- data_gene_protein$Accessions[match(rownames(scdata_sub),data_gene_protein$Genes)]
rownames(scdata_sub)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
res_cv <- lapply(group_stage_list,FUN = function(x){
  cv_calculate(adata = scdata_sub,group_stage = x,sample_num = 2)
})
data_tran <- do.call(rbind,res_cv) %>%
  mutate(omics = 'Transcription') %>%
  rename('Genes' = 'row_name')

scdata_pro$group_cell %>% unique()
scdata_sub <- subset(scdata_pro,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = protein_use)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
res_cv <- lapply(group_stage_list,FUN = function(x){
  cv_calculate(adata = scdata_sub,group_stage = x,sample_num = 2)
})
data_pro <- do.call(rbind,res_cv) %>%
  mutate(omics = 'Protein') %>%
  rename('Accessions' = 'row_name') %>%
  left_join(data_gene_protein,by = 'Accessions') %>%
  dplyr::select(-c('Accessions'))

data_plot <- rbind(data_pro,data_tran) %>% 
  group_by(Genes,stage) %>% 
  arrange(omics, .by_group = TRUE) %>%
  summarise(diff_cv = diff(CV_score), .groups = 'drop') %>% 
  pivot_wider(names_from = stage, values_from = diff_cv) %>%
  column_to_rownames('Genes') %>%
  rowwise() %>%
  na.omit() %>% 
  dplyr::select(c("Secondary","Early antral","Antral","Preovulatory"))

scdata_sub <- subset(scdata_tran,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = data_gene_protein$Genes)
gene_unselect_tran <- scdata_sub@assays$RNA$data_tpm %>% rowSums() %>% 
  as.data.frame() %>% rename_all(~c('Exp')) %>% 
  arrange(desc(Exp)) %>% rownames_to_column('Genes') %>% 
  slice(1:as.integer(n()*0.2),as.integer(n()*0.8):n()) %>% 
  pull(Genes)
scdata_sub <- subset(scdata_pro,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = protein_use)
gene_unselect_pro <- scdata_sub@assays$RNA$counts %>% rowSums() %>% 
  as.data.frame() %>% rename_all(~c('Exp')) %>% 
  arrange(desc(Exp)) %>% rownames_to_column('Genes') %>% 
  slice(1:as.integer(n()*0.2),as.integer(n()*0.8):n()) %>% 
  pull(Genes)
gene_unselect <- c(gene_unselect_tran,gene_unselect_pro) %>% unique()

data_plot <- data_pro %>% 
  left_join(data_tran,by = c('Genes','stage')) %>%
  na.omit() %>% 
  rename(
    cv_protein = 'CV_score.x',
    cv_tran = 'CV_score.y'
  ) %>%
  dplyr::select(c('Genes','cv_protein','cv_tran')) %>% 
  mutate(
    cv_protein_lg = log2(cv_protein+1),
    cv_tran_lg = log2(cv_tran+1),
    group = case_when(
      (cv_protein_lg > 1.3) & (cv_tran_lg > 1.3) ~ 'Double High Variability',
      (cv_protein_lg > 1.3) & (cv_tran_lg < 1.3) ~ 'High Transcription Variability',
      (cv_protein_lg < 1.3) & (cv_tran_lg > 1.3) ~ 'High Protein Variability',
      (cv_protein_lg < 1.3) & (cv_tran_lg < 1.3) ~ 'Double Low Variability'
    )
  ) %>% 
  mutate(
    group = factor(group,levels = c(
      'Double High Variability','High Transcription Variability',
      'High Protein Variability','Double Low Variability'))
  )
data_plot$group %>% table()
p1 <- ggplot(data_plot, aes(cv_tran_lg)) +
  geom_density(fill = "#d62728", alpha = 0.5) +

  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )+
  coord_flip()
p2 <- ggplot(data_plot, aes(cv_protein_lg)) +
  geom_density(fill = "#d62728", alpha = 0.5) +

  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = 'Density') +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    axis.text = element_blank(),
    axis.line = element_line(color = 'black')
  )
write.csv(data_plot,'./data_fig2g_GC.csv')
density_data <- density(data_plot$cv_tran_lg)




p3 <- ggplot(
  data = data_plot,
  aes(
    x = cv_protein_lg,
    y = cv_tran_lg,
    color = group
  )) +
  geom_point() +
  geom_hline(yintercept = 1.3,color = 'black',linetype = 'dashed',linewidth = 1.5) +
  geom_vline(xintercept = 1.3,color = 'black',linetype = 'dashed',linewidth = 1.5) +
  guides(
    color = guide_legend(
      title = 'Type of Coefficient of Variation',
      override.aes = list(size = 5))
  ) +
  scale_color_manual(values = paletteer_d("ggsci::default_nejm")) +
  labs(
    x = 'Protien Coefficient of Variation',
    y = 'Transcription Coefficient of Variation'
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  )
p_res <- p3 %>%
  insert_top(p2,height = 0.2) %>% 
  insert_right(p1,0.2)
p_res

data_plot <- rbind(data_pro,data_tran) %>% 
  group_by(Genes,stage) %>% 
  arrange(omics, .by_group = TRUE) %>%
  reframe(diff_cv = diff(CV_score), .groups = 'drop') %>% 
  na.omit()
tmp <- data_plot %>% group_by(Genes,stage) %>% 
  arrange(diff_cv,.by_group = TRUE)
data_plot %>% group_by(stage) %>% 
  summarise(median= median(diff_cv),mean = mean(diff_cv))
ggplot(data = data_plot,aes(x = diff_cv,fill = stage)) +
  geom_histogram(bins = 100,position = position_dodge(width = 0.8)) +
  theme_classic()
ggplot(data = data_plot,aes(x = diff_cv)) +
  geom_histogram(bins = 100) +
  facet_wrap(~stage,ncol = 1) +
  theme_classic()


data_pro %>% colnames()
data_plot <- data_pro
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  dplyr::reframe(
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup()

ggplot(data = data_plot_tmp, aes(x = group, y = Counts)) +
  geom_bar(
    stat = 'identity',
    position = 'identity',
    width = 1,
    alpha = 0.9,
    fill = my36colors[1]
  ) +
  geom_bar(
    data = data_plot_tmp[data_plot_tmp$Counts == max(data_plot_tmp$Counts), ],
    aes(x = group, y = Counts),
    stat = 'identity',
    position = 'identity',
    width = 1,
    fill = 'grey30'
  ) +
  labs(x = 'coefficients of variation', y = 'Counts') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5, 0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20, margin = margin(b = 6, unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )

legend_label <- data_plot %>% 
  group_by(stage) %>%
  summarise(CV_score = median(CV_score)) %>% 
  mutate(stage = factor(stage,levels=group_stage_list)) %>% 
  arrange(stage) %>% 
  reframe(legend = paste(stage,' [median(CV)=',round(CV_score,2),']',sep = '')) %>% 
  pull(legend)
names(legend_label) <- group_stage_list
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  group_by(stage) %>% 
  dplyr::reframe(
    stage = unique(stage),
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup() %>% 
  mutate(legend = legend_label[stage]) %>% 
  mutate(legend = factor(legend,levels = legend_label))
data_plot_line <- data_plot_tmp %>% 
  group_by(stage) %>% 
  slice_max(order_by = Counts, n = 1)
ggplot(data = data_plot_tmp,aes(x = group,y = Counts,fill = legend)) +
  geom_bar(stat = 'identity',position = 'identity',width = 1,alpha = 0.9) +
  labs(x = 'coefficients of variation',y = 'Counts') +
  scale_fill_manual(values = my36colors) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5,0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
ggplot(
  data = data_plot,
  aes(x = stage, y = CV_score, fill = stage)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width=0.2,
    show.legend = FALSE,
    size = 0.5,
    outlier.alpha = 0
  ) +
  scale_fill_manual(values = color_group_stage_cell %>% as.character()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size = 20,face = 'bold'),
    axis.title.x = element_blank(),
    axis.title = element_text(size=25,face = 'bold'),
    legend.position = 'none'
  )

data_plot <- data_tran
legend_label <- data_plot %>% 
  summarise(CV_score = median(CV_score))
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  dplyr::reframe(
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup()

ggplot(data = data_plot_tmp, aes(x = group, y = Counts)) +
  geom_bar(
    stat = 'identity',
    position = 'identity',
    width = 1,
    alpha = 0.9,
    fill = my36colors[1]
  ) +
  geom_bar(
    data = data_plot_tmp[data_plot_tmp$Counts == max(data_plot_tmp$Counts), ],
    aes(x = group, y = Counts),
    stat = 'identity',
    position = 'identity',
    width = 1,
    fill = 'grey30'
  ) +
  labs(x = 'coefficients of variation', y = 'Counts') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5, 0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20, margin = margin(b = 6, unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
legend_label <- data_plot %>% 
  group_by(stage) %>%
  summarise(CV_score = median(CV_score)) %>% 
  mutate(stage = factor(stage,levels=group_stage_list)) %>% 
  arrange(stage) %>% 
  reframe(legend = paste(stage,' [median(CV)=',round(CV_score,2),']',sep = '')) %>% 
  pull(legend)
names(legend_label) <- group_stage_list
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  group_by(stage) %>% 
  dplyr::reframe(
    stage = unique(stage),
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup() %>% 
  mutate(legend = legend_label[stage]) %>% 
  mutate(legend = factor(legend,levels = legend_label))
data_plot_line <- data_plot_tmp %>% 
  group_by(stage) %>% 
  slice_max(order_by = Counts, n = 1)
ggplot(data = data_plot_tmp,aes(x = group,y = Counts,fill = legend)) +
  geom_bar(stat = 'identity',position = 'identity',width = 1,alpha = 0.9) +
  labs(x = 'coefficients of variation',y = 'Counts') +
  scale_fill_manual(values = my36colors) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5,0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
ggplot(
  data = data_plot,
  aes(x = stage, y = CV_score, fill = stage)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width=0.2,
    show.legend = FALSE,
    size = 0.5,
    outlier.alpha = 0
  ) +
  scale_fill_manual(values = color_group_stage_cell %>% as.character()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size = 20,face = 'bold'),
    axis.title.x = element_blank(),
    axis.title = element_text(size=25,face = 'bold'),
    legend.position = 'none'
  )
scdata$group_cell %>% unique()
scdata_sub <- subset(scdata,subset = group_cell == 'Granulosa') %>% 
  subset(.,features = gene_select)
group_stage_list <- scdata_sub$group_stage %>% unique()
scdata_sub$group_stage %>% table()
res_cv <- lapply(group_stage_list,FUN = function(x){
  cv_calculate(adata = scdata_sub,group_stage = x,sample_num = 2)
})
data_plot <- do.call(rbind,res_cv)

p1 <- ggplot(
  data = data_plot,
  aes(x = stage, y = CV_score, fill = stage)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width=0.2,
    show.legend = FALSE,
    size = 0.5,
    outlier.alpha = 0
  ) +
  scale_fill_manual(values = my36colors) +
  labs(title = 'Granulosa') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size = 20,face = 'bold'),
    axis.title.x = element_blank(),
    axis.title = element_text(size=25,face = 'bold'),
    legend.position = 'none'
  )
p1
legend_label <- data_plot %>% 
  summarise(CV_score = median(CV_score))
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01) 
data_plot_tmp <- data_plot %>% 
  dplyr::reframe(
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup()

ggplot(data = data_plot_tmp, aes(x = group, y = Counts)) +
  geom_bar(
    stat = 'identity',
    position = 'identity',
    width = 1,
    alpha = 0.9,
    fill = my36colors[2]
  ) +
  geom_bar(
    data = data_plot_tmp[data_plot_tmp$Counts == max(data_plot_tmp$Counts), ],
    aes(x = group, y = Counts),
    stat = 'identity',
    position = 'identity',
    width = 1,
    fill = 'grey30'
  ) +
  annotate(
    geom = 'text',
    x = 90,
    y = 480,
    label = paste('median(CV) = ', round(legend_label,2)),
    size = 8
  ) +
  labs(x = 'coefficients of variation', y = 'Counts') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.5, 0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20, margin = margin(b = 6, unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
legend_label <- data_plot %>% 
  group_by(stage) %>%
  summarise(CV_score = median(CV_score)) %>% 
  mutate(stage = factor(stage,levels=group_stage_list)) %>% 
  arrange(stage) %>% 
  reframe(legend = paste(stage,' [median(CV)=',round(CV_score,2),']',sep = '')) %>% 
  pull(legend)
names(legend_label) <- group_stage_list
breaks <- seq(
  from = data_plot$CV_score %>% min(), 
  to = data_plot$CV_score %>% max(), 
  by = 0.01)
data_plot_tmp <- data_plot %>% 
  group_by(stage) %>% 
  dplyr::reframe(
    stage = unique(stage),
    Counts = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% as.numeric(),
    group = cut(CV_score, breaks = breaks, include.lowest = TRUE, right = FALSE) %>% table() %>% names()
  ) %>% 
  mutate(
    group = factor(group,levels = group %>% unique())
  ) %>% 
  filter(Counts>0) %>% 
  ungroup() %>% 
  mutate(legend = legend_label[stage]) %>% 
  mutate(legend = factor(legend,levels = legend_label))
data_plot_line <- data_plot_tmp %>% 
  group_by(stage) %>% 
  slice_max(order_by = Counts, n = 1)
p2 <- ggplot(data = data_plot_tmp,aes(x = group,y = Counts,fill = legend)) +
  geom_bar(stat = 'identity',position = 'identity',width = 1,alpha = 0.9) +
  labs(x = 'coefficients of variation',y = 'Counts') +
  scale_fill_manual(values = my36colors) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.7,0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
p2

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Seurat)
library(ggpubr)
library(paletteer)
my36colors <-c(
  '#B9DDF1FF','#7EAED3FF','#5081AEFF','#2A5783FF',
  "#FFBEB2FF","#F8826BFF","#E33E43FF","#AE123AFF"
) 




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
gene_use <- protein_list %>% 
  filter(protein %in% c('Q8VIK3','B1B0V2')) %>% 
  pull(Genes)
gene_use

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

scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
scdata_tran@meta.data <- meta_data_tran
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))

histone_data<-read.csv('./histone_gene_use.csv',sep=",",header = T,check.names = F) %>% 
  pull(names(.)[2])

gene_protein_coding <-
  read.csv('./gene_protein_coding.csv', row.names = 1)
head(gene_protein_coding)
vlnplot_pro <- function(plottype = 'first',scdata_pro,gene_select,celltype = 'Oocyte',mycolor_use =NULL, num_col = NULL){
  if(is.null(mycolor_use)){
    mycolor_use <- if(celltype =='Oocyte'){my36colors[1:4]}else{my36colors[5:8]}
  }else{
    mycolor_use <- mycolor_use#rep(mycolor_use,4)
  }
  gene_select <- gene_select[gene_select%in% rownames(scdata_pro)]
  if(is.null(num_col)){
    if(length(gene_select) <= 6){
      num_col <- length(gene_select)
    }
  }
  data_plot <- FetchData(object = scdata_pro,layer = 'count',clean = 'none',vars = c(gene_select,'group_stage_cell')) %>% 
    rownames_to_column('sample') %>% 
    reshape2::melt(id.var = c('sample','group_stage_cell'),
                   variable.name = 'gene',
                   value.name = 'expression') %>% 
    mutate(
      group = ifelse(substr(sample,2,2) == 'O','Oocyte','Granulosa'),
      gene = gene %>% as.character(),
      group_stage_cell = factor(group_stage_cell,levels = c(
        "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte", "Preovulatory: Oocyte",
        "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
      ))
    ) %>%
    as.data.frame() %>% 
    filter(group %in% celltype) %>% 
    mutate(gene = factor(gene,levels = gene_select))

  if(plottype == 'first'){
    p_pro <- ggplot(data = data_plot,mapping = aes(x = group_stage_cell,y = expression,fill = group_stage_cell)) +
      geom_boxplot(width = 0.8,show.legend = FALSE,outliers = TRUE) +
      geom_point(position = position_jitter(width = 0.2),alpha = 0.7,show.legend = FALSE) +
      facet_wrap(~gene,ncol = num_col,scales = 'free',strip.position = 'left') +
      scale_fill_manual(values = mycolor_use) +
      labs(x = '',y = '') +
      theme_classic() +
      theme(
        axis.text = element_blank(),
        axis.title = element_text(size = 20),
        legend.position = 'right',
        strip.text = element_text(size = 16),
        strip.placement = 'outside',
        strip.background = element_blank()
      )
  }else{
    p_pro <- ggboxplot(
      data_plot, x="group_stage_cell", y="expression",
      color ="group_stage_cell",#fill = 'group_stage_cell',
      width = 0.6,
      palette = mycolor_use,#paletteer_d("ggsci::nrc_npg"),
      add = "jitter",
      xlab = F,  bxp.errorbar=T,
      bxp.errorbar.width=0.6,
      size=1, outlier.shape=NA,legend = "none") +
      stat_compare_means(size = 5) +
      facet_wrap(~gene,ncol = num_col,scales = 'free',strip.position = 'left') +
      theme(
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1.2, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 20,face = 'italic',vjust = 0.7),
        strip.placement = 'outside'
      )
  }
  return(p_pro)
}


vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Granulosa',
  gene_select = gene_use,
  num_col = 6)
vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_use,
  num_col = 6)
ggsave('./protein_exp_Oocyte.png',width = 8,height = 7)
ggsave('./protein_exp_Oocyte.pdf',width = 24,height = 7)



gene_select <- c(
  'Wtap','Npm1','Bod1','Saraf','Cks2','Spry4','Gdf9','Dnmt1','Dnmt3a'
)
mycolor_use <- c(
  rep('#00A087FF',4),rep('#4DBBD5FF',4),rep('#3C5488FF',4)
)
p1 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select[1:2],
  mycolor_use = mycolor_use[1:4],
  num_col = 2)
p1
p2 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select[3:4],
  mycolor_use = mycolor_use[5:8],
  num_col = 2)
p2
p3 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select[5:9],
  mycolor_use = mycolor_use[9:12],
  num_col = 5)
p3
cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(2,2,5))
ggsave('./genes_pro_Oocyte.png',width = 36,height = 7)
ggsave('./genes_pro_Oocyte.pdf',width = 36,height = 7)
p1 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_tran,
  celltype = 'Oocyte',
  gene_select = gene_select[1:2],
  mycolor_use = mycolor_use[1:4],
  num_col = 2)
p1
p2 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_tran,
  celltype = 'Oocyte',
  gene_select = gene_select[3:4],
  mycolor_use = mycolor_use[5:8],
  num_col = 2)
p2
p3 <- vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_tran,
  celltype = 'Oocyte',
  gene_select = gene_select[5:9],
  mycolor_use = mycolor_use[9:12],
  num_col = 5)
p3
cowplot::plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(2,2,5))
ggsave('./genes_trans_Oocyte.png',width = 36,height = 7)
ggsave('./genes_trans_Oocyte.pdf',width = 36,height = 7)
