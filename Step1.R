library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(tibble)
library(eulerr)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
data <- read.csv('./other research.csv')


data_zhang <- data[,1] %>% .[nchar(.)>0] %>% 
  str_remove('2::') %>% str_split(';',2) %>% 
  as.data.frame() %>% .[1,] %>% as.character() %>%
  unique()
sample_select <- c(
  "SGC1","SGC2","SGC3","SGC4","SGC5","SGC6","SGC7","SGC8","SGC9",
  "SGC10","SGC11","SGC12","SO1","SO2","SO3","SO4","SO5","SO6","SO7",
  "SO8","SO9","SO10","SO11","SO12","AGC2","AGC3","AGC4","AGC5","AGC6","AGC7","AGC8","AGC9",
  "AGC10","AGC11","AGC12","AO2","AO3","AO4","AO5","AO6","AO7",
  "AO8","AO9","AO10","AO11","AO12"
)
data_protein <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  column_to_rownames('ProteinGroups') %>% 
  select(all_of(sample_select)) %>% 
  rownames_to_column('protein') %>% 
  rowwise() %>%
  mutate(sum = sum(c_across(-protein))) %>% 
  select(protein,sum) %>% 
  filter(sum>0) %>% 
  pull(protein)
data_plot <- list(
  `Zhang et al.` = data_zhang,
  `This study` = data_protein
)
library(venn)
names(data_plot)
venn::venn(
  x = data_plot,
  zcolor = 'style',
  borders = FALSE,
  ggplot = TRUE,
  box = FALSE,
  ilcs = 1,
  sncs = 1
) +
  theme(
    panel.border = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = 'pt')
  )

fig_1 <- plot(euler(data_plot),
     quantities = list(col="black", font=2, cex=1.5),
     labels = list(col="black",font=3,cex=2),
     edges = list(col=my36colors[1:4],lwd=5,lty=1),
     main = list(label=c("Secondary/Antral"),cex=2)
)
fig_1

data_li <- data[,2] %>% .[nchar(.)>0] %>% 
  str_remove('2::') %>% str_split(';',2) %>% 
  as.data.frame() %>% .[1,] %>% as.character() %>%
  unique()
data_sun <- data[,3] %>% .[nchar(.)>0] %>% 
  str_remove('2::') %>% str_split(';',2) %>% 
  as.data.frame() %>% .[1,] %>% as.character() %>%
  unique()
sample_select <- c(
  "PO1","PO2","PO3","PO4","PO5",
  "PO6","PO7","PO8","PO9"#,"PGC10","PO10"
)
data_protein <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  column_to_rownames('ProteinGroups') %>% 
  select(all_of(sample_select)) %>% 
  rownames_to_column('protein') %>% 
  rowwise() %>%
  mutate(sum = sum(c_across(-protein))) %>% 
  select(protein,sum) %>% 
  filter(sum>0) %>% 
  pull(protein)
data_plot <- list(
  `Li et al.` = data_li,
  `This study` = data_protein
)

fig_2 <- plot(euler(data_plot),
              quantities = list(col="black", font=2, cex=1.5),
              labels = list(col="black",font=3,cex=2),
              edges = list(col=my36colors[1:4],lwd=5,lty=1),
              main = list(label=c("Ocyte(Preovulatory)"),cex=2)
)
fig_2
data_plot <- list(
  `Sun et al.` = data_sun,
  `This study` = data_protein
)

fig_3 <- plot(euler(data_plot),
              quantities = list(col="black", font=2, cex=1.5),
              labels = list(col="black",font=3,cex=2),
              edges = list(col=my36colors[1:4],lwd=5,lty=1),
              main = list(label=c("Ocyte(Preovulatory)"),cex=2)
)
fig_3
cowplot::plot_grid(fig_1,fig_2,fig_3,ncol = 3)

ggsave(
  filename = './fig1c_venn.png',
  width = 24,height = 8
)
ggsave(
  filename = './fig1c_venn.pdf',
  width = 24,height = 8
)

data_plot <- list(
  `Li et al.` = data_li,
  sun.et.al = data_sun,
  `This study` = data_protein
)

fig_4 <- plot(euler(data_plot),
              quantities = list(col="black", font=2, cex=1.5),
              labels = list(col="black",font=3,cex=2),
              edges = list(col=my36colors[1:4],lwd=5,lty=1),
              main = list(label=c("Ocyte(Preovulatory)"),cex=2)
)
fig_4
cairo_pdf(
  filename = './fig1c_venn_20260103.png',
  width = 12,height = 12
)
fig_4

cairo_pdf(
  filename = './fig1c_venn_20260103.pdf',
  width = 12,height = 12
)
print(fig_4) 
dev.off()



setwd('./protein/')
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#7698b3","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#d6616b","#a55194","#ce6dbd","#e7ba52",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252"
) 
color_gc <- '#477BD1';color_oo <- '#E25247'
color_group_stage_cell <- c(
  '#2CA02C','#0088D4','#17BECF','#006BA4',
  '#F8B620','#FF800E','#D14B70','#ED444A'
)
names(color_group_stage_cell) <- c(
  "Secondary: Granulosa", "Early antral: Granulosa",
  "Antral: Granulosa", "Preovulatory: Granulosa",
  "Secondary: Oocyte","Early antral: Oocyte",
  "Antral: Oocyte", "Preovulatory: Oocyte"
)
data_protein <- read.csv('./data_protein_diann_mbr.csv')
data_clean <- read.csv('./data_protein_diann_mbr.csv',row.names = 1)

sample_level <- c(
  "SGC1","SGC2","SGC3","SGC4","SGC5","SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "SO1","SO2","SO3","SO4","SO5","SO6","SO7","SO8","SO9","SO10","SO11","SO12",
  "EGC2","EGC3","EGC4","EGC5","EGC6","EGC7",#"EGC1","EGC9","EGC10","EGC11","EGC12","EGC13",
  "EO2","EO3","EO4","EO5","EO6","EO7",#"EO1","EO9","EO10","EO11","EO12","EO13",
  "AGC2","AGC3","AGC4","AGC5","AGC7","AGC8","AGC9","AGC10","AGC11","AGC12",#"AGC6",
  "AO2","AO3","AO4","AO5","AO7","AO8","AO9","AO10","AO11","AO12",#"AO6",
  "PGC1","PGC2","PGC3","PGC4","PGC5","PGC6","PGC7","PGC8","PGC9",#"PGC10",
  "PO1","PO2","PO3","PO4","PO5","PO6","PO7","PO8","PO9"
)
data <- data_clean[,c(colnames(data_clean)[c(1, 2)],sample_level)]
group_stage_name <- c('Secondary','Early antral','Antral','Preovulatory')
names(group_stage_name) <- c('S','E','A','P')
group_cell_name <- c('Oocyte','Granulosa')
names(group_cell_name) <- c('O','G')
sample_meta_data <- 
  lapply(colnames(data)[-c(1, 2)], function(sample) {
    sum(!is.na(data[, sample]))
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1, 2)]) %>% 
  dplyr::select(sample) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))

protein_num <-
  lapply(colnames(data)[-c(1, 2)], function(sample) {
    sum(!is.na(data[, sample]))
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1, 2)]) %>% 
  dplyr::select(sample,num_pro) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))
stage_cell <- c('all','all_Granulosa','all_Oocyte',unique(protein_num$group_stage_cell))
data_plot <- lapply(stage_cell, function(stage_cell){
  if(stage_cell == 'all'){
    sample <- protein_num$sample
  }else if(stage_cell == 'all_Granulosa'){
    sample <- protein_num$sample[grepl('Granulosa',protein_num$group_stage_cell)] %>% 
      as.character()
  }else if(stage_cell == 'all_Oocyte'){
    sample <- protein_num$sample[grepl('Oocyte',protein_num$group_stage_cell)] %>% 
      as.character()
  }else{
    sample <- protein_num$sample[protein_num$group_stage_cell==stage_cell] %>% 
      as.character()
  }
  data_tmp <- data[,colnames(data) %in% sample]
  data_tmp[is.na(data_tmp)] <- 0
  pro_num <- sum(rowSums(data_tmp)>0)
  return(pro_num)
}) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(
    group_stage_cell = stage_cell,
    sample = stage_cell
  ) %>% 
  dplyr::select(group_stage_cell,sample,num_pro) %>% 
  rbind(protein_num[,c('group_stage_cell','sample','num_pro')])
data_plot$group_stage_cell[data_plot$sample %in% c('all_Granulosa', 'all_Oocyte')] <-
  'all'
data_plot$sample %>% unique()
sample_level <- c(
  "rep1","rep2","rep3","rep4","all","all_Granulosa","all_Oocyte",
  "rep7","rep8","rep9", "rep10",
  "Secondary: Granulosa","SGC1","SGC2","SGC3","SGC4","SGC5",
  "SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "Secondary: Oocyte","SO1","SO2","SO3","SO4","SO5","SO6",
  "SO7","SO8","SO9","SO10","SO11","SO12",
  "Early antral: Granulosa","EGC2","EGC3","EGC4","EGC5",
  "EGC6","EGC7",
  "Early antral: Oocyte","EO2","EO3","EO4","EO5","EO6","EO7",
  "Antral: Granulosa","AGC2","AGC3","AGC4","AGC5","AGC7",
  "AGC8","AGC9","AGC10","AGC11","AGC12",
  "Antral: Oocyte","AO2","AO3","AO4","AO5","AO7","AO8",
  "AO9","AO10","AO11","AO12",
  "Preovulatory: Granulosa","PGC1","PGC2","PGC3","PGC4",
  "PGC5","PGC6","PGC7","PGC8","PGC9",
  "Preovulatory: Oocyte","PO1","PO2","PO3","PO4","PO5",
  "PO6","PO7","PO8","PO9"
)
color_grey_sample <- c(
  "all","all_Granulosa","all_Oocyte",
  "Secondary: Granulosa","Secondary: Oocyte",
  "Early antral: Granulosa","Early antral: Oocyte", 
  "Antral: Granulosa", "Antral: Oocyte",
  "Preovulatory: Granulosa","Preovulatory: Oocyte")
color_use <- lapply(sample_level, function(sample){
  if(sample %in% color_grey_sample){
    return('lightslategrey')
  }else{
    sample_group <- sample_meta_data$group_stage_cell[sample_meta_data$sample == sample]
    color_tmp <- color_group_stage_cell[sample_group]
  }
  return(color_tmp)
}) %>% unlist() %>% as.character()
color_use[c(2,3)] <- c(color_gc,color_oo)
data_plot$sample <- factor(data_plot$sample,levels = sample_level)
data_plot$group_stage_cell <- factor(data_plot$group_stage_cell,levels = c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
))

label <- c('All','SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell],
    label_text = factor(label_text,levels = label %>% as.character())
  )
write.csv(data_plot,'./data_plot_protein-number.csv')

ggplot(
  data = data_plot,
  mapping = aes(x = label_text, y = num_pro, fill = sample)
) +
  geom_bar(
    stat = 'identity',
    position = position_dodge2(preserve = "single",padding = 0),
    color = 'white') +
  geom_hline(yintercept = 1500,linewidth = 1.5,color = 'grey80',linetype = 'longdash') +
  geom_hline(yintercept = 2500,linewidth = 1.5,color = 'grey80',linetype = 'longdash') +
  scale_y_continuous(
    labels = c(expression(italic(0)),
               expression(3%*%10^3),
               expression(6%*%10^3),
               expression(9%*%10^3),
               expression(12%*%10^3),
               expression(13%*%10^3)),
    expand = c(0,0),
    breaks = c(0,3%*%10^3,6%*%10^3,9%*%10^3,12%*%10^3,13%*%10^3),
    limits = c(0,13.5*10^3)
  ) +
  scale_fill_manual(values = color_use) +
  labs(y = 'No. of identified proteins') +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = 'none',
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      size = 18,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(
      size = 18, margin = margin(b = 6, unit = 'pt')
    ),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1b_protein.png',
  width = 12,height = 4
)
ggsave(
  filename = './fig1b_protein.pdf',
  width = 12,height = 4
)


protein_num <-
  lapply(colnames(data)[-c(1, 2)], function(sample) {
    sum(data[!is.na(data[, sample]),sample])
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1, 2)]) %>% 
  dplyr::select(sample,num_pro) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))
stage_cell <- c('all','all_Granulosa','all_Oocyte',unique(protein_num$group_stage_cell))
data_plot <- lapply(stage_cell, function(stage_cell){
  if(stage_cell == 'all'){
    sample <- protein_num$sample
  }else if(stage_cell == 'all_Granulosa'){
    sample <- protein_num$sample[grepl('Granulosa',protein_num$group_stage_cell)] %>% 
      as.character()
  }else if(stage_cell == 'all_Oocyte'){
    sample <- protein_num$sample[grepl('Oocyte',protein_num$group_stage_cell)] %>% 
      as.character()
  }else{
    sample <- protein_num$sample[protein_num$group_stage_cell==stage_cell] %>% 
      as.character()
  }
  data_tmp <- data[,colnames(data) %in% sample]
  data_tmp[is.na(data_tmp)] <- 0
  pro_num <- sum(rowSums(data_tmp)/length(sample))
  return(pro_num)
}) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(
    group_stage_cell = stage_cell,
    sample = stage_cell
  ) %>% 
  dplyr::select(group_stage_cell,sample,num_pro) %>% 
  rbind(protein_num[,c('group_stage_cell','sample','num_pro')])
data_plot$group_stage_cell[data_plot$sample %in% c('all_Granulosa', 'all_Oocyte')] <-
  'all'
data_plot$sample %>% unique()
sample_level <- c(
  "rep1","rep2","rep3","rep4","all","all_Granulosa","all_Oocyte",
  "rep7","rep8","rep9", "rep10",
  "Secondary: Granulosa","SGC1","SGC2","SGC3","SGC4","SGC5",
  "SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "Secondary: Oocyte","SO1","SO2","SO3","SO4","SO5","SO6","SO7",
  "SO8","SO9","SO10","SO11","SO12",
  "Early antral: Granulosa","EGC1","EGC2","EGC3","EGC4","EGC5",
  "EGC6","EGC7","EGC9","EGC10","EGC11","EGC12","EGC13",
  "Early antral: Oocyte","EO1","EO2","EO3","EO4","EO5","EO6","EO7",
  "EO9","EO10","EO11","EO12","EO13",
  "Antral: Granulosa","AGC2","AGC3","AGC4","AGC5","AGC6","AGC7",
  "AGC8","AGC9","AGC10","AGC11","AGC12",
  "Antral: Oocyte","AO2","AO3","AO4","AO5","AO6","AO7","AO8","AO9",
  "AO10","AO11","AO12",
  "Preovulatory: Granulosa","PGC1","PGC2","PGC3","PGC4","PGC5",
  "PGC6","PGC7","PGC8","PGC9","PGC10",
  "Preovulatory: Oocyte","PO1","PO2","PO3","PO4","PO5","PO6","PO7",
  "PO8","PO9"#,"PO10"#"AO1",
)
color_grey_sample <- c(
  "all","all_Granulosa","all_Oocyte",
  "Secondary: Granulosa","Secondary: Oocyte",
  "Early antral: Granulosa","Early antral: Oocyte", 
  "Antral: Granulosa", "Antral: Oocyte",
  "Preovulatory: Granulosa","Preovulatory: Oocyte")
color_use <- lapply(sample_level, function(sample){
  if(sample %in% color_grey_sample){
    return('lightslategrey')
  }else{
    sample_group <- sample_meta_data$group_stage_cell[sample_meta_data$sample == sample]
    color_tmp <- color_group_stage_cell[sample_group]
  }
  return(color_tmp)
}) %>% unlist() %>% as.character()
color_use[c(2,3)] <- c(color_gc,color_oo)
data_plot$sample <- factor(data_plot$sample,levels = sample_level)
data_plot$group_stage_cell <- factor(data_plot$group_stage_cell,levels = c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
))

label <- c('All','SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell],
    label_text = factor(label_text,levels = label %>% as.character())
  )
write.csv(data_plot,'./data_plot_protein-counts.csv')

ggplot(
  data = data_plot,
  mapping = aes(x = label_text, y = num_pro, fill = sample)
) +
  geom_bar(
    stat = 'identity',
    position = position_dodge2(preserve = "single",padding = 0),
    color = 'white') +
  scale_y_continuous(
    labels = c(expression(italic(0)),
               expression(6%*%10^7),
               expression(9%*%10^7),
               expression(12%*%10^7),
               expression(15%*%10^7),
               expression(18%*%10^7)),
    expand = c(0,0),
    breaks = c(0,6*10^7,9*10^7,12*10^7,15*10^7,18*10^7),
    limits = c(0,20*10^7)
  ) +
  scale_fill_manual(values = color_use) +
  labs(y = 'No. of identified proteins') +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = 'none',
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      size = 18,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(
      size = 18, margin = margin(b = 6, unit = 'pt')
    ),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1b_protein_counts.png',
  width = 18,height = 7
)
ggsave(
  filename = './fig1b_protein_counts.pdf',
  width = 18,height = 7
)

data_clean[is.na(data_clean)] <- 0
sample_select <- lapply(3:ncol(data_clean), function(sample_i){
  sample_data <- data_clean[,sample_i]
  if(sum(sample_data>0)>2000)return(colnames(data_clean)[sample_i])
}) %>% unlist()
data_clean <- data_clean %>% 
  dplyr::select(all_of(c(colnames(data_clean)[c(1,2)],sample_select)))
meta_data <- sample_meta_data
meta_data$group_stage_cell <- paste(meta_data$group_stage,meta_data$group_cell,sep = ': ')
group_stage_cell <- meta_data$group_stage_cell %>% unique()
data_plot <- lapply(group_stage_cell, function(group_use) {
  sample <-
    meta_data$sample[meta_data$group_stage_cell %in% group_use]
  data_tmp <- data_clean[, colnames(data_clean) %in% sample]
  protein_all <- apply(data_tmp, 2, function(sample) {
    data_clean$PG.ProteinGroups[sample > 0] %>% .[!is.na(.)]
  }) %>% unlist()
  protein_all_not_dup <- protein_all %>% unique()
  protein_all_dup <-
    protein_all %>% table() %>% as.data.frame() %>%
    rename_all( ~ c('protein', 'num')) %>% as.data.frame()
  data_res <- table(protein_all_dup$num) %>% as.data.frame() %>%
    rename_all( ~ c('sample_num', 'gene_num')) %>%
    mutate(group_stage_cell = group_use,
           duplicated_rate = gene_num / sum(gene_num) * 100)
  return(data_res)
}) %>% do.call(rbind, .) %>%
  mutate(
    group_stage_cell = factor(
      group_stage_cell,
      levels = c(
        "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
        "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
      )
    ),
    sample_num = factor(
      sample_num %>% as.numeric(),
      levels = 1:(sample_num %>% as.numeric() %>% max())# %>% rev()
    )
  ) %>%
  group_by(group_stage_cell) %>%
  mutate(
    tmp = sample_num %>% as.character() %>% as.numeric(),
    group_color = case_when(
      (tmp > max(tmp)*0.75) ~ '75%~100%',
      (tmp > max(tmp)*0.5) ~ '50%~75%',
      (tmp > max(tmp)*0.25) ~ '25%~50%',
      TRUE ~ '0~25%'
    )
  ) %>% 
  group_by(group_stage_cell,group_color) %>% 
  mutate(
    duplicated_rate_sum = sum(duplicated_rate)
  ) %>% 
  select(-sample_num,-tmp,-duplicated_rate,-gene_num) %>% distinct() %>% 
  group_by(group_stage_cell) %>%
  arrange(group_color) %>%
  mutate(label_loc = cumsum(duplicated_rate_sum) - duplicated_rate_sum / 2) %>% 
  dplyr::rename(duplicated_rate = duplicated_rate_sum)

label <- c('SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  ungroup() %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell] %>% as.character(),
    label_text = factor(label_text,levels = label %>% as.character())
  )
ggplot(
  data = data_plot,
  aes(x = label_text,y = -duplicated_rate)
) +
  geom_bar(
    aes(group = rev(group_color),fill = group_color),
    stat = 'identity',color = 'black') +
  geom_text(
    aes(
      x = label_text,
      y = -label_loc,
      label = paste(round(duplicated_rate, 1),'%',sep = '')
    ),
    color = 'black',
    size = 4
  ) +
  scale_fill_manual(
    values = paletteer_c("grDevices::Oslo", 50) %>% .[c(5,10,15,30)]
  ) +
  scale_y_continuous(
    expand = c(0,0.1),
    breaks = seq(-100,0,25),
    label = seq(0,100,25)
  ) +
  scale_x_discrete(expand = c(0,0.5)) +
  labs(y = 'repeatability (%)') +
  theme_classic() +
  theme(
    plot.margin = margin(l = 25,t = 10),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(angle = 0,size = 18,hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 24),
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )



dir.create('./figure1_20250508',showWarnings = FALSE)
ggsave(
  filename = './fig1c_protein.png',
  width = 8,height = 8
)
ggsave(
  filename = './fig1c_protein.pdf',
  width = 8,height = 8
)

write.csv(data_plot,'./data_fig1c_protein.csv',row.names = FALSE)


library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
library(purrr)
library(stringr)
marker_tabel <- readxl::read_xlsx('./marker gene.xlsx')
data_copy_num_sorted <-
  read.csv(
    "./data_copy_num_sorted.csv",
    row.names = 1)
colnames(data_copy_num_sorted)
data_gene_protein <- data_copy_num_sorted[,c('ProteinGroups','Genes')] %>% 
  rename_all(~c('Accessions','Genes'))
scdata <- readRDS('./scdata_use.rds')
meta_data <- scdata@meta.data
colnames(data_copy_num_sorted)
colnames(meta_data)
data_plot <- data_copy_num_sorted %>%
  dplyr::select(
    c(
      'ProteinGroups','Genes',
      meta_data$sample %>% as.character()
    )) %>% 
  pivot_longer(
    cols = names(.)[-c(1,2)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  left_join(
    meta_data %>% dplyr::select(c(group_stage,group_cell,sample)),
    by = 'sample'
  ) %>% 
  group_by(group_cell,ProteinGroups) %>% #group_stage,
  summarise(
    Genes = Genes[1],
    copy_number = mean(copy_number),
    .groups = 'keep'
  ) %>% 
  filter(copy_number>0) %>% 
  group_by(group_cell) %>% #group_stage
  arrange(desc(copy_number),.by_group = TRUE) %>%
  mutate(
    rank_num = row_number(),
  ) %>% 
  arrange(rank_num) %>% 
  mutate(
    cumulative_sum = cumsum(copy_number)/sum(copy_number)
  )
colnames(data_plot)

data_plot_sub <- data_plot

data_summary <- data_plot %>% 
  mutate(
    group_rank = case_when(
      cumulative_sum <= 0.25 ~ 'Q1',
      cumulative_sum <= 0.5 ~ 'Q2',
      cumulative_sum <= 0.75 ~ 'Q3',
      cumulative_sum <= 1 ~ 'Q4',
    )
  ) %>% 
  group_by(group_cell,group_rank) %>% 
  summarise(
    Counts = n(),.groups = 'drop'
  ) %>% 
  arrange(group_rank, factor(group_cell, levels = c("Oocyte", "Granulosa"))) %>%
  group_by(group_rank) %>%
  reframe(
    label = paste0(group_rank, "[", paste(group_cell, Counts, sep = ": ", collapse = ";"), "]")
  ) %>%
  distinct() %>% 
  mutate(
    x_loc = 7500,
    y_loc = seq(0.05,1,by = 0.25)
  )

data_plot_sub_label <- data_plot_sub %>% 
  group_by(group_cell) %>%
  arrange(rank_num) %>%
  dplyr::slice(1:10) %>% 
  mutate(rank_num = 1:10) 
colnames(data_plot_sub)
ggplot() + 
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0,
      ymax = 0.25
    ),
    fill = "#438E80",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.25,
      ymax = 0.5
    ),
    fill = "#336C61",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.5,
      ymax = 0.75
    ),
    fill = "#E69D70",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +  
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.75,
      ymax = 1
    ),
    fill = "#B65A20",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data_summary,
    mapping = aes(x = 10500,y = y_loc,label = label),
    size = 5,hjust = 1
  ) +
  geom_point(
    data = data_plot_sub,
    mapping = aes(
      x = rank_num, 
      y = cumulative_sum,
      color = group_cell)
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = cumulative_sum), 
    color = 'red',size = 3
  ) +
  geom_text(
    data = data_plot_sub_label[data_plot_sub_label$group_cell == 'Oocyte',],
    mapping = aes(x = 400,y = rank_num*0.05,label = Genes),
    size = 5,hjust = 0,color = '#865254'
  ) +
  geom_text(
    data = data_plot_sub_label[data_plot_sub_label$group_cell == 'Granulosa',],
    mapping = aes(x = 2000,y = rank_num*0.05,label = Genes),
    size = 5,hjust = 0,color = '#4d627c'
  ) +
  guides(color = guide_legend(
    title = 'Cell Type',
    override.aes = list(size = 4)
  )) +
  scale_color_manual(values = c('#4d627c','#865254')) +
  labs(
    x = 'Abundance Rank', 
    y = 'Cumulative Abundance[%]', 
    title = 'Protein') +
  scale_x_continuous(
    expand = c(0,100),
    labels = c(0,2500,5000,7500,9500,10850),
    breaks = c(0,2500,5000,7500,9500,10850)
    ) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
ggsave(
  filename = './fig1f_pro.png',
  width = 8,height = 12
)
ggsave(
  filename = './fig1f_pro.pdf',
  width = 8,height = 12
)

data_plot_sub <- data_plot
marker_select <- marker_tabel %>% 
  pull(`marker gene`) %>% 
  unique() %>% str_to_title()
data_plot_sub_label <- data_plot_sub %>% 
  left_join(
    marker_tabel %>% 
      dplyr::rename(Genes = `marker gene`),
    by = 'Genes'
    ) %>% 
  filter(
    !is.na(`cell type`)
  )
colnames(data_plot_sub)
ggplot(data = data_plot_sub) + 
  geom_point(
    aes(x = rank_num, y = copy_number,color = group_cell)
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = copy_number), 
    color = 'black',size = 3
  ) +
  geom_text_repel(
    data = data_plot_sub_label,
    mapping = aes(
      x = rank_num, 
      y = copy_number,
      color = group_cell,
      label = Genes
    ),
    show.legend = FALSE,
    force = 100,segment.curvature = 0.8,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = Inf,
    size = 6
  )+
  guides(color = guide_legend(
    title = 'Cell Type',
    override.aes = list(size = 4)
    )) +
  scale_color_manual(values = c('#865254','#4d627c')) +
  labs(
    x = 'Protein Number', 
    y = 'Log10(copy number + 1)', 
    title = 'Protein') +
  scale_x_continuous(expand = c(0,100)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
data_plot_sub <- data_plot %>% 
  filter(group_cell == 'Oocyte')
marker_select <- marker_tabel %>% 
  pull(`marker gene`) %>% 
  unique() %>% str_to_title()
data_plot_sub_label <- data_plot_sub %>% 
  filter(
    Genes %in% marker_select
  )
ggplot(data = data_plot_sub) + 
  geom_point(
    aes(x = rank_num, y = copy_number), 
    color = '#C34F73'
    ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = copy_number), 
    color = 'black',size = 3
  ) +
  geom_text_repel(
    data = data_plot_sub_label,
    mapping = aes(
      x = rank_num, 
      y = copy_number,
      label = Genes
    ),force_pull = 2,nudge_x = 10,size = 6
  )+
  labs(x = 'Protein Number', y = 'Log10(copy number + 1)', title = 'Oocyte') +
  scale_x_continuous(expand = c(0,100)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
ggsave(
  filename = './fig1d_protein_Oocyte.png',
  width = 6,height = 10
)
data_plot_sub <- data_plot %>% 
  filter(group_cell == 'Granulosa')
marker_select <- marker_tabel %>% 
  filter(`cell type` == 'GC') %>% 
  pull(`marker gene`) %>% 
  unique() %>% str_to_title()
data_plot_sub_label <- data_plot_sub %>% 
  filter(
    Genes %in% marker_select
  )
ggplot(data = data_plot_sub) + 
  geom_point(
    aes(x = rank_num, y = copy_number), 
    color = '#0088D4'
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = copy_number), 
    color = 'black',size = 3
  ) +
  geom_text_repel(
    data = data_plot_sub_label,
    mapping = aes(
      x = rank_num, 
      y = copy_number,
      label = Genes
    ),force_pull = 2,nudge_x = 10,size = 6
  )+
  labs(x = 'Protein Number', y = 'Log10(copy number + 1)', title = 'Granulosa') +
  scale_x_continuous(expand = c(0,100)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
ggsave(
  filename = './fig1d_protein_Granulosa.png',
  width = 6,height = 10
)



library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
fun_go <- function(data_plot,group_label,group_cell = 'Oocyte'){
  print(group_label)
  genes <- data_plot %>% 
    filter(expression_group == group_label) %>% 
    pull(Genes)
  gene <- str_to_title(genes)

  gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO of",group_label),font.size = 20)
  file_name <- paste(
    "./GO_",
    group_cell,'_',group_label,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")

  data_go <- ego@result
  file_name <- paste(
    "./GO_",
    group_cell,'_',group_label,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
  return(TRUE)
}
group_labels <- data_plot_o$expression_group %>% unique()
lapply(group_labels, function(group_label) {
  fun_go(data_plot = data_plot_o,
         group_label = group_label,
         group_cell = 'Oocyte')
})
lapply(group_labels, function(group_label) {
  fun_go(data_plot = data_plot_g,
         group_label = group_label,
         group_cell = 'Granulosa')
})


data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
data_tf_list <- c(data_raw_1$Symbol,data_raw_2$tf) %>% unique

library(stringr)
library(msigdbr)
library(dplyr)
gene_list <- msigdbr(species = 'Mus musculus')
gene_list$gs_subcat %>% unique()
data_sm_pathway <- read.csv('./KEGG_Signal molecule.txt') %>% 
  pull(x)
data_sm_pathway %in% (gene_list$gs_name %>% unique())
colnames(gene_list)
gene_select <- gene_list %>% 
  filter(gs_name %in% data_sm_pathway) %>% 
  pull(gene_symbol) %>% unique()
data_sm_list <- gene_list %>% 
  filter(gs_name %in% data_sm_pathway) %>% 
  pull(gene_symbol) %>% unique()
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
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
    mutate(group = 'All'),
  data_plot_pro %>% 
    filter(group_cell == 'Oocyte',Genes %in% data_tf_list) %>% 
    mutate(group = 'TFs'),
  data_plot_pro %>% 
    filter(group_cell == 'Oocyte',Genes %in% data_sm_list) %>% 
    mutate(group = 'SMs')
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
    legend.position = c(0.1,0.9),
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
  filename = './fig1e_proteins_Oocyte.png',
  width = 6,height = 6
)
data_plot_pro_Granulosa <- data_plot_pro %>% 
  filter(group_cell == 'Granulosa')

data_plot <- rbind(
  data_plot_pro %>% 
    filter(group_cell == 'Granulosa') %>% 
    mutate(group = 'All'),
  data_plot_pro %>% 
    filter(group_cell == 'Granulosa',Genes %in% data_tf_list) %>% 
    mutate(group = 'TFs'),
  data_plot_pro %>% 
    filter(group_cell == 'Granulosa',Genes %in% data_sm_list) %>% 
    mutate(group = 'SMs')
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
    legend.position = c(0.1,0.9),
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
  filename = './fig1e_proteins_Granulosa.png',
  width = 6,height = 6
)



library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
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


scdata_tran <- readRDS('./scdata_tpm.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))



colnames(meta_data_tran)
colnames(data_tran)
data_plot_tran <- data_tran %>% 
  rownames_to_column('Genes') %>%  
  reshape2::melt(
    id.vars = 'Genes',
    measure.vars = colnames(.)[2:ncol(.)],
    variable.name = 'sample',
    value.name = 'copy_num'
  ) %>% 
  dplyr::filter(copy_num>0) %>%
  left_join(meta_data_tran,by = 'sample') %>% 
  mutate(
    log10_cp = log10(copy_num),
    group_omics = 'transcription'
  ) %>% 
  dplyr::select(
    Genes,sample,copy_num,group_cell,
    group_stage,group_sample,log10_cp, group_omics
  )

colnames(meta_data_pro)
colnames(data_pro)
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
  dplyr::select(
    Genes,sample,copy_num,group_cell,
    group_stage,group_sample,log10_cp, group_omics
  )
data_plot <- rbind(data_plot_tran,data_plot_pro) %>% 
  mutate(
    group_omics_cell = paste(group_omics,group_cell,sep = '_')
  ) %>% 
  mutate(group_omics_cell = factor(group_omics_cell,levels = unique(group_omics_cell)))


pal <- c('#A3D2E2','#86C0CB','#6da0ac','#4d7179')
pal <- c('#CF9698', '#B76265', '#BF7578', '#865254')
pal <- c('#A4BDDA', '#7ea1cc', '#6682a4','#4d627c')
color_use <- c(
  '#865254','#4d627c','#CF9698','#A4BDDA'
)
tmp <- ggplot(data_plot, aes(x = log10_cp,fill = group_omics_cell)) + 
  geom_density(alpha = 0.9,color = 'grey96',linewidth = 1.5)
density_data <- ggplot_build(tmp)$data[[1]]

data_plot %>% colnames()
data_median <- data_plot %>% 
  group_by(group_omics_cell) %>% 
  summarise(median_value = median(log10_cp)) %>% 
  rowwise() %>% 
  mutate(
    y = density_data %>% 
      filter(group == as.integer(as.factor(group_omics_cell))) %>% 
      filter(abs(x - median_value) == min(abs(x - median_value))) %>% 
      pull(y)
  )
data_median
data_median_label <- data_median %>% 
  as.data.frame() %>% 
  mutate(
    label_x = c(0.71,0.5,5.1,6.22),
    label = paste('Median:',median_value %>% round(digits = 2))
  )
data_median_label
data_range <- data_plot %>%
  group_by(group_omics_cell) %>%
  summarise(
    lower_90 = quantile(log10_cp, 0.025),
    upper_90 = quantile(log10_cp, 0.975)
  ) %>%

  mutate(
    log10_span = paste(
      (upper_90 - lower_90) %>% round(digits = 1),
      'OM',sep = ' '
    )
  )
data_range$y_loc <- c(0.05,0.1,0.05,0.1)
data_range
ggplot(data_plot, aes(x = log10_cp,fill = group_omics_cell)) + 
  geom_density(alpha = 0.9,color = 'grey96',linewidth = 1.5) +  
  geom_segment(
    data = data_range,
    mapping = aes(
      x = lower_90,xend = upper_90,y = y_loc,yend = y_loc
      ),
    color = 'grey80',
    arrow = arrow(
      angle = 20,length = unit(3,units = 'pt'),
      type = "closed", ends = "both"
    )
  ) +
  geom_text(
    data = data_range,
    mapping = aes(
      x = (upper_90+lower_90)/2,y = y_loc + 0.015,label = log10_span
    ),
    color = 'black',
  ) +
  geom_segment(
    data = data_median, 
    aes(x = median_value, xend = median_value,
        y = 0, yend = y),
    linetype = "dashed", size = 1,
    color = 'white',show.legend = FALSE) +
  geom_point(
    data = data_median,
    aes(x = median_value, y = y),
    size = 3,color = 'black',show.legend = FALSE) +
  scale_fill_manual(values = color_use) +
  geom_text_repel(
    data = data_median_label,
    aes(x = label_x, y = 0.2,label = label),
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0.5,
    force = 90,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = Inf
  ) +
  labs(title = "Protein & Transcription", x = "Log10(abundance)", y = "Density") +
  scale_y_continuous(
    expand = c(0,0),
    breaks = seq(0,0.6,0.1),
    limits = c(0,0.6)
  ) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0.5,color = 'white'))) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.2,0.9),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20,margin = margin(b = 6,unit = 'pt')),
    strip.clip = 'on',
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1F_proteins_transcript.png',
  width = 12,height = 6
)

ggsave(
  filename = './fig1F_proteins_transcript.pdf',
  width = 12,height = 6
)




library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(paletteer)
setwd('./protein/')
dir()
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
df <- read.csv('./follicleresult/report.tsv',sep = '\t')
df<-df[df$Proteotypic == 1 & df$Protein.Q.Value <= 0.01,]
colnames(df)
data_protein<-df[,c('File.Name',"Run" ,"Precursor.Id","PG.Quantity")]
data_protein<-data_protein[!duplicated(data_protein),]
sample_level <- c(
  "SGC1","SGC2","SGC3","SGC4","SGC5","SGC6","SGC7","SGC8","SGC9",
  "SGC10","SGC11","SGC12","SO1","SO2","SO3","SO4","SO5","SO6","SO7",
  "SO8","SO9","SO10","SO11","SO12","EGC1","EGC2","EGC3","EGC4",
  "EGC5","EGC6","EGC7","EGC9","EGC10","EGC11","EGC12","EGC13","EO1",
  "EO2","EO3","EO4","EO5","EO6","EO7","EO9","EO10","EO11","EO12",
  "EO13","AGC2","AGC3","AGC4","AGC5","AGC6","AGC7","AGC8","AGC9",
  "AGC10","AGC11","AGC12","AO2","AO3","AO4","AO5","AO6","AO7",
  "AO8","AO9","AO10","AO11","AO12","PGC1","PGC2","PGC3","PGC4","PGC5",
  "PGC6","PGC7","PGC8","PGC9","PGC10","PO1","PO2","PO3","PO4","PO5",
  "PO6","PO7","PO8","PO9","PO10"#"AO1",
)

data_protein$File.Name[1:3]
colnames(data_protein)
data_protein_res<-df[,c('File.Name',"Run" ,"Precursor.Id","PG.Quantity")] %>% 
  .[!duplicated(.),] %>% 
  mutate(
    section = File.Name %>% 
      str_split(pattern = '_',n=16) %>% 
      do.call(rbind,.) %>% 
      .[,11]
  ) %>% 
  select(section,Precursor.Id,PG.Quantity) %>% 
  mutate(section = str_replace(section,'P0','PO')) %>% 
  filter(section != 'AO1') %>% 
  pivot_wider(names_from = section,values_from = PG.Quantity) %>% 
  rename("PG.ProteinGroups" = Precursor.Id)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#7698b3","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#d6616b","#a55194","#ce6dbd","#e7ba52",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252"
) 
data_clean <- data_protein_res

sample_level <- c(
  "SGC1","SGC2","SGC3","SGC4","SGC5","SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "SO1","SO2","SO3","SO4","SO5","SO6","SO7","SO8","SO9","SO10","SO11","SO12",
  "EGC2","EGC3","EGC4","EGC5","EGC6","EGC7",#"EGC1","EGC9","EGC10","EGC11","EGC12","EGC13",
  "EO2","EO3","EO4","EO5","EO6","EO7",#"EO1","EO9","EO10","EO11","EO12","EO13",
  "AGC2","AGC3","AGC4","AGC5","AGC7","AGC8","AGC9","AGC10","AGC11","AGC12",#"AGC6",
  "AO2","AO3","AO4","AO5","AO7","AO8","AO9","AO10","AO11","AO12",#"AO6",
  "PGC1","PGC2","PGC3","PGC4","PGC5","PGC6","PGC7","PGC8","PGC9",#"PGC10",
  "PO1","PO2","PO3","PO4","PO5","PO6","PO7","PO8","PO9"
)
data <- data_clean[,c(colnames(data_clean)[c(1, 2)],sample_level)]

sample_level <- c(
  "SGC1","SGC2","SGC3","SGC4","SGC5","SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "SO1","SO2","SO3","SO4","SO5","SO6","SO7","SO8","SO9","SO10","SO11","SO12",
  "EGC2","EGC3","EGC4","EGC5","EGC6","EGC7",#"EGC1","EGC9","EGC10","EGC11","EGC12","EGC13",
  "EO2","EO3","EO4","EO5","EO6","EO7",#"EO1","EO9","EO10","EO11","EO12","EO13",
  "AGC2","AGC3","AGC4","AGC5","AGC7","AGC8","AGC9","AGC10","AGC11","AGC12",#"AGC6",
  "AO2","AO3","AO4","AO5","AO7","AO8","AO9","AO10","AO11","AO12",#"AO6",
  "PGC1","PGC2","PGC3","PGC4","PGC5","PGC6","PGC7","PGC8","PGC9",#"PGC10",
  "PO1","PO2","PO3","PO4","PO5","PO6","PO7","PO8","PO9"
)
data <- data_clean[,c(colnames(data_clean)[c(1, 2)],sample_level)]
group_stage_name <- c('Secondary','Early antral','Antral','Preovulatory')
names(group_stage_name) <- c('S','E','A','P')
group_cell_name <- c('Oocyte','Granulosa')
names(group_cell_name) <- c('O','G')
sample_meta_data <- 
  lapply(colnames(data)[-c(1, 2)], function(sample) {
    sum(!is.na(data[, sample]))
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1, 2)]) %>% 
  dplyr::select(sample) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))

protein_num <-
  lapply(colnames(data)[-c(1, 2)], function(sample) {
    sum(!is.na(data[, sample]))
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1, 2)]) %>% 
  dplyr::select(sample,num_pro) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))
stage_cell <- c('all','all_Granulosa','all_Oocyte',unique(protein_num$group_stage_cell))
data_plot <- lapply(stage_cell, function(stage_cell){
  if(stage_cell == 'all'){
    pro_num <- protein_num$num_pro %>% sum() %>% {./nrow(protein_num)}
  }else if(stage_cell == 'all_Granulosa'){
    pro_num <- protein_num %>% 
      filter(group_cell == 'Granulosa') %>% 
      {sum(.$num_pro)/nrow(.)}
  }else if(stage_cell == 'all_Oocyte'){
    pro_num <- protein_num %>% 
      filter(group_cell == 'Oocyte') %>% 
      {sum(.$num_pro)/nrow(.)}
  }else{
    pro_num <- protein_num %>% 
      filter(group_stage_cell == stage_cell) %>% 
      {sum(.$num_pro)/nrow(.)}
  }
  return(pro_num)
}) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(
    group_stage_cell = stage_cell,
    sample = stage_cell
  ) %>% 
  dplyr::select(group_stage_cell,sample,num_pro) %>% 
  rbind(protein_num[,c('group_stage_cell','sample','num_pro')])
data_plot$group_stage_cell[data_plot$sample %in% c('all_Granulosa', 'all_Oocyte')] <-
  'all'
data_plot$sample %>% unique()
sample_level <- c(
  "rep1","rep2","rep3","rep4","all","all_Granulosa","all_Oocyte",
  "rep7","rep8","rep9", "rep10",
  "Secondary: Granulosa","SGC1","SGC2","SGC3","SGC4","SGC5",
  "SGC6","SGC7","SGC8","SGC9","SGC10","SGC11","SGC12",
  "Secondary: Oocyte","SO1","SO2","SO3","SO4","SO5","SO6",
  "SO7","SO8","SO9","SO10","SO11","SO12",
  "Early antral: Granulosa","EGC2","EGC3","EGC4","EGC5",
  "EGC6","EGC7",
  "Early antral: Oocyte","EO2","EO3","EO4","EO5","EO6","EO7",
  "Antral: Granulosa","AGC2","AGC3","AGC4","AGC5","AGC7",
  "AGC8","AGC9","AGC10","AGC11","AGC12",
  "Antral: Oocyte","AO2","AO3","AO4","AO5","AO7","AO8",
  "AO9","AO10","AO11","AO12",
  "Preovulatory: Granulosa","PGC1","PGC2","PGC3","PGC4",
  "PGC5","PGC6","PGC7","PGC8","PGC9",
  "Preovulatory: Oocyte","PO1","PO2","PO3","PO4","PO5",
  "PO6","PO7","PO8","PO9"
)
color_grey_sample <- c(
  "all","all_Granulosa","all_Oocyte",
  "Secondary: Granulosa","Secondary: Oocyte",
  "Early antral: Granulosa","Early antral: Oocyte", 
  "Antral: Granulosa", "Antral: Oocyte",
  "Preovulatory: Granulosa","Preovulatory: Oocyte")
color_use <- lapply(sample_level, function(sample){
  if(sample %in% color_grey_sample){
    return('lightslategrey')
  }else{
    sample_group <- sample_meta_data$group_stage_cell[sample_meta_data$sample == sample]
    color_tmp <- color_group_stage_cell[sample_group]
  }
  return(color_tmp)
}) %>% unlist() %>% as.character()
color_use[c(2,3)] <- c(color_gc,color_oo)
data_plot$sample <- factor(data_plot$sample,levels = sample_level)
data_plot$group_stage_cell <- factor(data_plot$group_stage_cell,levels = c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
))

label <- c('Mean','SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "All","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell],
    label_text = factor(label_text,levels = label %>% as.character())
  )

data_plot$label_text %>% table()
ggplot(
  data = data_plot,
  mapping = aes(x = label_text, y = num_pro, fill = sample)
) +
  geom_bar(
    stat = 'identity',
    position = position_dodge2(preserve = "single",padding = 0),
    color = 'white') +
  scale_y_continuous(
    labels = c(expression(italic(0)),
               expression(1%*%10^4),
               expression(2%*%10^4),
               expression(3%*%10^4),
               expression(4%*%10^4),
               expression(5%*%10^4)),
    expand = c(0,0),
    breaks = c(0,1%*%10^4,2%*%10^4,3%*%10^4,4%*%10^4,5%*%10^4),
    limits = c(0,5.1*10^4)
  ) +
  scale_fill_manual(values = color_use) +
  labs(y = 'No. of identified peptides') +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = 'none',
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      size = 18,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(
      size = 18, margin = margin(b = 6, unit = 'pt')
    ),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1b_peptide_counts.png',
  width = 18,height = 7
)
ggsave(
  filename = './fig1b_peptide_counts.pdf',
  width = 18,height = 7
)



library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Seurat)
library(RColorBrewer)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bd2","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
color_gc <- '#477BD1';color_oo <- '#E25247'
color_group_stage_cell <- c(
  '#2CA02C','#0088D4','#17BECF','#006BA4',
  '#F8B620','#FF800E','#D14B70','#ED444A'
)
names(color_group_stage_cell) <- c(
  "Secondary: Granulosa", "Early antral: Granulosa",
  "Antral: Granulosa", "Preovulatory: Granulosa",
  "Secondary: Oocyte","Early antral: Oocyte",
  "Antral: Oocyte", "Preovulatory: Oocyte"
)

meta_data <- read.csv('./data_meta.csv',row.names = 1)
data_process <- read.csv('./data_fragment.csv',row.names = 1)

sample_level <- rownames(meta_data)
data <- data_process[,sample_level] %>% 
  rownames_to_column('genes')
group_stage_name <- c('Secondary','Early antral','Antral','Preovulatory')
names(group_stage_name) <- c('S','E','A','P')
group_cell_name <- c('Oocyte','Granulosa')
names(group_cell_name) <- c('O','G')
sample_meta_data <- 
  lapply(colnames(data)[-c(1)], function(sample) {
    sum(data[, sample]>0)
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1)]) %>% 
  dplyr::select(sample) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))

protein_num <-
  lapply(colnames(data)[-c(1)], function(sample) {
    sum(data[, sample]>0)
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1)]) %>% 
  dplyr::select(sample,num_pro) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))


stage_cell <- c('all','all_Granulosa','all_Oocyte',unique(protein_num$group_stage_cell))
data_plot <- lapply(stage_cell, function(stage_cell){
  if(stage_cell == 'all'){
    sample <- protein_num$sample
  }else if(stage_cell == 'all_Granulosa'){
    sample <- protein_num$sample[grepl('Granulosa',protein_num$group_stage_cell)] %>% 
      as.character()
  }else if(stage_cell == 'all_Oocyte'){
    sample <- protein_num$sample[grepl('Oocyte',protein_num$group_stage_cell)] %>% 
      as.character()
  }else{
    sample <- protein_num$sample[protein_num$group_stage_cell==stage_cell] %>% 
      as.character()
  }
  data_tmp <- data[,colnames(data) %in% sample]
  data_tmp[is.na(data_tmp)] <- 0
  pro_num <- sum(rowSums(data_tmp)>0)
  return(pro_num)
}) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(
    group_stage_cell = stage_cell,
    sample = stage_cell
  ) %>% 
  dplyr::select(group_stage_cell,sample,num_pro) %>% 
  rbind(protein_num[,c('group_stage_cell','sample','num_pro')])
data_plot$group_stage_cell[data_plot$sample %in% c(
  'all_Granulosa','all_Oocyte'
)] <- 'all'
data_plot$sample %>% unique()
sample_level <- c(
  "rep1","rep2","rep3","rep4","all","all_Granulosa","all_Oocyte",
  "rep7","rep8","rep9", "rep10",
  "Secondary: Granulosa","SGC3","SGC4","SGC5","SGC7","SGC8","SGC11","SGC21","SGC24",
  "Secondary: Oocyte","SO3","SO4","SO5","SO7","SO8","SO11","SO21","SO24",
  "Early antral: Granulosa","EGC6","EGC7","EGC21","EGC23","EGC24","EGC25","EGC26",
  "Early antral: Oocyte","EO6","EO7","EO21","EO23","EO24","EO25","EO26",
  "Antral: Granulosa","AGC2","AGC4","AGC5","AGC9","AGC10","AGC24",
  "Antral: Oocyte","AO2","AO4","AO5","AO9","AO10","AO24",
  "Preovulatory: Granulosa","PGC1","PGC2","PGC3","PGC4","PGC5","PGC6","PGC9",
  "PGC10","PGC11","PGC22","PGC23",
  "Preovulatory: Oocyte","PO1","PO2","PO3","PO4","PO5","PO6","PO9","PO10",
  "PO11","PO22","PO23"
)
sample_level <- sample_level[sample_level %in% (data_plot$sample %>% unique())]
color_grey_sample <- c(
  "all","all_Granulosa","all_Oocyte",
  "Secondary: Granulosa","Secondary: Oocyte",
  "Early antral: Granulosa","Early antral: Oocyte", 
  "Antral: Granulosa", "Antral: Oocyte",
  "Preovulatory: Granulosa","Preovulatory: Oocyte")
color_use <- lapply(sample_level, function(sample){
  if(sample %in% color_grey_sample){
    return('lightslategrey')
  }else{
    sample_group <- sample_meta_data$group_stage_cell[sample_meta_data$sample == sample]
    color_tmp <- color_group_stage_cell[sample_group]
  }
  return(color_tmp)
}) %>% unlist() %>% as.character()
color_use[c(2,3)] <- c(color_gc,color_oo)
data_plot$sample <- factor(data_plot$sample,levels = sample_level)
data_plot$group_stage_cell <- factor(data_plot$group_stage_cell,levels = c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
))

label <- c('All','SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell],
    label_text = factor(label_text,levels = label %>% as.character())
  )
write.csv(data_plot,'./data_plot_trans_gene-number.csv')
ggplot(
  data = data_plot,
  mapping = aes(x = label_text, y = num_pro, fill = sample)
) +
  geom_bar(
    stat = 'identity',
    position = position_dodge2(preserve = "single",padding = 0),
    color = 'white') +
  geom_hline(yintercept = 11000,linewidth = 1.5,color = 'grey80',linetype = 'longdash') +
  geom_hline(yintercept = 8000,linewidth = 1.5,color = 'grey80',linetype = 'longdash') +
  scale_y_continuous(
    labels = c(expression(italic(0)),
               expression(1%*%10^4),
               expression(2%*%10^4),
               expression(3%*%10^4),
               expression(4%*%10^4),
               expression(4.5%*%10^4)),
    expand = c(0,0),
    breaks = c(0,1*10^4,2*10^4,3*10^4,4*10^4,4.5*10^4),
  ) +
  scale_fill_manual(values = color_use) +
  labs(y = 'No. of identified genes') +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = 'none',
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      size = 18,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(
      size = 18, margin = margin(b = 6, unit = 'pt')
      ),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1b_tran.png',
  width = 12,height = 4
)
ggsave(
  filename = './fig1b_tran_0206.pdf',
  width = 12,height = 4
)

sample_level <- rownames(meta_data)
data <- data_process[,sample_level] %>% 
  rownames_to_column('genes')

group_stage_name <- c('Secondary','Early antral','Antral','Preovulatory')
names(group_stage_name) <- c('S','E','A','P')
group_cell_name <- c('Oocyte','Granulosa')
names(group_cell_name) <- c('O','G')
sample_meta_data <- 
  lapply(colnames(data)[-c(1)], function(sample) {
    sum(data[, sample]>0)
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1)]) %>% 
  dplyr::select(sample) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))


protein_num <-
  lapply(colnames(data)[-c(1)], function(sample) {
    sum(data[, sample])
  }) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(sample = colnames(data)[-c(1)]) %>% 
  dplyr::select(sample,num_pro) %>% 
  mutate(
    group_cell = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,2,2)}) %>% unlist()
  ) %>% 
  mutate(
    group_cell = group_cell_name[group_cell]
  ) %>% 
  mutate(
    group_stage = sample %>% lapply(.,function(x){x[1] %>% as.character() %>% substr(.,1,1)}) %>% unlist()
  ) %>% 
  mutate(
    group_stage = group_stage_name[group_stage] %>% factor(.,levels = group_stage_name)
  ) %>% 
  mutate(
    group_sample = sample %>% str_remove('O') %>% str_remove('GC')
  ) %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_sample = factor(group_sample,levels = group_sample %>% unique())) %>% 
  mutate(sample = factor(sample,levels = sample_level))


stage_cell <- c('all','all_Granulosa','all_Oocyte',unique(protein_num$group_stage_cell))
data_plot <- lapply(stage_cell, function(stage_cell){
  if(stage_cell == 'all'){
    sample <- protein_num$sample
  }else if(stage_cell == 'all_Granulosa'){
    sample <- protein_num$sample[grepl('Granulosa',protein_num$group_stage_cell)] %>% 
      as.character()
  }else if(stage_cell == 'all_Oocyte'){
    sample <- protein_num$sample[grepl('Oocyte',protein_num$group_stage_cell)] %>% 
      as.character()
  }else{
    sample <- protein_num$sample[protein_num$group_stage_cell==stage_cell] %>% 
      as.character()
  }
  data_tmp <- data[,colnames(data) %in% sample]
  data_tmp[is.na(data_tmp)] <- 0
  pro_num <- sum(rowSums(data_tmp)/length(sample))
  return(pro_num)
}) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  rename_all(~c('num_pro')) %>% 
  mutate(
    group_stage_cell = stage_cell,
    sample = stage_cell
  ) %>% 
  dplyr::select(group_stage_cell,sample,num_pro) %>% 
  rbind(protein_num[,c('group_stage_cell','sample','num_pro')])
data_plot$group_stage_cell[data_plot$sample %in% c(
  'all_Granulosa','all_Oocyte'
)] <- 'all'
data_plot$sample %>% unique()
sample_level <- c(
  "rep1","rep2","rep3","rep4","all","all_Granulosa","all_Oocyte",
  "rep7","rep8","rep9", "rep10",
  "Secondary: Granulosa","SGC3","SGC4","SGC5","SGC7","SGC8","SGC11","SGC21","SGC24",
  "Secondary: Oocyte","SO3","SO4","SO5","SO7","SO8","SO11","SO21","SO24",
  "Early antral: Granulosa","EGC6","EGC7","EGC21","EGC23","EGC24","EGC25","EGC26",
  "Early antral: Oocyte","EO6","EO7","EO21","EO23","EO24","EO25","EO26",
  "Antral: Granulosa","AGC2","AGC4","AGC5","AGC9","AGC10","AGC24",
  "Antral: Oocyte","AO2","AO4","AO5","AO9","AO10","AO24",
  "Preovulatory: Granulosa","PGC1","PGC2","PGC3","PGC4","PGC5","PGC6","PGC9",
  "PGC10","PGC11","PGC22","PGC23",
  "Preovulatory: Oocyte","PO1","PO2","PO3","PO4","PO5","PO6","PO9","PO10",
  "PO11","PO22","PO23"
)
sample_level <- sample_level[sample_level %in% (data_plot$sample %>% unique())]
color_grey_sample <- c(
  "all","all_Granulosa","all_Oocyte",
  "Secondary: Granulosa","Secondary: Oocyte",
  "Early antral: Granulosa","Early antral: Oocyte", 
  "Antral: Granulosa", "Antral: Oocyte",
  "Preovulatory: Granulosa","Preovulatory: Oocyte")
color_use <- lapply(sample_level, function(sample){
  if(sample %in% color_grey_sample){
    return('lightslategrey')
  }else{
    sample_group <- sample_meta_data$group_stage_cell[sample_meta_data$sample == sample]
    color_tmp <- color_group_stage_cell[sample_group]
  }
  return(color_tmp)
}) %>% unlist() %>% as.character()
color_use[c(2,3)] <- c(color_gc,color_oo)
data_plot$sample <- factor(data_plot$sample,levels = sample_level)
data_plot$group_stage_cell <- factor(data_plot$group_stage_cell,levels = c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
))

label <- c('Mean','SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "all","Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell],
    label_text = factor(label_text,levels = label %>% as.character())
  )
write.csv(data_plot,'./data_plot_trans_gene-fragments.csv')

ggplot(
  data = data_plot,
  mapping = aes(x = label_text, y = num_pro, fill = sample)
) +
  geom_bar(
    stat = 'identity',
    position = position_dodge2(preserve = "single",padding = 0),
    color = 'white') +
  scale_y_continuous(
    labels = c(expression(italic(0)),
               expression(1%*%10^7),
               expression(1.5%*%10^7),
               expression(2%*%10^7),
               expression(2.5%*%10^7),
               expression(3%*%10^7)),
    expand = c(0,0),
    breaks = c(0,1*10^7,1.5*10^7,2*10^7,2.5*10^7,3*10^7),
    limits = c(0,3.2*10^7)
  ) +
  scale_fill_manual(values = color_use) +
  labs(y = 'Expression(fragments)') +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    legend.position = 'none',
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      size = 18,
      hjust = 0.5
    ),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(
      size = 18, margin = margin(b = 6, unit = 'pt')
    ),
    strip.clip = 'on',
    strip.background = element_rect(
      fill = NA,
      linewidth = 0,
      colour = NA
    ),
    strip.placement = 'outside'
  )
ggsave(
  filename = './fig1b_tran_fragments.pdf',
  width = 12,height = 4
)
ggsave(
  filename = './fig1b_tran_fragments_0206.pdf',
  width = 12,height = 4
)

colnames(meta_data)
meta_data$group_stage_cell <- paste(meta_data$group_stage,meta_data$group_cell,sep = ': ')
group_stage_cell <- meta_data$group_stage_cell %>% unique()
data_plot <- lapply(group_stage_cell, function(group_use) {
  sample <- meta_data %>% 
    filter(group_stage_cell %in% group_use) %>% 
    pull(sample)
  data_tmp <- data_process[, colnames(data_process) %in% sample]
  protein_all <- apply(data_tmp, 2, function(sample) {
    rownames(data_tmp)[sample > 0] %>% .[!is.na(.)]
  }) %>% unlist()
  protein_all_not_dup <- protein_all %>% unique()
  protein_all_dup <-
    protein_all %>% table() %>% as.data.frame() %>%
    rename_all( ~ c('protein', 'num')) %>% as.data.frame()
  data_res <- table(protein_all_dup$num) %>% 
    as.data.frame() %>%
    rename_all( ~ c('sample_num', 'gene_num')) %>%
    mutate(group_stage_cell = group_use,
           duplicated_rate = gene_num / sum(gene_num) * 100)
  return(data_res)
}) %>% do.call(rbind, .) %>%
  mutate(
    group_stage_cell = factor(
      group_stage_cell,levels = c(
        "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
        "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
      )),
    sample_num = factor(
      sample_num %>% as.numeric(),
      levels = 1:(sample_num %>% as.numeric() %>% max())# %>% rev()
    )
  ) %>% 
  group_by(group_stage_cell) %>% 
  mutate(
    tmp = sample_num %>% as.character() %>% as.numeric(),
    group_color = case_when(
      (tmp > max(tmp)*0.75) ~ '75%~100%',
      (tmp > max(tmp)*0.5) ~ '50%~75%',
      (tmp > max(tmp)*0.25) ~ '25%~50%',
      TRUE ~ '0~25%'
    )
  ) %>% 
  group_by(group_stage_cell,group_color) %>% 
  mutate(
    duplicated_rate_sum = sum(duplicated_rate)
  ) %>% 
  select(-sample_num,-tmp,-duplicated_rate,-gene_num) %>% distinct() %>% 
  group_by(group_stage_cell) %>%
  arrange(group_color) %>%
  mutate(label_loc = cumsum(duplicated_rate_sum) - duplicated_rate_sum / 2) %>% 
  dplyr::rename(duplicated_rate = duplicated_rate_sum)

label <- c('SO','EO','AO','PO','SG','EG','AG','PG')
names(label) <- c(
  "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
  "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
)
data_plot <- data_plot %>% 
  ungroup() %>% 
  mutate(
    label_text = label[data_plot$group_stage_cell] %>% as.character(),
    label_text = factor(label_text,levels = label %>% as.character())
  )
library(paletteer)
ggplot(
  data = data_plot,
  aes(x = label_text,y = -duplicated_rate)
) +
  geom_bar(
    aes(group = rev(group_color),fill = group_color),
    stat = 'identity',color = 'black') +
  geom_text(
    aes(
      x = label_text,
      y = -label_loc,
      label = paste(round(duplicated_rate, 1),'%',sep = '')
    ),
    color = 'black',
    size = 4
  ) +
  scale_fill_manual(
    values = paletteer_c("grDevices::Oslo", 50) %>% .[c(5,10,15,30)]
  ) +
  scale_y_continuous(
    expand = c(0,0.1),
    breaks = seq(-100,0,25),
    label = seq(0,100,25)
    ) +
  scale_x_discrete(expand = c(0,0.5)) +
  labs(y = 'repeatability (%)') +
  theme_classic() +
  theme(
    plot.margin = margin(l = 25,t = 10),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0,size = 18,hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 24),
    strip.background = element_rect(fill = NA,linewidth = 0,colour = NA),
    strip.placement = 'outside'
  )



dir.create('./figure1_20250508',showWarnings = FALSE)
ggsave(
  filename = './fig1c_transcript.png',
  width = 8,height = 8
)
ggsave(
  filename = './fig1c_transcript.pdf',
  width = 8,height = 8
)

write.csv(data_plot,'./data_fig1c_transcript.csv',row.names = FALSE)


library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
library(purrr)
library(stringr)
marker_tabel <- readxl::read_xlsx('./marker gene.xlsx')
scdata <- readRDS('./scdata_tpm.rds')
meta_data <- scdata@meta.data
data_copy_num_sorted <- scdata@assays$RNA$data_tpm %>%
  rownames_to_column('Genes')
colnames(data_copy_num_sorted)
colnames(meta_data)
data_plot <- data_copy_num_sorted %>%
  dplyr::select(
    c(
      'Genes',
      meta_data$sample %>% as.character()
    )) %>% 
  pivot_longer(
    cols = names(.)[-c(1)],
    names_to = 'sample',
    values_to = 'copy_number') %>% 
  left_join(
    meta_data %>% dplyr::select(c(group_stage,group_cell,sample)),
    by = 'sample'
  ) %>% 
  group_by(group_cell,Genes) %>% 
  summarise(
    Genes = Genes[1],
    copy_number = mean(copy_number),
    .groups = 'keep'
  ) %>% 
  filter(copy_number>0) %>% 
  group_by(group_cell) %>% #group_stage
  arrange(desc(copy_number),.by_group = TRUE) %>%
  mutate(
    rank_num = row_number(),
    copy_number = log10(copy_number+1)
  ) %>% 
  arrange(rank_num) %>% 
  mutate(
    cumulative_sum = cumsum(copy_number)/sum(copy_number)
  )
colnames(data_plot)

data_plot_sub <- data_plot

data_summary <- data_plot %>% 
  mutate(
    group_rank = case_when(
      cumulative_sum <= 0.25 ~ 'Q1',
      cumulative_sum <= 0.5 ~ 'Q2',
      cumulative_sum <= 0.75 ~ 'Q3',
      cumulative_sum <= 1 ~ 'Q4',
    )
  ) %>% 
  group_by(group_cell,group_rank) %>% 
  summarise(
    Counts = n(),.groups = 'drop'
  ) %>% 
  arrange(group_rank, factor(group_cell, levels = c("Oocyte", "Granulosa"))) %>%
  group_by(group_rank) %>%
  reframe(
    label = paste0(group_rank, "[", paste(group_cell, Counts, sep = ": ", collapse = ";"), "]")
  ) %>%
  distinct() %>% 
  mutate(
    x_loc = 7500,
    y_loc = seq(0.05,1,by = 0.25)
  )

data_plot_sub_label <- data_plot_sub %>% 
  group_by(group_cell) %>%
  arrange(rank_num) %>%
  dplyr::slice(1:10) %>% 
  mutate(rank_num = 1:10) 
colnames(data_plot_sub)
ggplot() + 
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0,
      ymax = 0.25
    ),
    fill = "#438E80",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.25,
      ymax = 0.5
    ),
    fill = "#336C61",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.5,
      ymax = 0.75
    ),
    fill = "#E69D70",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +  
  geom_rect(
    aes(
      xmin = data_plot_sub$rank_num %>% min(),
      xmax = data_plot_sub$rank_num %>% max(),
      ymin = 0.75,
      ymax = 1
    ),
    fill = "#B65A20",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data_summary,
    mapping = aes(x = 39500,y = y_loc,label = label),#
    size = 5,hjust = 1
  ) +
  geom_point(
    data = data_plot_sub,
    mapping = aes(
      x = rank_num,
      y = cumulative_sum,
      color = group_cell)
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = cumulative_sum),
    color = 'red',size = 3
  ) +
  geom_text(
    data = data_plot_sub_label[data_plot_sub_label$group_cell == 'Oocyte',],
    mapping = aes(x = 5400,y = rank_num*0.05,label = Genes),
    size = 5,hjust = 0,color = '#865254'
  ) +
  geom_text(
    data = data_plot_sub_label[data_plot_sub_label$group_cell == 'Granulosa',],
    mapping = aes(x = 11000,y = rank_num*0.05,label = Genes),
    size = 5,hjust = 0,color = '#4d627c'
  ) +
  guides(color = guide_legend(
    title = 'Cell Type',
    override.aes = list(size = 4)
  )) +
  scale_color_manual(values = c('#4d627c','#865254')) +
  labs(
    x = 'Abundance Rank', 
    y = 'Cumulative Abundance[%]', 
    title = 'Transcription') +
  scale_x_continuous(
    expand = c(0,100),
  ) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )

ggsave(
  filename = './fig1f_trans.png',
  width = 8,height = 12
)
ggsave(
  filename = './fig1f_trans.pdf',
  width = 8,height = 12
)
data_plot_sub <- data_plot %>% 
  filter(group_cell == 'Oocyte')
marker_select <- marker_tabel %>% 
  filter(`cell type` == 'oocyte') %>% 
  pull(`marker gene`) %>% 
  unique() %>% str_to_title()
data_plot_sub_label <- data_plot_sub %>% 
  filter(
    Genes %in% marker_select
  )
ggplot(data = data_plot_sub) + 
  geom_point(
    aes(x = rank_num, y = copy_number), 
    color = '#C34F73'
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = copy_number), 
    color = 'black',size = 3
  ) +
  geom_text_repel(
    data = data_plot_sub_label,
    mapping = aes(
      x = rank_num, 
      y = copy_number,
      label = Genes
    ),force_pull = 2,nudge_x = 10,size = 6
  )+
  labs(x = 'Gene Rank', y = 'Log10(TPM + 1)', title = 'Oocyte') +
  scale_x_continuous(expand = c(0,1000)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    plot.margin = margin(r = 16),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
ggsave(
  filename = './fig1d_transcript_Oocyte.png',
  width = 6,height = 10
)
data_plot_sub <- data_plot %>% 
  filter(group_cell == 'Granulosa')
marker_select <- marker_tabel %>% 
  filter(`cell type` == 'GC') %>% 
  pull(`marker gene`) %>% 
  unique() %>% str_to_title()
data_plot_sub_label <- data_plot_sub %>% 
  filter(
    Genes %in% marker_select
  )
ggplot(data = data_plot_sub) + 
  geom_point(
    aes(x = rank_num, y = copy_number), 
    color = '#0088D4'
  ) +
  geom_point(
    data = data_plot_sub_label,
    mapping = aes(x = rank_num, y = copy_number), 
    color = 'black',size = 3
  ) +
  geom_text_repel(
    data = data_plot_sub_label,
    mapping = aes(
      x = rank_num, 
      y = copy_number,
      label = Genes
    ),force_pull = 2,nudge_x = 10,size = 6
  )+
  labs(x = 'Gene Rank', y = 'Log10(TPM + 1)', title = 'Granulosa') +
  scale_x_continuous(expand = c(0,1000)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(
      length = unit(0.25, units = 'mm'), angle = 30
    )),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    plot.title = element_text(size = 24,hjust = 0.5,vjust = 0.5)
  )
ggsave(
  filename = './fig1d_transcript_Granulosa.png',
  width = 6,height = 10
)

scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
data_tran <- scdata_tran@assays$RNA$counts %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]

library(reshape)
library(ggplot2)
library(paletteer)
library(ggpubr)
data_nor_tran_melt<-reshape2::melt(
  data_tran %>% rownames_to_column('Genes'),
  id.vars ="Genes",
  variable.name = 'sample',
  value.name = 'tpm') %>% 
  mutate(
    tpm = as.numeric(tpm)
  ) %>% 
  mutate(
    tpm_log10 = log10(tpm+1)
  ) %>% 
  left_join(
    meta_data_tran %>% dplyr::select(sample, group_stage_cell),
    by = 'sample') %>% 
  mutate(
    group_stage_cell = factor(
      group_stage_cell,
      levels = meta_data_tran$group_stage_cell %>% unique()
    )
  )
ggboxplot(
  data = data_nor_tran_melt[data_nor_tran_melt$tpm_log10>0,],
  x = "sample",
  y = "tpm_log10",outlier.shape = NA,
  fill = "group_stage_cell",
  xlab = "Sample",
  ylab = "log10(TPM+1)",
  title = "Transcription"
)+ 
  scale_y_continuous(limits = c(-0.5,2)) +
  scale_fill_manual(values = my36colors)+
  theme(
    axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5),
    legend.position="top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )
ggsave('./transcript_normalization_tpm.png',width = 18,height = 7)
ggsave('./transcript_normalization_tpm.pdf',width = 18,height = 7)


data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
data_tf_list <- c(data_raw_1$Symbol,data_raw_2$tf) %>% unique


data_sm_list <- read.csv('./maternal_gene.txt',header = FALSE)
human <- readRDS('D:/jiaoxi/ssgsea/gene_tran/human.rds')
mouse <- readRDS('D:/jiaoxi/ssgsea/gene_tran/mouse.rds')
data_sm_list <-
  biomaRt::getLDS(
    attributes = "hgnc_symbol",
    filters = "hgnc_symbol",
    values = data_sm_list,
    mart = human,
    attributesL = "mgi_symbol",
    martL = mouse,
    uniqueRows = TRUE
  )
saveRDS(data_sm_list,'./data_sm_list.rds')
data_sm_list <- data_sm_list$MGI.symbol

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
setwd('./')


scdata_tran <- readRDS('./scdata_tpm.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))
data_plot_tran <- data_tran %>% 
  rownames_to_column('Genes') %>%  
  reshape2::melt(
    id.vars = 'Genes',
    measure.vars = colnames(.)[2:ncol(.)],
    variable.name = 'sample',
    value.name = 'copy_num'
  ) %>% 
  filter(copy_num>0) %>%
  left_join(meta_data_tran,by = 'sample') %>% 
  mutate(
    log10_cp = log10(copy_num),
    group_omics = 'transcription'
  ) %>% 
  dplyr::select(
    Genes,sample,copy_num,group_cell,
    group_stage,group_sample,log10_cp, group_omics
  )
data_plot_trans_Oo <- data_plot_tran %>% 
  filter(group_cell == 'Oocyte')

data_plot <- rbind(
  data_plot_trans_Oo %>% 
    mutate(group = 'All'),
  data_plot_trans_Oo %>% 
    filter(group_cell == 'Oocyte',Genes %in% data_tf_list) %>% 
    mutate(group = 'TFs'),
  data_plot_trans_Oo %>% 
    filter(group_cell == 'Oocyte',Genes %in% data_sm_list) %>% 
    mutate(group = 'SMs')
)

ggplot() +   
  geom_histogram(
    data = data_plot_trans_Oo,
    mapping = aes(x = log10_cp,y = ..density..),
    bins = 100, alpha = 0.3,color = 'white',
    position = 'identity'
  ) +
  geom_histogram(#geom_density
    data = data_plot, 
    mapping = aes(x = log10_cp,color = group,y = ..density..),
    bins = 100,fill = 'white',alpha = 0.8,linewidth = 1,key_glyph = "path"
  ) +
  geom_vline(
    xintercept = mean(data_plot_trans_Oo$log10_cp),
    color = 'black',linetype = 'longdash',
    linewidth = 2
  ) +
  annotate(
    geom = 'text',
    x = (mean(data_plot_trans_Oo$log10_cp)-0.75),
    y= 0.02,label = paste('Mean:',mean(data_plot_trans_Oo$log10_cp) %>% round(.,digits = 2)),
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
    legend.position = c(0.1,0.9),
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
  filename = './fig1e_trans_Oocyte.png',
  width = 10,height = 9
)
data_plot_tran_Granulosa <- data_plot_tran %>% 
  filter(group_cell == 'Granulosa')

data_plot <- rbind(
  data_plot_tran_Granulosa %>% 
    mutate(group = 'All'),
  data_plot_tran_Granulosa %>% 
    filter(group_cell == 'Granulosa',Genes %in% data_tf_list) %>% 
    mutate(group = 'TFs'),
  data_plot_tran_Granulosa %>% 
    filter(group_cell == 'Granulosa',Genes %in% data_sm_list) %>% 
    mutate(group = 'SMs')
)

ggplot() +   
  geom_histogram(
    data = data_plot_tran_Granulosa,
    mapping = aes(x = log10_cp,y = ..density..),
    bins = 100, alpha = 0.3,color = 'white',
    position = 'identity'
  ) +
  geom_histogram(#geom_density
    data = data_plot, 
    mapping = aes(x = log10_cp,color = group,y = ..density..),
    bins = 100,fill = 'white',alpha = 0.8,linewidth = 1,key_glyph = "path"
  ) +
  geom_vline(
    xintercept = median(data_plot_tran_Granulosa$log10_cp),
    color = 'black',linetype = 'longdash',linewidth = 2
  ) +
  annotate(
    geom = 'text',
    x = (median(data_plot_tran_Granulosa$log10_cp)-0.8),
    y= 0.04,label = paste('Median:',median(data_plot_tran_Granulosa$log10_cp) %>% round(.,digits = 2)),
    size = 6
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
    legend.position = c(0.1,0.9),
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
  filename = './fig1e_trans_Granulosa.png',
  width = 10,height = 9
)

