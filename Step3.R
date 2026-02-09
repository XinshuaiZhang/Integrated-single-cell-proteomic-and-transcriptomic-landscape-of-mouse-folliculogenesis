rm(list = ls())

load('./protein_gc.RData')
names(acore.list) <- 1:16
tmp <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('group','Accessions','MEM.SHIP')) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  dplyr::select(-'Accessions') %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
colnames(data_plot)
data_acorelist <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('cluster','Accessions','MEM.SHIP'))
data_plot_use <- data_plot %>% 
  rename('Accessions' = Gene_symbol) %>% 
  mutate(cluster = cluster %>% as.character()) %>% 
  left_join(
    data_acorelist,
    by = c('cluster','Accessions')
  ) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup() %>% 
  dplyr::select(-c('Accessions')) %>% 
  rename('Gene_symbol' = 'Genes')


res_mfuzz_pro <- read.csv('./res_mfuzz_protein_gc.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_gc.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","5","8","11","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("6",'11','13','15','16'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3",'5',"10"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_cluster_raw <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% stack() %>% 
  rename_all(~c('Gene_symbol','cluster'))
data_plot_use <- data_plot_use %>% 
  dplyr::select(-c(cluster)) %>% 
  left_join(data_cluster_raw,by = 'Gene_symbol',relationship ="many-to-many")
data_plot_use$cluster %>% table()


plot_single_cluster=function(plottype = 'second',data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster %in% cluster,]
  if(plottype == 'first'){
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
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }else if(plottype == 'second'){
    p <- ggplot(data_plot, aes(x=time, y=(values))) + 
      geom_line(aes(group=Gene_symbol),color = 'grey90',linewidth=1,alpha=1) +
      geom_smooth(
        aes(x=as.numeric(time), y=values),linewidth=2,
        method = "loess",formula = 'y ~ x', span = 0.8,se = FALSE, color = '#486789') +
      scale_linetype_manual(values = 'dashed') + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }
  return(p)
}

path <- './figure4_cluster_lineplot_GC'
dir.create(path,showWarnings = FALSE)
for(cluster_use in (data_plot_use$cluster %>% levels())){
  print(cluster_use)
  path_save <- paste(path,'/',cluster_use,'.pdf',sep = '')
  pdf(path_save,width = 8,height = 7)
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  dev.off()
  path_save <- paste(path,'/',cluster_use,'.png',sep = '')
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  ggsave(path_save,width = 9,height = 7)
}
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:7){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)#,heights = 1
  return(p)
}
path_save <- paste(path,'/','cluster_all','.png',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)
path_save <- paste(path,'/','cluster_all','.pdf',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)

plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 14),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:12){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,16,1)[-seq(1,16,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
path <- './figures4'
dir.create(path,showWarnings = FALSE)
path_save <- paste(path,'/','cluster_all_Granulosa_protein','.pdf',sep = '')
ggsave(plot = p_mfuzz,width = 14,height = 8,dpi = 300,filename = path_save)


library(stringr)
library(tibble)
library(tidyr)

load('./protein_o.RData')
names(acore.list) <- 1:16
tmp <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('group','Accessions','MEM.SHIP')) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  dplyr::select(-'Accessions') %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
colnames(data_plot)
data_acorelist <- acore.list %>% 
  do.call(rbind,.) %>% 
  tibble::rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('cluster','Accessions','MEM.SHIP'))
data_plot_use <- data_plot %>% 
  dplyr::rename('Accessions' = Gene_symbol) %>% 
  mutate(cluster = cluster %>% as.character()) %>% 
  left_join(
    data_acorelist,
    by = c('cluster','Accessions')
  ) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup() %>% 
  dplyr::select(-c('Accessions')) %>% 
  mutate(
    cluster = case_when(
      cluster %in% c("1","3","6","7","9","10","11","13") ~ 'pro_inc',
      cluster %in% c('12') ~ 'pro_dec'
    )
  ) %>% 
  dplyr::rename('Gene_symbol' = 'Genes')


res_mfuzz_pro <- read.csv('./res_mfuzz_protein_o.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_o.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","6","7","9","10","11","13"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('4'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3","8","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("1","4","5","11","12","15","16"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_cluster_raw <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% stack() %>% 
  rename_all(~c('Gene_symbol','cluster'))
data_plot_use <- data_plot_use %>% 
  dplyr::select(-c(cluster)) %>% 
  left_join(data_cluster_raw,by = 'Gene_symbol',relationship ="many-to-many")
data_plot_use$cluster %>% table()


plot_single_cluster=function(plottype = 'second',data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster %in% cluster,]
  if(plottype == 'first'){
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
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }else if(plottype == 'second'){
    p <- ggplot(data_plot, aes(x=time, y=(values))) + 
      geom_line(aes(group=Gene_symbol),color = 'grey90',linewidth=1,alpha=1) +
      geom_smooth(
        aes(x=as.numeric(time), y=values),linewidth=2,
        method = "loess",formula = 'y ~ x', span = 0.8,se = FALSE, color = '#486789') +
      scale_linetype_manual(values = 'dashed') + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }
  return(p)
}

path <- './figure4_cluster_lineplot'
dir.create(path,showWarnings = FALSE)
for(cluster_use in (data_plot_use$cluster %>% levels())){
  print(cluster_use)
  path_save <- paste(path,'/',cluster_use,'.pdf',sep = '')
  pdf(path_save,width = 8,height = 7)
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  dev.off()
  path_save <- paste(path,'/',cluster_use,'.png',sep = '')
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  ggsave(path_save,width = 9,height = 7)
}
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:7){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)#,heights = 1
  return(p)
}
path_save <- paste(path,'/','cluster_all','.png',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)
path_save <- paste(path,'/','cluster_all','.pdf',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 14),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:12){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,16,1)[-seq(1,16,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
path <- './figures4'
dir.create(path,showWarnings = FALSE)
path_save <- paste(path,'/','cluster_all_Oocyte_protein','.pdf',sep = '')
ggsave(plot = p_mfuzz,width = 14,height = 8,dpi = 300,filename = path_save)


rm(list = ls())

load('./trans_gc.RData')
names(acore.list) <- 1:16
tmp <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('group','Genes','MEM.SHIP')) %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
colnames(data_plot)
data_acorelist <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('cluster','Genes','MEM.SHIP'))
data_plot_use <- data_plot %>% 
  rename('Genes' = Gene_symbol) %>% 
  mutate(cluster = cluster %>% as.character()) %>% 
  left_join(
    data_acorelist,
    by = c('cluster','Genes')
  ) %>%  
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup() %>%
  rename('Gene_symbol' = 'Genes')


res_mfuzz_pro <- read.csv('./res_mfuzz_protein_gc.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_gc.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","5","8","11","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("6",'11','13','15','16'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3",'5',"10"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)

data_cluster_raw <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% stack() %>% 
  rename_all(~c('Gene_symbol','cluster'))
data_plot_use <- data_plot_use %>% 
  dplyr::select(-c(cluster)) %>% 
  left_join(data_cluster_raw,by = 'Gene_symbol',relationship ="many-to-many")
data_plot_use$cluster %>% table()


plot_single_cluster=function(plottype = 'second',data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster %in% cluster,]
  if(plottype == 'first'){
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
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }else if(plottype == 'second'){
    p <- ggplot(data_plot, aes(x=time, y=(values))) + 
      geom_line(aes(group=Gene_symbol),color = 'grey90',linewidth=1,alpha=1) +
      geom_smooth(
        aes(x=as.numeric(time), y=values),linewidth=2,
        method = "loess",formula = 'y ~ x', span = 0.8,se = FALSE, color = '#486789') +
      scale_linetype_manual(values = 'dashed') + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }
  return(p)
}

path <- './figure4_cluster_lineplot_trans_GC'
dir.create(path,showWarnings = FALSE)
for(cluster_use in (data_plot_use$cluster %>% levels())){
  print(cluster_use)
  path_save <- paste(path,'/',cluster_use,'.pdf',sep = '')
  pdf(path_save,width = 8,height = 7)
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  dev.off()
  path_save <- paste(path,'/',cluster_use,'.png',sep = '')
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  ggsave(path_save,width = 9,height = 7)
}
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:7){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)#,heights = 1
  return(p)
}
path_save <- paste(path,'/','cluster_all','.png',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)
path_save <- paste(path,'/','cluster_all','.pdf',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)

plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
  colo <- colorRampPalette(rev(c("#305E95","grey92")))(30)
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
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 14),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:12){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,16,1)[-seq(1,16,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
path <- './figures4'
dir.create(path,showWarnings = FALSE)
path_save <- paste(path,'/','cluster_all_Granulosa_transcript','.pdf',sep = '')
ggsave(plot = p_mfuzz,width = 14,height = 8,dpi = 300,filename = path_save)

rm(list = ls())

load('./trans_o.RData')
names(acore.list) <- 1:16
tmp <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('group','Genes','MEM.SHIP')) %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
colnames(data_plot)
data_acorelist <- acore.list %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('group') %>% 
  mutate(group = group %>% str_split('\\.',simplify = TRUE) %>% .[,1]) %>% 
  rename_all(~c('cluster','Genes','MEM.SHIP'))
data_plot_use <- data_plot %>% 
  rename('Genes' = Gene_symbol) %>% 
  mutate(cluster = cluster %>% as.character()) %>% 
  left_join(
    data_acorelist,
    by = c('cluster','Genes')
  ) %>%  
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup() %>%
  rename('Gene_symbol' = 'Genes')


res_mfuzz_pro <- read.csv('./res_mfuzz_protein_o.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_o.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","6","7","9","10","11","13"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('4'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3","8","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("1","4","5","11","12","15","16"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_cluster_raw <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% stack() %>% 
  rename_all(~c('Gene_symbol','cluster'))
data_plot_use <- data_plot_use %>% 
  dplyr::select(-c(cluster)) %>% 
  left_join(data_cluster_raw,by = 'Gene_symbol',relationship ="many-to-many")
data_plot_use$cluster %>% table()


plot_single_cluster=function(plottype = 'second',data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster %in% cluster,]
  if(plottype == 'first'){
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
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }else if(plottype == 'second'){
    p <- ggplot(data_plot, aes(x=time, y=(values))) + 
      geom_line(aes(group=Gene_symbol),color = 'grey90',linewidth=1,alpha=1) +
      geom_smooth(
        aes(x=as.numeric(time), y=values),linewidth=2,
        method = "loess",formula = 'y ~ x', span = 0.8,se = FALSE, color = '#486789') +
      scale_linetype_manual(values = 'dashed') + 
      theme(panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
            plot.title = element_text(hjust = 0.5,size = 22),
            legend.position = 'none')+
      labs(x = 'Stage',y = 'Relative Expression',title = paste('Cluster ',cluster,sep = ''))
  }
  return(p)
}

path <- './figure4_cluster_lineplot_trans_Oocyte'
dir.create(path,showWarnings = FALSE)
for(cluster_use in (data_plot_use$cluster %>% levels())){
  print(cluster_use)
  path_save <- paste(path,'/',cluster_use,'.pdf',sep = '')
  pdf(path_save,width = 8,height = 7)
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  dev.off()
  path_save <- paste(path,'/',cluster_use,'.png',sep = '')
  plot_single_cluster(data_plot=data_plot_use,cluster = cluster_use)
  ggsave(path_save,width = 9,height = 7)
}
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:7){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)#,heights = 1
  return(p)
}
path_save <- paste(path,'/','cluster_all','.png',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)
path_save <- paste(path,'/','cluster_all','.pdf',sep = '')
plot_mfuzz(data_plot_use)
ggsave(path_save,width = 9,height = 38)

plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
  colo <- colorRampPalette(rev(c("#305E95","grey92")))(30)
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
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1,size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 14),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:12){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,16,1)[-seq(1,16,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
path <- './figures4'
dir.create(path,showWarnings = FALSE)
path_save <- paste(path,'/','cluster_all_Oocyte_transcript','.pdf',sep = '')
ggsave(plot = p_mfuzz,width = 14,height = 8,dpi = 300,filename = path_save)


library(colorRamp2)
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
replace_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  

  non_outliers <- x[x >= lower_bound & x <= upper_bound]
  tmp <- x[x <= lower_bound | x >= upper_bound]
  if (length(non_outliers) > 0) {

    x[x > upper_bound] <- max(non_outliers, na.rm = TRUE)
  }
  return(x)
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
write.csv(data_singlestage,'./data_exp_gene_pro.csv')
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
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))

histone_data<-read.csv('./histone_gene_use.csv',sep=",",header = T,check.names = F) %>% 
  pull(names(.)[2])



gene_protein_coding <-
  read.csv('./gene_protein_coding.csv', row.names = 1)
head(gene_protein_coding)



res_mfuzz_pro <- read.csv('./res_mfuzz_protein_gc.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_gc.csv',row.names = 1)
value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","5","8","11","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("6",'11','13','15','16'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3",'5',"10"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
)
data_gene_use <- stack(data_res) %>% 
  rename_all(~c('Genes','group')) %>% 
  mutate(group = factor(group,levels = names(data_res)))
data_gene_use$group %>% table()



scdata_pro$group_cell %>% unique()
scdata_pro_sub <- subset(scdata_pro,subset = group_cell == 'Granulosa')
data_pro_use <- FetchData(object = scdata_pro_sub,vars = c('group_stage_cell',data_gene_use$Genes))
data_plot <- data_pro_use %>% 
  group_by(group_stage_cell) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  column_to_rownames('group_stage_cell') %>% 
  dplyr::filter(rowSums(.) > 0) %>%
  {log10(.+1)} %>% {scale(.)} %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('Genes') %>% 
  tidyr::complete(Genes = data_gene_use$Genes, fill = list(cluster = NA)) %>%
  mutate(Genes = factor(Genes,levels = data_gene_use$Genes[data_gene_use$Genes %in% Genes])) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
min_value <- min(data_plot,na.rm = TRUE)
data_plot <- data_plot %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
gene_navalue <- data_plot %>% filter(is.na(rowSums(.))) %>% rownames()
data_plot[rownames(data_plot) == 'Lsm11',]
data_plot[is.na(data_plot)] <- min(data_plot,na.rm = TRUE)
data_protein_use <- data_gene_use %>% 
  filter(Genes %in% rownames(data_plot))
library(ComplexHeatmap)
color_group_stage_cell <- c(
  '#F8B620','#FF800E','#ED444A','#D14B70'
)
group <- data.frame(
  group_stage_cell = colnames(data_plot) %>% str_remove(': Granulosa') %>% factor(.,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )),
  rowname = colnames(data_plot)
) %>% 
  column_to_rownames('rowname')
names(color_group_stage_cell) <- group$group_stage_cell
topanno = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  show_legend = FALSE,
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
draw(topanno)
group_row <- data.frame(
  group = data_protein_use$group,
  rowname = rownames(data_plot)
) %>% 
  column_to_rownames('rowname')
leftanno = rowAnnotation(
  df = group_row, # Row annotation
  border = TRUE,
  show_annotation_name = FALSE,show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 24),
  simple_anno_size = unit(20, "points"),
  gap = unit(12, "points")
)
draw(leftanno)
ht <- Heatmap(
  matrix = data_plot %>% as.matrix(),
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
ht = draw(ht)

gene_order <- row_order(ht) %>% 
  stack()#ht@row_names_param[["labels"]]
gene_order <- rownames(data_plot) %>% 
  .[gene_order$values]
data_plot_order <- data_plot %>% 
  rownames_to_column('Genes') %>% 
  mutate(
    Genes = factor(Genes,levels = gene_order)
  ) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
data_plot_order[rownames(data_plot_order) %in% gene_navalue,] <- NA
col_fun_pro <- colorRamp2(
  breaks = c(min(data_plot_order,na.rm = TRUE), 0, max(data_plot_order,na.rm = TRUE)),
  colors = c("#716CAC", "white", "#c2473b")
)
ht <- Heatmap(
  matrix = data_plot_order %>% as.matrix(),
  col = col_fun_pro,
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = 'Protein',
  column_title_gp = gpar(fontsize = 14, fontface = "bold", hjust = 1),
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht,
  heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,
  padding = unit(c(10, 10, 10, 30), "mm"))


scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Granulosa')
data_tran_use <- FetchData(object = scdata_tran_sub,layer = 'counts',vars = c('group_stage_cell',data_gene_use$Genes))
data_plot_tran <- data_tran_use %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    "Secondary_Granulosa","Early antral_Granulosa","Antral_Granulosa",'Preovulatory_Granulosa'
  ))) %>% 
  group_by(group_stage_cell) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  column_to_rownames('group_stage_cell') %>% 
  {log10(.+1)} %>%  {scale(.)} %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  tidyr::complete(Genes = data_gene_use$Genes, fill = list(cluster = NA)) %>%
  mutate(Genes = factor(Genes,levels = gene_order)) %>% 
  arrange(Genes) %>% column_to_rownames('Genes')
min_value <- min(data_plot_tran,na.rm = TRUE)
data_plot_tran <- data_plot_tran %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
data_tran_use <- data_gene_use %>%
  filter(Genes %in% rownames(data_plot_tran))
topanno_tran = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
col_fun_tran <- colorRamp2(
  breaks = c(min(data_plot_tran), 0, max(data_plot_tran)),
  colors = c("#0D7291", "white", "#FC7136")
)
ht_tran <- Heatmap(
  matrix = data_plot_tran %>% as.matrix(),
  col = col_fun_tran,
  top_annotation = topanno_tran,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_tran_use$group,
  row_title = NULL,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht_tran,heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))





ht_legend <- Legend(
  col_fun = col_fun_pro, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_legend)
ht_tran_legend <- Legend(
  col_fun = col_fun_tran, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_tran_legend)
ht_grob1 <- grid.grabExpr(draw(
  ht,
  heatmap_legend_side = NULL
))
ht_grob2 <- grid.grabExpr(draw(
  ht_tran,
  annotation_legend_list = list(ht_legend,ht_tran_legend), 
  merge_legends = TRUE, 
  column_title = "Transcript",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", just = "center")
))


library(gridExtra)
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))

pdf('./cluster_of_8type_heatmap_GC.pdf')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()

pdf('./cluster_of_8type_heatmap_GC.pdf')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()

library(colorRamp2)
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
replace_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  

  non_outliers <- x[x >= lower_bound & x <= upper_bound]
  tmp <- x[x <= lower_bound | x >= upper_bound]
  if (length(non_outliers) > 0) {

    x[x > upper_bound] <- max(non_outliers, na.rm = TRUE)
  }
  return(x)
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

scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
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



res_mfuzz_pro <- read.csv('./res_mfuzz_protein_o.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_o.csv',row.names = 1)
value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","6","7","9","10","11","13"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3","8","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("1","4","5","11","12","15","16"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
)
data_gene_use <- stack(data_res) %>% 
  rename_all(~c('Genes','group')) %>% 
  mutate(group = factor(group,levels = names(data_res)))
data_gene_use$group %>% table()



scdata_pro$group_stage_cell
scdata_pro_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte')
data_pro_use <- FetchData(object = scdata_pro_sub,vars = c('group_stage_cell',data_gene_use$Genes))
data_plot <- data_pro_use %>% 
  group_by(group_stage_cell) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  column_to_rownames('group_stage_cell') %>% 
  dplyr::filter(rowSums(.) > 0) %>%
  {log10(.+1)} %>% {scale(.)} %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('Genes') %>% 
  tidyr::complete(Genes = data_gene_use$Genes, fill = list(cluster = NA)) %>%
  mutate(Genes = factor(Genes,levels = data_gene_use$Genes[data_gene_use$Genes %in% Genes])) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
min_value <- min(data_plot,na.rm = TRUE)
data_plot <- data_plot %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
gene_navalue <- data_plot %>% filter(is.na(rowSums(.))) %>% rownames()
data_plot[rownames(data_plot) == 'Lsm11',]
data_plot[is.na(data_plot)] <- min(data_plot,na.rm = TRUE)
data_protein_use <- data_gene_use %>% 
  filter(Genes %in% rownames(data_plot))
library(ComplexHeatmap)
color_group_stage_cell <- c(
  '#F8B620','#FF800E','#ED444A','#D14B70'
)
group <- data.frame(
  group_stage_cell = colnames(data_plot) %>% str_remove(': Oocyte') %>% factor(.,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )),
  rowname = colnames(data_plot)
) %>% 
  column_to_rownames('rowname')
names(color_group_stage_cell) <- group$group_stage_cell
topanno = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  show_legend = FALSE,
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
draw(topanno)
group_row <- data.frame(
  group = data_protein_use$group,
  rowname = rownames(data_plot)
) %>% 
  column_to_rownames('rowname')
leftanno = rowAnnotation(
  df = group_row, # Row annotation
  border = TRUE,
  show_annotation_name = FALSE,show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 24),
  simple_anno_size = unit(20, "points"),
  gap = unit(12, "points")
)
draw(leftanno)
ht <- Heatmap(
  matrix = data_plot %>% as.matrix(),
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
ht = draw(ht)

gene_order <- row_order(ht) %>% 
  stack()#ht@row_names_param[["labels"]]
gene_order <- rownames(data_plot) %>% 
  .[gene_order$values]
data_plot_order <- data_plot %>% 
  rownames_to_column('Genes') %>% 
  mutate(
    Genes = factor(Genes,levels = gene_order)
  ) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
data_plot_order[rownames(data_plot_order) %in% gene_navalue,] <- NA
col_fun_pro <- colorRamp2(
  breaks = c(min(data_plot_order,na.rm = TRUE), 0, max(data_plot_order,na.rm = TRUE)),
  colors = c("#716CAC", "white", "#c2473b")
)
ht <- Heatmap(
  matrix = data_plot_order %>% as.matrix(),
  col = col_fun_pro,
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = 'Protein',
  column_title_gp = gpar(fontsize = 14, fontface = "bold", hjust = 1),
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
ht
draw(
  ht,
  heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))


scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte')
data_tran_use <- FetchData(object = scdata_tran_sub,layer = 'counts',vars = c('group_stage_cell',data_gene_use$Genes))
data_plot_tran <- data_tran_use %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = c(
    "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte",'Preovulatory_Oocyte'
  ))) %>% 
  group_by(group_stage_cell) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  column_to_rownames('group_stage_cell') %>% 
  {log10(.+1)} %>%  {scale(.)} %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  tidyr::complete(Genes = data_gene_use$Genes, fill = list(cluster = NA)) %>%
  mutate(Genes = factor(Genes,levels = gene_order)) %>% 
  arrange(Genes) %>% column_to_rownames('Genes')
min_value <- min(data_plot_tran,na.rm = TRUE)
data_plot_tran <- data_plot_tran %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
data_tran_use <- data_gene_use %>%
  filter(Genes %in% rownames(data_plot_tran))
topanno_tran = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
col_fun_tran <- colorRamp2(
  breaks = c(min(data_plot_tran), 0, max(data_plot_tran)),
  colors = c("#0D7291", "white", "#FC7136")
)
ht_tran <- Heatmap(
  matrix = data_plot_tran %>% as.matrix(),
  col = col_fun_tran,
  top_annotation = topanno_tran,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_tran_use$group,
  row_title = NULL,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht_tran,heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))




ht_legend <- Legend(
  col_fun = col_fun_pro, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_legend)
ht_grob1 <- grid.grabExpr(draw(
  ht,
  heatmap_legend_side = NULL
))
ht_grob2 <- grid.grabExpr(draw(
  ht_tran,
  annotation_legend_list = list(ht_legend), 
  merge_legends = TRUE, 
  column_title = "Transcript",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", just = "center")
))


library(gridExtra)
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))

pdf('./cluster_of_8type_heatmap_Oocyte.pdf')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()
pdf('./cluster_of_8type_heatmap_Oocyte.pdf')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()

data_raw <- read.csv('./figure4_c.csv')
tmp <- lapply(1:nrow(data_raw),FUN = function(x){
  tmp <- data_raw$geneID[x] %>% str_split('/') %>% unlist()
})
names(tmp) <- paste(data_raw$Description,data_raw$group,sep = '--')
library(tidyr)
data_filter <- stack(tmp) %>% 
  rename_all(~c('Genes','pathwayname')) %>% 
  separate(col = pathwayname, into = c("pathwayname", "group"), sep = "--") %>% 
  group_by(Genes,group) %>% 
  filter(n()==1) %>% 
  group_by(pathwayname,group) %>% 
  reframe(counts = n(),gene_use = paste(Genes,collapse = ';'),.groups = 'drop') %>% 
  left_join(data_raw %>% rename(pathwayname = Description),by = c('pathwayname','group')) %>% 
  arrange(group)
write.csv(data_filter,'./figure4_c_filter.csv')
data_filter_gene <- stack(tmp) %>% 
  rename_all(~c('Genes','pathwayname')) %>% 
  separate(col = pathwayname, into = c("pathwayname", "group"), sep = "--") %>% 
  group_by(Genes,group) %>% 
  reframe(pathway_counts = n(),pathway = paste(pathwayname,collapse = ';')) %>% 
  arrange(group,desc(pathway_counts),Genes)
write.csv(data_filter_gene,'./figure4_c_filter_gene.csv')

library(statmod)
library(ggplot2)
library(tibble)
library(stringr)
library(dplyr)
library(reshape2)
library(stringi)
library(Seurat)
library(RColorBrewer)
library(pheatmap)
library(monocle3)
library(Seurat)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
)


data_monoce2_loc <- read.csv('./monocle2_plot_data.csv',row.names = 1) %>% 
  mutate(
    group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary: Oocyte","Early antral: Oocyte",
      "Antral: Oocyte","Preovulatory: Oocyte"
    ))
  )
colnames(data_monoce2_loc)

p <- ggplot(
  data = data_monoce2_loc,
  mapping = aes(x = data_dim_1, y = data_dim_2)
  ) +
  geom_smooth(
    method = "loess", 
    se = FALSE,
    size = 3,color = 'lightslategrey',
    ) +
  geom_point(aes(color = group_stage_cell),size = 4) +
  guides(color = guide_legend(title = 'Stage')) +
  scale_color_manual(values = c('#FFAD72','#F76D5E','#D82632','#A50021')) +
  labs(x = 'Component 1',y = 'Component 2') +
  theme_classic()
p

smooth_data <- ggplot_build(p)$data[[1]]
start_point <- smooth_data[1, c("x", "y")]
second_point <- smooth_data[2, c("x", "y")]


p + annotate(
  "segment",
  x = second_point$x,
  y = second_point$y,
  xend = start_point$x,
  yend = start_point$y,
  arrow = arrow(type = "closed", length = unit(0.2, "inches")),
  color = 'lightslategrey',
  size = 1.5
)
ggsave(
  filename = './result_0725_Fig1&4_/',
  width = 10,height = 6
)

library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(tibble)

data_scale <- read.csv('./data_copy_num_sorted.csv',row.names = 1)
data_gene_protein <- data_scale[,c('ProteinGroups','Genes')] %>% 
  dplyr::rename(Accessions = ProteinGroups)


data_p <- read.csv(
  './data_p_kruskal_test_Oocyte.csv',
  row.names = 1) %>% 
  dplyr::select(c('Genes','pvalue_Permutation_pro')) %>% 
  dplyr::filter(pvalue_Permutation_pro<0.05)
gene_select_p <- data_gene_protein %>% 
  dplyr::filter(Genes %in% data_p$Genes) %>% 
  pull(Accessions)


scdata <- readRDS('./scdata_use_filter_outliers.rds')

FindVariableFeatures(scdata)
scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata <- subset(scdata,subset = group_cell == 'Oocyte')
tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rowwise() %>%
  summarise(
    num_positive = sum(c_across(everything()) > 0)
    ) %>% 
  mutate(
    protein_group = rownames(scdata)
  )
gene_select <- tmp$protein_group[tmp$num_positive>2]

rownames(scdata) %>% length()
scdata <- subset(scdata,features = gene_select)
scdata <- subset(scdata,features = gene_select_p)
rownames(scdata) %>% length()
scdata@meta.data %>% colnames()
scdata_sce <- as.SingleCellExperiment(scdata)
scdata_sce <- prepSCE(scdata_sce,
                      kid = "group_cell",
                      gid = "group_stage",#group——id，
                      sid = "sample",#sample_id, 
                      drop=T)
pb <- aggregateData(scdata_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()

boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 16
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')
write.csv(data_plot,'./data_plot_protein.csv')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
p_mfuzz
ggsave(width = 16,height = 8,dpi = 300,filename = './fig_all.png')


data_plot_gene <- data_plot %>% 
  as.data.frame() %>% 
  dplyr::filter(cluster == 5)
gene_use <- acore.list[[5]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:50) %>% pull(NAME)
library(tibble)
data_tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column(var = 'Accessions') %>% 
  filter(Accessions %in% gene_use) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  dplyr::select(-c(Accessions)) %>% 
  column_to_rownames(var = 'Genes')
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  mutate(group_stage_cell = factor(
    group_stage_cell,levels = c(
      "Secondary_Granulosa", "Early antral_Granulosa",
      "Antral_Granulosa", "Preovulatory_Granulosa"
    )
  )) %>%
  pull(group_stage_cell) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  pull(group_stage) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
annotation_col = data.frame(
  Group = scdata$group_stage
)
rownames(annotation_col) = colnames(data_tmp)
ann_colors = list(Group = c(
  'Secondary' = '#C34F73',
  'Early antral' = '#4B8600',
  'Antral' = '#4DBBD5',
  'Preovulatory' = '#3C5488'
))
p_pheatmap <- pheatmap::pheatmap(
  data_tmp, 
  scale = 'row',
  display_numbers = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  gaps_col = break_heatmap,
  annotation_names_row = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  annotation_colors = ann_colors,
  na_col = 'white',
  border_color = 'white',
  fontsize = 20, fontsize_row = 10,
  angle_col = 45,
  main = 'cluster 5(Top 50 genes)'
)
p_pheatmap_ggplot <- cowplot::plot_grid(p_pheatmap$gtable)
p_pheatmap_ggplot
ggsave(width = 16,height = 10,bg = 'white',filename = './p_pheatmap_5.png')


scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
gene_select <- acore.list[[11]] %>% 
  filter(MEM.SHIP>0.7) %>% 
  filter(!(grepl('Gm',NAME))) %>% 
  pull(NAME)
DotPlot(
  object = scdata,
  scale = TRUE,
  features = gene_select,
  group.by = 'group_stage'
) +
  coord_flip() +
  scale_color_gradient(low = '#DAE5F0FF',high = '#26456EFF') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.text = element_text(size = 20)
  )
ProNames <- acore.list[[11]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:12) %>% pull(NAME)
data_plot <- FetchData(
  object = scdata,layer = 'counts',
  clean = 'none',vars = c(ProNames,'group_stage_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_stage_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  as.data.frame() %>%
  mutate(gene = protein)

ggplot(data = data_plot,mapping = aes(x = group_stage_cell,y = expression,fill = group_stage_cell)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 4,scales = 'free_y',strip.position = 'left') +
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



library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
acore.list <- acore(dat,cl=cl,min.acore=0.2)
genes<-acore.list[[5]] %>% 
  rename(Accessions = NAME) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  pull(Genes)
gene<-str_to_title(genes)
gene
keys(org.Mm.eg.db, keytype="SYMBOL") %>% head()

gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
p_go <-
  dotplot(ego,
          showCategory = 15,
          font.size = 20)
ggsave("./GO_Cluster5.png", dpi=300, width=6, height = 6, units = "in")
barplot(ego, showCategory=20,title="EnrichmentGO")
data_go <- ego@result
write.csv(data_go,'./data_go_cluster5.csv')
p_tmp_1 <- cowplot::plot_grid(
  p_pheatmap_ggplot,
  p_go,
  rel_widths = c(1.2, 1),
  axis = 't',
  align = 'v'
)
p_res <- cowplot::plot_grid(
  p_mfuzz,
  p_tmp_1,ncol = 1,
  rel_heights = c(1, 1.5),
  align = 'v'
)
ggsave(
  plot = p_res,
  filename = './p_res_pro_o.png',
  width = 20,bg = 'white',
  height = 24
)


save.image('./protein_o.RData')
saveRDS(cl,'./protein_o_cl.rds')
saveRDS(dat,'./protein_o_dat.rds')
res_mfuzz <- lapply(1:16, function(i){
  data_mfuzz <- acore.list[[i]] %>% 
    mutate(group = i)
}) %>% do.call(rbind,.)  %>% 
  rename(Accessions = NAME) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  left_join(data_p,by = 'Genes')
write.csv(res_mfuzz,'./res_mfuzz_protein_o.csv')


acore.list <- acore(dat,cl=cl,min.acore=0.2)
dir.create('./res')
for(i in 1:length(acore.list)){
  library(clusterProfiler)
  print(i)
  genes<-acore.list[[i]] %>% 
    rename(Accessions = NAME) %>% 
    left_join(data_gene_protein,by = 'Accessions') %>% 
    pull(Genes)
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  file_name <- paste(
    "./GO_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")

  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}

acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot_pro <- tmp %>% 
  filter(
    cluster %in% c(5,9)
  ) %>% 
  dplyr::rename(Accessions = 'Gene_symbol') %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  dplyr::select(-c('Accessions')) %>% 
  dplyr::rename('Gene_symbol' = 'Genes') %>% 
  group_by(Gene_symbol) %>% 
  summarise(
    across(names(.)[c(1,2,3,4)],.fns = sum),
    cluster = first(cluster))
data_plot_pro$cluster %>% table()


ha = HeatmapAnnotation(
  bar = data_plot_pro %>% colnames() %>% .[2:5],
  col = list(
    bar = c(
      "Secondary" = "#C34F73","Early antral" = "#0088D4", 
      "Antral" = "#4B8600", 'Preovulatory'= '#5F559B')
    ),
  show_annotation_name = FALSE
  )
ha = HeatmapAnnotation(
  foo = anno_block(
    gp = gpar(fill = c('#C34F73','#0088D4','#4B8600','#5F559B')),
    labels = data_plot_pro %>% colnames() %>% .[2:5],
    labels_gp = gpar(col = "white", fontsize = 10)
    )
  )
p1 <- Heatmap(
  data_plot_pro %>% 
    select(-c('Gene_symbol','cluster')) %>% 
    as.matrix(),
  name = "mat",
  top_annotation = ha,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  row_split = data_plot_pro$cluster %>% as.character(),
  column_split = 1:4,
  column_title = NULL,
  row_title = NULL,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  show_heatmap_legend = FALSE
)
p1


data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')
data_plot_gene <- data_plot %>% 
  as.data.frame() %>% 
  dplyr::filter(cluster %in% c(5,9))
gene_use <- c(
  acore.list[[5]] %>%
    arrange(desc(MEM.SHIP)) %>%
    pull(NAME),
  acore.list[[9]] %>%
    arrange(desc(MEM.SHIP)) %>%
    pull(NAME)
) %>% unique()

library(tibble)
data_tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column(var = 'Accessions') %>% 
  filter(Accessions %in% gene_use) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  mutate(
    Accessions = factor(Accessions,levels = gene_use %>% unique())
  ) %>% 
  arrange(Accessions) %>% 
  mutate(
    Genes = factor(Genes,levels = Genes %>% unique())
  ) %>% 
  dplyr::select(-c(Accessions)) %>% 
  group_by(Genes) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  column_to_rownames(var = 'Genes')
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  mutate(group_stage_cell = factor(
    group_stage_cell,levels = c(
      "Secondary_Granulosa", "Early antral_Granulosa",
      "Antral_Granulosa", "Preovulatory_Granulosa"
    )
  )) %>%
  pull(group_stage_cell) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  pull(group_stage) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
annotation_col = data.frame(
  Group = scdata$group_stage
)
rownames(annotation_col) = colnames(data_tmp)
ann_colors = list(Group = c(
  'Secondary' = '#C34F73',
  'Early antral' = '#4B8600',
  'Antral' = '#4DBBD5',
  'Preovulatory' = '#3C5488'
))
p_protein <- pheatmap::pheatmap(
  data_tmp, 
  scale = 'row',
  display_numbers = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  gaps_col = break_heatmap,
  annotation_names_row = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_legend = FALSE,
  annotation_colors = ann_colors,
  legend = FALSE,
  na_col = 'white',
  border_color = 'white',
  fontsize = 20, 
  fontsize_row = 10,
  angle_col = 45,
  main = 'Protein'
)
p_protein_ggplot <- cowplot::plot_grid(p_protein$gtable)
p_protein_ggplot



scdata_tran <- readRDS('./scdata_tpm.rds') %>% 
  subset(.,subset = group_cell == 'Oocyte')
data_tmp_tran <- scdata_tran@assays$RNA$counts %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'Genes') %>% 
  filter(Genes %in% (data_tmp %>% rownames())) %>% 
  column_to_rownames(var = 'Genes') %>% 
  filter(
    rowSums(.)>0
  )
break_heatmap <- scdata_tran@meta.data %>% 
  arrange(group_cell)  %>% 
  mutate(group_stage_cell = factor(
    group_stage_cell,levels = c(
      "Secondary_Granulosa", "Early antral_Granulosa",
      "Antral_Granulosa", "Preovulatory_Granulosa"
    )
  )) %>%
  pull(group_stage_cell) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
break_heatmap <- scdata_tran@meta.data %>% 
  arrange(group_cell)  %>% 
  pull(group_stage) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
annotation_col = data.frame(
  Group = scdata_tran$group_stage
)
rownames(annotation_col) = colnames(data_tmp_tran)
ann_colors = list(Group = c(
  'Secondary' = '#C34F73',
  'Early antral' = '#4B8600',
  'Antral' = '#4DBBD5',
  'Preovulatory' = '#3C5488'
))
p_tran <- pheatmap::pheatmap(
  data_tmp_tran %>% as.matrix(), 
  scale = 'row',
  display_numbers = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE, 
  gaps_col = break_heatmap,
  annotation_names_row = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_legend = FALSE,
  annotation_colors = ann_colors,
  legend = FALSE,
  na_col = 'white',
  border_color = 'white',
  fontsize = 20, 
  fontsize_row = 10,
  angle_col = 45,
  main = 'mRNA'
)
p_tran_ggplot <- cowplot::plot_grid(p_tran$gtable)
p_tran_ggplot
p_protein_ggplot + p_tran_ggplot

library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(tibble)

data_scale <- read.csv('./data_copy_num_sorted.csv',row.names = 1)
data_gene_protein <- data_scale[,c('ProteinGroups','Genes')] %>% 
  rename(Accessions = ProteinGroups)

data_p <- read.csv(
  './data_p_kruskal_test_Granulosa.csv',
  row.names = 1)
gene_select_p <- data_gene_protein %>% 
  dplyr::filter(Genes %in% data_p$Genes) %>% 
  pull(Accessions)
scdata <- readRDS('./scdata_use_filteroutliers.rds')
FindVariableFeatures(scdata)
scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata$group_cell %>% unique()
scdata <- subset(scdata,subset = group_cell == 'Granulosa')
tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rowwise() %>%
  summarise(
    num_positive = sum(c_across(everything()) > 0)
  ) %>% 
  mutate(
    protein_group = rownames(scdata)
  )
gene_select <- tmp$protein_group[tmp$num_positive>2]
rownames(scdata) %>% length()
scdata <- subset(scdata,features = gene_select)
scdata <- subset(scdata,features = gene_select_p)
rownames(scdata) %>% length()
scdata@meta.data %>% colnames()
scdata_sce <- as.SingleCellExperiment(scdata)
scdata_sce <- prepSCE(scdata_sce,
                      kid = "group_cell",
                      gid = "group_stage",#group——id，
                      sid = "sample",#sample_id, 
                      drop=T)
pb <- aggregateData(scdata_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()

boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 16
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
ggsave(width = 16,height = 8,dpi = 300,filename = './fig_all_GC.png')


data_plot_gene <- data_plot %>% 
  as.data.frame() %>% 
  dplyr::filter(cluster == 4)
gene_use <- acore.list[[4]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:50) %>% pull(NAME)
library(tibble)
data_tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column(var = 'Accessions') %>% 
  filter(Accessions %in% gene_use) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  dplyr::select(-c(Accessions)) %>% 
  column_to_rownames(var = 'Genes')
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  mutate(group_stage_cell = factor(
    group_stage_cell,levels = c(
      "Secondary_Granulosa", "Early antral_Granulosa",
      "Antral_Granulosa", "Preovulatory_Granulosa"
    )
  )) %>%
  pull(group_stage_cell) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  pull(group_stage) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
annotation_col = data.frame(
  Group = scdata$group_stage
)
rownames(annotation_col) = colnames(data_tmp)
ann_colors = list(Group = c(
  'Secondary' = '#C34F73',
  'Early antral' = '#4B8600',
  'Antral' = '#4DBBD5',
  'Preovulatory' = '#3C5488'
))
p_pheatmap <- pheatmap::pheatmap(
  data_tmp, 
  scale = 'row',
  display_numbers = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  gaps_col = break_heatmap,
  annotation_names_row = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  annotation_colors = ann_colors,
  na_col = 'white',
  border_color = 'white',
  fontsize = 20, fontsize_row = 10,
  angle_col = 45,
  main = 'cluster 5(Top 50 genes)'
)
p_pheatmap_ggplot <- cowplot::plot_grid(p_pheatmap$gtable)
p_pheatmap_ggplot
ggsave(width = 16,height = 10,bg = 'white',filename = './p_pheatmap_3_GC.png')


scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
gene_select <- acore.list[[4]] %>% 
  filter(MEM.SHIP>0.7) %>% 
  filter(!(grepl('Gm',NAME))) %>% 
  pull(NAME)
DotPlot(
  object = scdata,
  scale = TRUE,
  features = gene_select,
  group.by = 'group_stage'
) +
  coord_flip() +
  scale_color_gradient(low = '#DAE5F0FF',high = '#26456EFF') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.text = element_text(size = 20)
  )
ProNames <- acore.list[[4]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:12) %>% pull(NAME)
data_plot <- FetchData(
  object = scdata,layer = 'counts',
  clean = 'none',vars = c(ProNames,'group_stage_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_stage_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  as.data.frame() %>%
  mutate(gene = protein)

ggplot(data = data_plot,mapping = aes(x = group_stage_cell,y = expression,fill = group_stage_cell)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 4,scales = 'free_y',strip.position = 'left') +
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



library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
acore.list <- acore(dat,cl=cl,min.acore=0.2)
genes<-acore.list[[5]] %>% 
  rename(Accessions = NAME) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  pull(Genes)
gene<-str_to_title(genes)
gene
keys(org.Mm.eg.db, keytype="SYMBOL") %>% head()

gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
p_go <-
  dotplot(ego,
          showCategory = 15,
          font.size = 20)
ggsave("./GO_Cluster5.png", dpi=300, width=6, height = 6, units = "in")
barplot(ego, showCategory=20,title="EnrichmentGO")
data_go <- ego@result
write.csv(data_go,'./data_go_cluster5.csv')
p_tmp_1 <- cowplot::plot_grid(
  p_pheatmap_ggplot,
  p_go,
  rel_widths = c(1.2, 1),
  axis = 't',
  align = 'v'
)
p_res <- cowplot::plot_grid(
  p_mfuzz,
  p_tmp_1,ncol = 1,
  rel_heights = c(1, 1.5),
  align = 'v'
)
ggsave(
  plot = p_res,
  filename = './p_res_pro_o.png',
  width = 20,bg = 'white',
  height = 24
)


save.image('./protein_gc.RData')
res_mfuzz <- lapply(1:16, function(i){
  data_mfuzz <- acore.list[[i]] %>% 
    mutate(group = i)
}) %>% do.call(rbind,.)  %>% 
  rename(Accessions = NAME) %>% 
  left_join(data_gene_protein,by = 'Accessions') %>% 
  left_join(data_p %>% dplyr::select(-group),by = 'Genes')
write.csv(res_mfuzz,'./res_mfuzz_protein_gc.csv')
res_mfuzz$group %>% table()

acore.list <- acore(dat,cl=cl,min.acore=0.2)
dir.create('./res_gc')
for(i in 1:length(acore.list)){
  library(clusterProfiler)
  print(i)
  genes<-acore.list[[i]] %>% 
    rename(Accessions = NAME) %>% 
    left_join(data_gene_protein,by = 'Accessions') %>% 
    pull(Genes)
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  file_name <- paste(
    "./GO_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")
  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}

library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)

gene_use <- acore.list[[5]] %>%
  arrange(desc(MEM.SHIP)) %>%
  pull(NAME)
gene_select <- data_gene_protein$Genes[data_gene_protein$Accessions %in% gene_use]
data_TPM <- read.csv('./data_TPM.csv',row.names = 1)
scdata_tran <- CreateSeuratObject(
  data_TPM, 
  min.cells = 0, 
  min.features = 300)
scdata_tran@meta.data <- read.csv('./data_meta.csv',row.names = 1)
FindVariableFeatures(scdata_tran)
scdata_tran <- NormalizeData(
  scdata_tran, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata_tran$group_stage <- factor(
  scdata_tran$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata_tran <- subset(scdata_tran,subset = group_cell == 'Oocyte')
tmp <- scdata_tran@assays$RNA$counts %>% as.data.frame() %>% 
  rowwise() %>%
  summarise(
    num_positive = sum(c_across(everything()) > 0)
  ) %>% 
  mutate(
    protein_group = rownames(scdata_tran)
  )
rownames(scdata_tran) %>% length()
scdata_tran <- subset(scdata_tran,features = gene_select)
rownames(scdata_tran) %>% length()

scdata_tran@meta.data %>% colnames()
scdata_tran_sce <- as.SingleCellExperiment(scdata_tran)
scdata_tran_sce <- prepSCE(scdata_tran_sce,
                      kid = "group_cell",
                      gid = "group_stage",#group——id，
                      sid = "sample",#sample_id, 
                      drop=T)
pb <- aggregateData(scdata_tran_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()


boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 12
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
p_mfuzz

save.image('./protein_o_5_tran.RData')
saveRDS(cl,'./protein_o_5_tran_cl.rds')
saveRDS(dat,'./protein_o_5_tran_dat.rds')


acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp <- dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
colnames(tmp)
data_plot_tran <- tmp %>% 
  mutate(
    cluster = factor(
      cluster,levels = c(3,11,4,5,1,9,8,12,6,10,2,7)
    ),
    cluster_modify = case_when(
      cluster %in% c(3,11,4,5) ~ 'type II',
      cluster %in% c(1,9) ~ 'type III',
      cluster %in% c(8,12,6) ~ 'type IIII',
      cluster %in% c(10,2,7) ~ 'type I',
    )
  ) %>% 
  arrange(desc(cluster))
annotation_col = data.frame(
  Group = data_plot_tran$cluster_modify
)
rownames(annotation_col) = rownames(data_plot_tran)
ComplexHeatmap::pheatmap(
  mat = data_plot_tran %>% 
    select(-c('Gene_symbol','cluster','cluster_modify')),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = F,
  annotation_row = annotation_col,
  annotation_names_row = FALSE,
  show_rownames = FALSE,
  show_colnames = T,
  fontsize_col = 13,
  fontsize_row = 12,
  legend = TRUE,
  angle_col = '0',
  heatmap_legend_param = list(title = 'Relative Expression',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)


acore.list <- acore(dat,cl=cl,min.acore=0)
dir.create('./res_cluster5_tran')
for(i in 1:length(acore.list)){
  library(clusterProfiler)
  print(i)
  genes<-acore.list[[i]] %>% 
    rename(Accessions = NAME) %>% 
    pull(Accessions)
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  file_name <- paste(
    "./GO_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")

  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}


library(ComplexHeatmap)
gene_exp <- dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
colnames(tmp)
data_plot_tran <- tmp %>% 
  mutate(
    cluster = factor(
      cluster,levels = c(3,11,4,5,1,9,8,12,6,10,2,7)
    ),
    cluster_modify = case_when(
      cluster %in% c(3,11,4,5) ~ 'type I',
      cluster %in% c(1,9) ~ 'type II',
      cluster %in% c(8,12,6) ~ 'type III',
      cluster %in% c(10,2,7) ~ 'type IIII',
    ) %>% 
      factor(.,levels = c('type I','type II','type III','type IIII'))
  ) %>% 
  arrange(cluster_modify)
random_text = function(n) {
  sapply(1:n, function(i) {
    paste0(sample(letters, sample(4:10, 1)), collapse = "")
  })
}
text_list = list(
  text1 = random_text(4),
  text2 = random_text(4),
  text3 = random_text(4),
  text4 = random_text(4)
)

ha_top = HeatmapAnnotation(
  foo_top = anno_block(
    gp = gpar(fill = c('#C34F73','#0088D4','#4B8600','#5F559B')),
    labels = data_plot_pro %>% colnames() %>% .[2:5],
    labels_gp = gpar(col = "white", fontsize = 10)
  )
)
ha = rowAnnotation(
  foo = anno_empty(
    border = FALSE,
    width = max_text_width(unlist(text_list)) + unit(4, "mm"))
)
Heatmap(
  data_plot_tran %>% 
    select(-c('Gene_symbol','cluster','cluster_modify')) %>% 
    as.matrix(),
  name = "mat",
  row_split = data_plot_tran$cluster_modify %>% as.character(),
  top_annotation = ha_top,
  right_annotation = ha,
  column_split = 1:4,
  column_title = NULL,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(title = 'Relative Expression',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
color_use <- c('#C34F73','#0088D4','#4B8600','#5F559B')
for(i in 1:4) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = color_use[i], col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}







library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)


load('./protein_o.RData')

gene_use <- acore.list[[9]] %>%
  arrange(desc(MEM.SHIP)) %>%
  pull(NAME)
gene_select <- data_gene_protein$Genes[data_gene_protein$Accessions %in% gene_use]
data_TPM <- read.csv('./data_TPM.csv',row.names = 1)
scdata_tran <- CreateSeuratObject(
  data_TPM, 
  min.cells = 0, 
  min.features = 300)
scdata_tran@meta.data <- read.csv('./data_meta.csv',row.names = 1)
FindVariableFeatures(scdata_tran)
scdata_tran <- NormalizeData(
  scdata_tran, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata_tran$group_stage <- factor(
  scdata_tran$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata_tran <- subset(scdata_tran,subset = group_cell == 'Oocyte')
rownames(scdata_tran) %>% length()
scdata_tran <- subset(scdata_tran,features = gene_select)
rownames(scdata_tran) %>% length()

scdata_tran@meta.data %>% colnames()
scdata_tran_sce <- as.SingleCellExperiment(scdata_tran)
scdata_tran_sce <- prepSCE(scdata_tran_sce,
                           kid = "group_cell",
                           gid = "group_stage",#group——id，
                           sid = "sample",#sample_id, 
                           drop=T)
pb <- aggregateData(scdata_tran_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()


boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 12
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)
p_mfuzz

save.image('./protein_o_9_tran.RData')
saveRDS(cl,'./protein_o_9_tran_cl.rds')
saveRDS(dat,'./protein_o_9_tran_dat.rds')


acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp <- dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
colnames(tmp)
data_plot_tran <- tmp %>% 
  mutate(
    cluster = factor(
      cluster,levels = c(9,12,6,1,3,11,2,8,10,4,5,7)
    ),
    cluster_modify = case_when(
      cluster %in% c(9,12,6,1,3,11) ~ 'type II',
      cluster %in% c(2) ~ 'type III',
      cluster %in% c(8,10) ~ 'type IIII',
      cluster %in% c(4,5,7) ~ 'type I',
    )
  ) %>% 
  arrange(desc(cluster))
annotation_col = data.frame(
  Group = data_plot_tran$cluster
)
rownames(annotation_col) = rownames(data_plot_tran)
ComplexHeatmap::pheatmap(
  mat = data_plot_tran %>% 
    select(-c('Gene_symbol','cluster','cluster_modify')),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = F,
  annotation_row = annotation_col,
  annotation_names_row = FALSE,
  show_rownames = FALSE,
  show_colnames = T,
  fontsize_col = 13,
  fontsize_row = 12,
  legend = TRUE,
  angle_col = '0',
  heatmap_legend_param = list(title = 'Relative Expression',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)


acore.list <- acore(dat,cl=cl,min.acore=0)
unlink('./res_cluster9_tran')
dir.create('./res_cluster9_tran')
for(i in 1:length(acore.list)){
  library(clusterProfiler)
  print(i)
  genes<-acore.list[[i]] %>% 
    rename(Accessions = NAME) %>% 
    pull(Accessions)
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  file_name <- paste(
    "./GO_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")

  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}


library(ComplexHeatmap)
gene_exp <- dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
colnames(tmp)
data_plot_tran <- tmp %>% 
  mutate(
    cluster = factor(
      cluster,levels = c(9,12,6,1,3,11,2,8,10,4,5,7)
    ),
    cluster_modify = case_when(
      cluster %in% c(9,12,6,1,3,11) ~ 'type I',
      cluster %in% c(2) ~ 'type II',
      cluster %in% c(8,10) ~ 'type III',
      cluster %in% c(4,5,7) ~ 'type IIII',
    ) %>% 
      factor(.,levels = c('type I','type II','type III','type IIII'))
  ) %>% 
  arrange(cluster_modify)
random_text = function(n) {
  sapply(1:n, function(i) {
    paste0(sample(letters, sample(4:10, 1)), collapse = "")
  })
}
text_list = list(
  text1 = random_text(4),
  text2 = random_text(4),
  text3 = random_text(4),
  text4 = random_text(4)
)

ha_top = HeatmapAnnotation(
  foo_top = anno_block(
    gp = gpar(fill = c('#C34F73','#0088D4','#4B8600','#5F559B')),
    labels = data_plot_pro %>% colnames() %>% .[2:5],
    labels_gp = gpar(col = "white", fontsize = 10)
  )
)
ha = rowAnnotation(
  foo = anno_empty(
    border = FALSE,
    width = max_text_width(unlist(text_list)) + unit(4, "mm"))
)
Heatmap(
  data_plot_tran %>% 
    select(-c('Gene_symbol','cluster','cluster_modify')) %>% 
    as.matrix(),
  name = "mat",
  row_split = data_plot_tran$cluster_modify %>% as.character(),
  bottom_annotation = ha_top,
  right_annotation = ha,
  column_split = 1:4,
  column_title = NULL,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(title = 'Relative Expression',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
color_use <- c('#C34F73','#0088D4','#4B8600','#5F559B')
for(i in 1:4) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = color_use[i], col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}







library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@meta.data <- read.csv('./data_meta.csv',row.names = 1)
FindVariableFeatures(scdata)
scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata <- subset(scdata,subset = group_cell == 'Oocyte')
tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rowwise() %>%
  summarise(
    num_positive = sum(c_across(everything()) > 0)
  ) %>% 
  mutate(
    protein_group = rownames(scdata)
  )
gene_select <- tmp$protein_group[tmp$num_positive>2]
rownames(scdata) %>% length()
scdata <- subset(scdata,features = gene_select)
rownames(scdata) %>% length()

data_p <- read.csv(
  './data_p_kruskal_test_Oocyte.csv',
  row.names = 1) %>% 
  dplyr::select(c('Genes','pvalue_Permutation_tran','group')) %>% 
  dplyr::filter(pvalue_Permutation_tran<0.05)
gene_select_p <- data_p %>% 
  pull(Genes)
scdata <- subset(scdata,features = gene_select_p)
rownames(scdata) %>% length()

scdata@meta.data %>% colnames()
scdata_sce <- as.SingleCellExperiment(scdata)
scdata_sce <- prepSCE(scdata_sce,
                      kid = "group_cell",
                      gid = "group_stage",#group——id，
                      sid = "sample",#sample_id, 
                      drop=T)
pb <- aggregateData(scdata_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()


boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 16
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)#
ggsave(width = 16,height = 8,dpi = 300,filename = './fig_all.png')


data_plot_gene <- data_plot %>% 
  as.data.frame() %>% 
  dplyr::filter(cluster == 6)
gene_use <- acore.list[[6]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:50) %>% pull(NAME)
library(tibble)
data_tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rownames_to_column(var = 'genes') %>% 
  filter(genes %in% gene_use) %>% 
  column_to_rownames(var = 'genes')
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  mutate(group_stage_cell = factor(
    group_stage_cell,levels = c(
      "Secondary_Granulosa", "Early antral_Granulosa",
      "Antral_Granulosa", "Preovulatory_Granulosa"
    )
  )) %>%
  pull(group_stage_cell) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
break_heatmap <- scdata@meta.data %>% 
  arrange(group_cell)  %>% 
  pull(group_stage) %>% 
  table() %>% as.numeric()
break_heatmap <- lapply(1:length(break_heatmap), function(x){sum(break_heatmap[1:x])}) %>% unlist()
annotation_col = data.frame(
  Group = scdata$group_stage
)
rownames(annotation_col) = colnames(data_tmp)
ann_colors = list(Group = c(
  'Secondary' = '#C34F73',
  'Early antral' = '#4B8600',
  'Antral' = '#4DBBD5',
  'Preovulatory' = '#3C5488'
))
p_pheatmap <- pheatmap::pheatmap(
  data_tmp, 
  scale = 'row',
  display_numbers = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  gaps_col = break_heatmap,
  annotation_names_row = FALSE,
  annotation_col = annotation_col,
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  annotation_colors = ann_colors,
  na_col = 'white',
  border_color = 'white',
  fontsize = 20, fontsize_row = 10,
  angle_col = 45,
  main = 'cluster 6(Top 50 genes)'
)
p_pheatmap_ggplot <- cowplot::plot_grid(p_pheatmap$gtable)
p_pheatmap_ggplot

scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
gene_select <- acore.list[[6]] %>% 
  filter(MEM.SHIP>0.7) %>% 
  filter(!(grepl('Gm',NAME))) %>% 
  pull(NAME)
DotPlot(
  object = scdata,
  scale = TRUE,
  features = gene_select,
  group.by = 'group_stage'
) +
  coord_flip() +
  scale_color_gradient(low = '#DAE5F0FF',high = '#26456EFF') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.text = element_text(size = 20)
  )
ProNames <- acore.list[[6]] %>% 
  arrange(desc(MEM.SHIP)) %>% 
  slice(1:12) %>% pull(NAME)
data_plot <- FetchData(
  object = scdata,layer = 'counts',
  clean = 'none',vars = c(ProNames,'group_stage_cell')) %>% 
  rownames_to_column('sample') %>% 
  reshape2::melt(id.var = c('sample','group_stage_cell'),
                 variable.name = 'protein',
                 value.name = 'expression') %>% 
  as.data.frame() %>%
  mutate(gene = protein)

ggplot(data = data_plot,mapping = aes(x = group_stage_cell,y = expression,fill = group_stage_cell)) +
  geom_violin(scale = 'width',alpha = 0.8,trim = FALSE) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.7) +
  facet_wrap(~gene,ncol = 4,scales = 'free_y',strip.position = 'left') +
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



library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
acore.list <- acore(dat,cl=cl,min.acore=0.2)
genes<-acore.list[[6]][,1]
gene<-str_to_title(genes)
gene
keys(org.Mm.eg.db, keytype="SYMBOL") %>% head()

gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
p_go <-
  dotplot(ego,
          showCategory = 15,
          font.size = 20)
barplot(ego, showCategory=20,title="EnrichmentGO")
data_go <- ego@result
write.csv(data_go,'./data_go_cluster6.csv')

p_tmp_1 <- cowplot::plot_grid(
  p_pheatmap_ggplot,
  p_go,
  rel_widths = c(1.2, 1),
  axis = 't',
  align = 'v'
)
p_res <- cowplot::plot_grid(
  p_mfuzz,
  p_tmp_1,ncol = 1,
  rel_heights = c(1, 1.5),
  align = 'v'
)
ggsave(
  plot = p_res,
  filename = './p_res_pro_gc.png',
  width = 20,bg = 'white',
  height = 24
)

save.image('./trans_o.RData')
load('./trans_o.RData')
res_mfuzz <- lapply(1:16, function(i){
  data_mfuzz <- acore.list[[i]] %>% 
    mutate(group = i)
}) %>% do.call(rbind,.) %>%  
  rename(Genes = NAME,cluster = group) %>% 
  left_join(data_p,by = 'Genes') %>% 
  arrange(cluster,group)
write.csv(res_mfuzz,'./res_mfuzz_tran_o.csv')

acore.list <- acore(dat,cl=cl,min.acore=0.2)
dir.create('./res')
for(i in 1:length(acore.list)){
  print(i)
  genes<-acore.list[[i]] %>% 
    pull(NAME)
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  file_name <- paste(
    "./GO_Cluster_",
    i,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")

  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    i,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@meta.data <- read.csv('./data_meta.csv',row.names = 1)
FindVariableFeatures(scdata)
scdata <- NormalizeData(
  scdata, 
  assay = 'RNA',
  normalization.method = "LogNormalize", 
  scale.factor = 10000)

scdata$group_stage <- factor(
  scdata$group_stage,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )
)
scdata$group_cell %>% unique()
scdata <- subset(scdata,subset = group_cell == 'Granulosa')
tmp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  rowwise() %>%
  summarise(
    num_positive = sum(c_across(everything()) > 0)
  ) %>% 
  mutate(
    protein_group = rownames(scdata)
  )
gene_select <- tmp$protein_group[tmp$num_positive>2]
rownames(scdata) %>% length()
scdata <- subset(scdata,features = gene_select)
rownames(scdata) %>% length()

data_p <- read.csv(
  './data_p_kruskal_test_Granulosa.csv',
  row.names = 1) %>% 
  dplyr::select(c('Genes','pvalue_Permutation_tran','group')) %>% 
  dplyr::filter(pvalue_Permutation_tran<0.05)
gene_select_p <- data_p %>% 
  pull(Genes)
scdata <- subset(scdata,features = gene_select_p)
rownames(scdata) %>% length()

scdata@meta.data %>% colnames()
scdata_sce <- as.SingleCellExperiment(scdata)
scdata_sce <- prepSCE(scdata_sce,
                      kid = "group_cell",
                      gid = "group_stage",#group——id，
                      sid = "sample",#sample_id, 
                      drop=T)
pb <- aggregateData(scdata_sce,
                    assay = "counts", 
                    fun = "mean",
                    by = c("group_id"))
DEGs_exp_averp <- pb@assays@data@listData[[1]] %>% 
  as.data.frame() %>% 
  dplyr::select(Secondary,`Early antral`,Antral,Preovulatory) %>% 
  as.matrix()


boxplot(DEGs_exp_averp)
dat <- new(
  'ExpressionSet',
  exprs = DEGs_exp_averp)
dat <- filter.NA(dat, thres = 0.25)
dat <- fill.NA(dat, mode = 'mean')
dat <- filter.std(dat, min.std = 0)
dat <- standardise(dat)
c <- 16
m <- mestimate(dat)

set.seed(1234)
cl <- mfuzz(dat, c = c, m = m)
library(RColorBrewer)
Color <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



acore.list <- acore(dat,cl=cl,min.acore=0)
gene_exp<-dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
data_plot<-melt(tmp,id.vars=c('Gene_symbol','cluster'),variable.name = 'time',
                measure.vars=measure_vars,value.name = 'values')


plot_single_cluster=function(data_plot=data_plot,cluster=2){
  data_plot=data_plot[data_plot$cluster==cluster,]
  tmp <- acore.list[[cluster]]
  for(gene in tmp$NAME){
    data_plot$MEM.SHIP[data_plot$Gene_symbol == gene] <- tmp$MEM.SHIP[tmp$NAME==gene]
  }
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
          axis.text.x = element_text(angle = 45,hjust = 1,size = 20),
          plot.title = element_text(hjust = 0.5,size = 22),
          legend.position = 'none')+
    labs(x = 'Stage',y = 'Expression',title = paste('Cluster ',cluster,sep = ''))
  return(p)
}
plot_single_cluster(data_plot=data_plot,cluster = 1)
plot_mfuzz <- function(data_plot=data_plot,plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  features<-sort(unique(data_plot$cluster))
  plot_list <- purrr::map(features, function(x) plot_single_cluster(data_plot=data_plot,cluster=x))
  for(i in 1:8){
    plot_list[[i]]<- plot_list[[i]] +
      theme(axis.text.x = element_blank())+
      labs(x='')
  }
  for(i in seq(1,12,1)[-seq(1,12,4)]){
    plot_list[[i]]<- plot_list[[i]] +
      labs(y='')
  }
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)#,heights = 1
  return(p)
}
p_mfuzz <- plot_mfuzz(data_plot)#
ggsave(width = 16,height = 8,dpi = 300,filename = './fig_all.png')


save.image('./trans_gc.RData')
load('./trans_gc.RData')
res_mfuzz <- lapply(1:16, function(i){
  data_mfuzz <- acore.list[[i]] %>% 
    mutate(group = i)
}) %>% do.call(rbind,.) %>%  
  rename(Genes = NAME,cluster = group) %>% 
  left_join(data_p,by = 'Genes') %>% 
  arrange(cluster,group)
write.csv(res_mfuzz,'./res_mfuzz_tran_gc.csv')

library(Mfuzz)
library(Seurat)
library(muscat)
library(SingleCellExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(msigdbr)


load('./protein_o.RData')
protein_cluster <- data.frame(
  'Accessions' = cl[["cluster"]] %>% names(),
  'cluster' = cl[["cluster"]]
)

data_Maternal <- read.csv('./maternal_gene.txt',header = FALSE)
acore.list




data_maternal_use <- data_gene_protein %>% 
  filter(Genes %in% (data_Maternal$V1 %>% str_to_title())) %>% 
  left_join(
    protein_cluster,by = 'Accessions'
  )
data_maternal_use$cluster %>% table()
p_mfuzz
protein_use <- data_maternal_use %>% 
  filter(cluster %in% c(1,7,8)) %>% 
  pull(Accessions)
data_plot <- FetchData(object = scdata,vars = c(protein_use,'group_stage')) %>% 
  pivot_longer(cols = all_of(protein_use),names_to = 'Accessions',values_to = 'copy_number') %>% 
  left_join(
    data_maternal_use,by = 'Accessions',multiple = 'first'
  )
data_plot$group_stage %>% unique()
my_comparisons <- list(
  c('Secondary','Early antral'),
  c('Secondary','antral'),
  c('Secondary','Preovulatory'),
  c('Early antral','Antral'),
  c('Early antral','Preovulatory'),
  c('Antral','Preovulatory')
)
ggplot(data = data_plot,aes(x = group_stage,y = copy_number)) +
  geom_violin(
    trim = F,
    size = 0.2,
    show.legend = F,
    width = 0.6
  ) + 
  geom_boxplot(aes(fill = group_stage),width = 0.1,show.legend = FALSE) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))),
    comparisons = list(
      c('Secondary','Preovulatory')
    ),
    label.y = 6.7
  ) +
  labs(y = NULL, x = NULL) +
  scale_fill_manual(values = c('#C34F73','#0088D4','#4B8600','#5F559B')) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      colour = "black",
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.line = element_line(size = 0.2, color = "black"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length = unit(.5, "mm")
  )


gene_set <- msigdbr(species = 'Mus musculus')
colnames(gene_set)
gene_select <- gene_set %>% 
  dplyr::select(gs_name,gene_symbol) %>% 
  filter(grepl('OOCYTE_MEIOSIS',gs_name)) %>% 
  pull(gene_symbol) %>% unique()

data_meiosis_use <- data_gene_protein %>% 
  filter(Genes %in% gene_select) %>% 
  left_join(
    protein_cluster,by = 'Accessions'
  )
data_meiosis_use$Genes
data_meiosis_use$cluster %>% table()
p_mfuzz
protein_use <- data_meiosis_use %>% 
  filter(cluster %in% c(1,5)) %>% 
  pull(Accessions)


library(ComplexHeatmap)
gene_exp <- dat@assayData$exprs
gene_cluster<-as.data.frame(cl$cluster)
tmp<-as.data.frame(gene_exp)
tmp$Gene_symbol<-rownames(tmp)
tmp$cluster<-apply(tmp, 1, function(x){
  gene_cluster[rownames(gene_cluster)==x[length(x)],1]
})
measure_vars<-colnames(tmp)[c(-length(colnames(tmp)),-length(colnames(tmp))+1)]
colnames(tmp)
data_plot <- tmp %>% 
  mutate(
    cluster = factor(
      cluster,levels = c(9,12,6,1,3,11,2,8,10,4,5,7)
    ),
    cluster_modify = case_when(
      cluster %in% c(9,12,6,1,3,11) ~ 'type I',
      cluster %in% c(2) ~ 'type II',
      cluster %in% c(8,10) ~ 'type III',
      cluster %in% c(4,5,7) ~ 'type IIII',
    ) %>% 
      factor(.,levels = c('type I','type II','type III','type IIII'))
  ) %>% 
  arrange(cluster_modify) %>% 
  filter(
    Gene_symbol %in% data_meiosis_use$Accessions
  )
random_text = function(n) {
  sapply(1:n, function(i) {
    paste0(sample(letters, sample(4:10, 1)), collapse = "")
  })
}
text_list = list(
  text1 = random_text(4),
  text2 = random_text(4),
  text3 = random_text(4),
  text4 = random_text(4)
)

ha_top = HeatmapAnnotation(
  foo_top = anno_block(
    gp = gpar(fill = c('#C34F73','#0088D4','#4B8600','#5F559B')),
    labels = data_plot_pro %>% colnames() %>% .[2:5],
    labels_gp = gpar(col = "white", fontsize = 10)
  )
)
ha = rowAnnotation(
  foo = anno_empty(
    border = FALSE,
    width = max_text_width(unlist(text_list)) + unit(4, "mm"))
)
Heatmap(
  data_plot %>% 
    dplyr::select(-c('Gene_symbol','cluster','cluster_modify')) %>% 
    as.matrix(),
  name = "mat",
  row_split = 4,
  top_annotation = ha_top,
  right_annotation = ha,
  column_split = 1:4,
  column_title = NULL,
  cluster_rows = TRUE,
  cluster_columns = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_names_centered = TRUE,
  heatmap_legend_param = list(title = 'Relative Expression',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
color_use <- c('#C34F73','#0088D4','#4B8600','#5F559B')
for(i in 1:4) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = color_use[i], col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}


data_plot <- FetchData(object = scdata,vars = c(protein_use,'group_stage')) %>% 
  pivot_longer(cols = all_of(protein_use),names_to = 'Accessions',values_to = 'copy_number') %>% 
  left_join(
    data_meiosis_use,by = 'Accessions',multiple = 'first'
  )
data_plot$group_stage %>% unique()
my_comparisons <- list(
  c('Secondary','Early antral'),
  c('Secondary','antral'),
  c('Secondary','Preovulatory'),
  c('Early antral','Antral'),
  c('Early antral','Preovulatory'),
  c('Antral','Preovulatory')
)
ggplot(data = data_plot,aes(x = group_stage,y = copy_number)) +
  geom_violin(
    trim = F,
    size = 0.2,
    show.legend = F,
    width = 0.6
  ) + 
  geom_boxplot(aes(fill = group_stage),width = 0.1,show.legend = FALSE) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))),
    comparisons = list(
      c('Early antral','Antral'),
      c('Early antral','Preovulatory'),
      c('Secondary','Preovulatory')
    ),
  ) +
  labs(y = NULL, x = NULL) +
  scale_fill_manual(values = c('#C34F73','#0088D4','#4B8600','#5F559B')) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      colour = "black",
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.line = element_line(size = 0.2, color = "black"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length = unit(.5, "mm")
  )





my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bscdatad22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
res_mfuzz_pro <- read.csv('./res_mfuzz_protein_gc.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_gc.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","5","8","11","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('7'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("6",'11','13','15','16'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3",'5',"10"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_plot <- list(
  'pro_inc' = pro_inc %>% unique(),
  'pro_dec' = pro_dec %>% unique(),
  'tran_dec' = tran_dec %>% unique(),
  'tran_inc' = tran_inc %>% unique()
)
color_use <- c(
  '#C87E80','#AE9EC2','#F9BD90','#8DA8C3',#'#FDF0C4',
  rep('#FDF0C4',10)
)
p_venn <- plot(
  euler(data_plot, shape = "ellipse", order = TRUE),
  quantities = list(
    type = c("counts"),
    col = "black",
    font = 2,
    cex = 1
  ),
  fills = list(fill = color_use),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 2
  ),
  edges = list(col = color_use, lwd = 2, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
p_venn
pdf('./GC_venn.pdf',width = 8,height = 6)
p_venn
dev.off()

data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
)

matrix_data <- do.call(rbind, lapply(names(data_res), function(name) {
  data.frame('List_Name' = name, 'Element' = data_res[[name]])
}))
matrix_data$List_Name %>% table()
write.csv(matrix_data,'./data_venn_GC.csv')

data_plot <- matrix_data$List_Name %>% table() %>% 
  as.data.frame() %>% 
  rename_all(~c('Cluster','Counts')) %>% 
  mutate(Cluster = factor(Cluster,levels = names(data_res))) %>% 
  arrange(Cluster) %>% 
  mutate(
    lable = paste('Cluster ',1:n(),'\n(n=',Counts,')',sep = '')
  )
ggplot(data = data_plot,aes(x = Cluster,y = Counts,fill = Cluster)) +
  geom_bar(stat = 'identity',width = .7) +
  geom_text(aes(y = Counts+300,label = lable),vjust = 1) +
  scale_fill_manual(values = my36colors) +
  scale_y_continuous(expand = c(0,10)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = 'black'),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
  )
ggsave('./cluster_of_8type_gc.png',width = 6,height = 3,bg = 'white')
ggsave('./cluster_of_8type_gc.pdf',width = 6,height = 3,bg = 'white')



library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)

data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc_without_tran_inc' = setdiff(pro_inc, tran_inc),
  'pro_dec_without_tran_dec' = setdiff(pro_dec, tran_dec),
  'tran_inc_without_pro_inc' = setdiff(tran_inc,pro_inc),
  'tran_dec_without_pro_dec' = setdiff(tran_dec,pro_dec)
)
dir.create('./go_mfuzz_6_GC')
for(i in 1:length(data_res)){
  library(clusterProfiler)
  print(i)
  genes<-data_res[[i]]
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  cluster_name <- names(data_res)[i]
  file_name <- paste(
    "./GO_Cluster_",
    cluster_name,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")
  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    cluster_name,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}

library(dplyr)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
res_mfuzz_pro <- read.csv('./res_mfuzz_protein_o.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_o.csv',row.names = 1)

value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","6","7","9","10","11","13"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('4'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3","8","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("1","4","5","11","12","15","16"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc_without_tran_inc' = setdiff(pro_inc, tran_inc),
  'tran_inc_without_pro_inc' = setdiff(tran_inc,pro_inc)
)




library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
data(geneList)
for(class_use in names(data_res)){
  print(class_use)
  gene <- data_res[[class_use]] %>% 
    clusterProfiler::bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") %>% 
    dplyr::distinct(.,SYMBOL,.keep_all=TRUE)
  library(org.Mm.eg.db)
  ego <- enrichGO(
    gene = gene$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",

    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  tmp <- ego@result %>% 
    filter(
      p.adjust<0.05
    )
  ego <- pairwise_termsim(ego,showCategory = nrow(tmp))

  emapplot(ego, showCategory = nrow(tmp)) 
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_all.pdf',sep = ''
  )
  ggsave(
    filename = file_path,
    width = 50,height = 50,limitsize = FALSE
  )
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_all.png',sep = ''
  )
  ggsave(
    './ego_onlypro.png',
    width = 50,height = 50,limitsize = FALSE
  )
  emapplot(ego, showCategory = 300) 
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_300.pdf',sep = ''
  )
  ggsave(
    filename = file_path,
    width = 50,height = 50,limitsize = FALSE
  )
  emapplot(ego, showCategory = 200) 
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_200.pdf',sep = ''
  )
  ggsave(
    filename = file_path,
    width = 50,height = 50,limitsize = FALSE
  )
  emapplot(ego, showCategory = 100) 
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_100.pdf',sep = ''
  )
  ggsave(
    filename = file_path,
    width = 20,height = 20,limitsize = FALSE
  )
  emapplot(ego, showCategory = 50) 
  file_path <- paste(
    './figure4_venn/',
    class_use,
    'ego_onlypro_50.pdf',sep = ''
  )
  ggsave(
    filename = file_path,
    width = 50,height = 50,limitsize = FALSE
  )
}






ekegg <- enrichKEGG(
  gene = gene$ENTREZID,
  organism = 'mmu',
  pvalueCutoff = 0.05
)

ekegg <- pairwise_termsim(ekegg)
emapplot(ekegg, showCategory = 20)

library(enrichplot)
barplot(edo, showCategory=20)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)
heatplot(edox, foldChange=geneList)
emapplot(edo)

edo <- pairwise_termsim(edo)

emapplot(edo)

library(eulerr)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
res_mfuzz_pro <- read.csv('./res_mfuzz_protein_o.csv',row.names = 1)
colnames(res_mfuzz_pro)
res_mfuzz_pro <- res_mfuzz_pro %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
res_mfuzz_tran <- read.csv('./res_mfuzz_tran_o.csv',row.names = 1)
res_mfuzz_tran <- res_mfuzz_tran %>% 
  group_by(Genes) %>% 
  slice_max(order_by = MEM.SHIP, n = 1) %>%
  ungroup()
value_cuf <- 0
pro_inc <- res_mfuzz_pro %>% 
  filter(group %in% c("1","3","6","7","9","10","11","13"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
pro_dec <- res_mfuzz_pro %>% 
  filter(group %in% c('12'),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tmp <- pro_inc[pro_inc %in% pro_dec]
tmp
res_mfuzz_pro %>% 
  filter(group %in% c('4'),MEM.SHIP>value_cuf) %>% 
  filter(
    Genes %in% tmp
  ) %>% 
  arrange(desc(MEM.SHIP))
tran_inc <- res_mfuzz_tran %>% 
  filter(cluster %in% c("3","8","14"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
tran_dec <- res_mfuzz_tran %>% 
  filter(cluster %in% c("1","4","5","11","12","15","16"),MEM.SHIP>value_cuf) %>% 
  pull(Genes)
data_plot <- list(
  'pro_inc' = pro_inc %>% unique(),
  'pro_dec' = pro_dec %>% unique(),
  'tran_inc' = tran_inc %>% unique(),
  'tran_dec' = tran_dec %>% unique()
)
color_use <- c(
  '#C87E80','#AE9EC2','#F9BD90','#8DA8C3',#'#FDF0C4',
  rep('#FDF0C4',10)
)
p_venn <- plot(
  euler(data_plot, shape = "ellipse", order = TRUE),
  quantities = list(
    type = c("counts",'percent'),
    col = "black",
    font = 2,
    cex = 1
  ),
  fills = list(fill = color_use,aplha = 0.3),
  labels = list(
    labels = names(data_plot),
    col = "black",
    font = 3,
    cex = 2
  ),
  edges = list(col = color_use, lwd = 2, lty = 1),
  legend = list(
    labels = names(data_plot),
    font = 1,
    cex = 2,
    side = "right",
    x = 10
  )
)
pdf('./Oo_venn.pdf',width = 8,height = 6)
p_venn
dev.off()
data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
)



matrix_data <- do.call(rbind, lapply(names(data_res), function(name) {
  data.frame('List_Name' = name, 'Element' = data_res[[name]])
}))
matrix_data$List_Name %>% table()
write.csv(matrix_data,'./data_venn.csv')

data_plot <- matrix_data$List_Name %>% table() %>% 
  as.data.frame() %>% 
  rename_all(~c('Cluster','Counts')) %>% 
  mutate(Cluster = factor(Cluster,levels = names(data_res))) %>% 
  arrange(Cluster) %>% 
  mutate(
    lable = paste('Cluster ',1:n(),'\n(n=',Counts,')',sep = '')
  )
ggplot(data = data_plot,aes(x = Cluster,y = Counts,fill = Cluster)) +
  geom_bar(stat = 'identity',width = .7) +
  geom_text(aes(y = Counts+300,label = lable),vjust = 1) +
  scale_fill_manual(values = my36colors) +
  scale_y_continuous(expand = c(0,10)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = 'black'),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = 'none'
  )
ggsave('./cluster_of_8type.png',width = 6,height = 3,bg = 'white')
ggsave('./cluster_of_8type.pdf',width = 6,height = 3,bg = 'white')


library(org.Mm.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)

data_res <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc_without_tran_inc' = setdiff(pro_inc, tran_inc),
  'pro_dec_without_tran_dec' = setdiff(pro_dec, tran_dec),
  'tran_inc_without_pro_inc' = setdiff(tran_inc,pro_inc),
  'tran_dec_without_pro_dec' = setdiff(tran_dec,pro_dec)
)
dir.create('./go_mfuzz_6')
for(i in 1:length(data_res)){
  library(clusterProfiler)
  print(i)
  genes<-data_res[[i]]
  gene<-str_to_title(genes)

  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 

  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

  ego <- enrichGO(OrgDb="org.Mm.eg.db", gene = gene$ENTREZID, ont = "ALL", pvalueCutoff = 1, readable= TRUE)
  dotplot(ego,showCategory=30,title=paste("Enrichment GO Cluster",i),font.size = 20)
  cluster_name <- names(data_res)[i]
  file_name <- paste(
    "./GO_Cluster_",
    cluster_name,'.png',sep = ''
  )
  ggsave(file_name, dpi=300, width=12, height = 25, units = "in")
  data_go <- ego@result
  file_name <- paste(
    "./GO_Cluster_",
    cluster_name,'.csv',sep = ''
  )
  write.csv(data_go,file_name)
}



data_pathway_raw <- read.csv(
  './pathway_select.csv'
  ) %>% 
  tidyr::separate(col = GeneRatio,into = c('use','all'),sep = '/') %>% 
  mutate(
    generatio_use = as.numeric(use)/as.numeric(all),
    p.adjlog = -log10(p.adjust)
    ) %>% 
  dplyr::select(-c('use','all')) %>% 
  arrange(group,desc(generatio_use))
dir_path <- './figure4_pointplot_Oocyte/'
dir.create(dir_path,showWarnings = FALSE)
plot_point <- function(group_use,dir_path,width_use,height_use,bg_color){
  data_pathway_use <- data_pathway_raw %>% 
    filter(group == group_use) %>% 
    arrange(p.adjlog) %>% 
    mutate(
      Description = Description %>% str_wrap(width = 40),
      Description = factor(Description,levels= Description)
    )
  size_min <- data_pathway_use$generatio_use %>% min() %>% round(digits = 2)
  size_max <- data_pathway_use$generatio_use %>% max() %>%
    {.*100} %>% floor() %>% {./100}
  break_size <- c(
    size_min,size_min+(size_max-size_min)/2,size_max
  )
  size_min <- data_pathway_use$p.adjlog %>% min() %>% round(digits = 2)
  size_max <- data_pathway_use$p.adjlog %>% max() %>%
    {.*100} %>% floor() %>% {./100}
  break_p <- c(
    size_min-0.1,size_min+(size_max-size_min)/2,size_max-0.1
  )
  p <- ggplot(data = data_pathway_use,
              aes(
                x = p.adjlog,
                y = Description,
                size = generatio_use,
                color = p.adjlog
              )) +
    geom_point() +
    scale_size_continuous(range = c(3, 7),breaks = break_size) +
    guides(size  = guide_legend(title = 'Gene Ratio')) +
    scale_color_gradient(low = '#82A0D2',high = '#3975B9',guide = guide_colorbar(title = '-Log(p.adjust)'),breaks = break_p) +
    labs(x = '-Log(p.adjust)', y = '', title = 'GSEA of Oocyte') +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA,linewidth = 0.5),
      panel.background = element_rect(fill = 'white'),

      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      plot.title = element_blank(),
      axis.line =element_blank(),
      axis.text = element_text(size = 13,face = 'bold'),
      axis.text.x = element_text(angle = 0,hjust = 0.5,size = 13)
    )
  file_path <- paste(
    dir_path,
    group_use,'.png',sep = ''
  )
  ggsave(filename = file_path,plot = p,width = width_use,height = height_use,bg = bg_color)
  file_path <- paste(
    dir_path,
    group_use,'.pdf',sep = ''
  )
  ggsave(filename = file_path,plot = p,width = width_use,height = height_use,bg = bg_color)
  return(p)
}
data_pathway_raw$group %>% table()
plot_point(group_use = 'increase',dir_path = dir_path,width_use = 7.5,height_use = 5,bg_color = '#F3F3F1')
plot_point(group_use = 'transcript',dir_path = dir_path,width_use = 7.5,height_use = 5.5,bg_color = '#F5DEDB')
plot_point(group_use = 'protein',dir_path = dir_path,width_use = 7.5,height_use = 10.5,bg_color = '#D8EBF8')



library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Seurat)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
setwd('./')
replace_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  

  non_outliers <- x[x >= lower_bound & x <= upper_bound]
  tmp <- x[x <= lower_bound | x >= upper_bound]
  if (length(non_outliers) > 0) {

    x[x > upper_bound] <- max(non_outliers, na.rm = TRUE)
  }
  return(x)
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

scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
data_tran <- scdata_tran@assays$RNA$data_tpm %>% as.data.frame()
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



`%notin%` <- Negate(f = `%in%`)
tmp <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(group !='protein') %>% 
  pull(Description)
pathway_select <- c(
  "Golgi vesicle transport",
  "establishment of protein localization to organelle",
  "generation of precursor metabolites and energy",
  "proteasome-mediated ubiquitin-dependent protein catabolic process",
  "autophagy","actin filament organization",
  "response to oxidative stress",
  "response to peptide hormone"
)
data_pathway_raw <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(group =='protein', Description %in% pathway_select) %>% 
  arrange(desc(-log10(p.adjust)))
data_pathway_gene <- data_pathway_raw %>% 
  pull(geneID) %>% 
  str_split(pattern = '/')
names(data_pathway_gene) <- data_pathway_raw$Description
data_gene_use <- stack(data_pathway_gene) %>% 
  rename_all(~c('Genes','group')) %>% 
  mutate(group = factor(group,levels = names(data_pathway_gene)))
data_gene_use$group %>% table()



scdata_pro$group_stage_cell
scdata_pro_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte')
data_plot_raw <- lapply(data_pathway_gene, function(gene_use){
  data_pro_use <- FetchData(object = scdata_pro_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_pro_use %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    dplyr::filter(rowSums(.) > 0) %>%
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes')
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = data_pathway_raw$Description)
  ) %>% 
  arrange(pathway_name,Genes)
data_plot_raw$pathway_name %>% unique()
data_plot <- data_plot_raw %>% dplyr::select(-c(pathway_name,Genes))

min_value <- min(data_plot,na.rm = TRUE)
data_plot <- data_plot %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
gene_navalue <- data_plot %>% filter(is.na(rowSums(.))) %>% rownames()
data_plot[is.na(data_plot)] <- min(data_plot,na.rm = TRUE)
data_protein_use <- data_gene_use
library(ComplexHeatmap)
color_group_stage_cell <- c(
  '#F8B620','#FF800E','#ED444A','#D14B70'
)
group <- data.frame(
  group_stage_cell = colnames(data_plot) %>% str_remove(': Oocyte') %>% factor(.,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )),
  rowname = colnames(data_plot)
) %>% 
  column_to_rownames('rowname')
names(color_group_stage_cell) <- group$group_stage_cell
topanno = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  show_legend = FALSE,
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
draw(topanno)
group_row <- data.frame(
  group = data_protein_use$group %>% str_wrap(width = 30,whitespace_only = TRUE),
  rowname = rownames(data_plot)
) %>% 
  column_to_rownames('rowname')
color_group_pathway <- my36colors[1:length(group_row$group %>% unique())]
names(color_group_pathway) <- group_row$group %>% unique()
leftanno = rowAnnotation(
  df = group_row, # Row annotation
  border = FALSE,
  show_annotation_name = FALSE,show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 24),
  simple_anno_size = unit(20, "points"),
  gap = unit(12, "points"),
  col = list(
    group = color_group_pathway
  )
)
draw(leftanno)
ht <- Heatmap(
  matrix = data_plot %>% as.matrix(),
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
ht = draw(ht)

gene_order <- row_order(ht) %>% 
  stack()#ht@row_names_param[["labels"]]
gene_order <- rownames(data_plot) %>% 
  .[gene_order$values]
gene_order_name <- data_plot_raw %>%
  dplyr::select(c(Genes,pathway_name)) %>% 
  .[gene_order %>% as.numeric(),] %>% unstack()
data_plot_order <- data_plot %>% 
  rownames_to_column('Genes') %>% 
  mutate(
    Genes = factor(Genes,levels = gene_order)
  ) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
data_plot_order[rownames(data_plot_order) %in% gene_navalue,] <- NA
col_fun_pro <- colorRamp2(
  breaks = c(min(data_plot_order), 0, max(data_plot_order)),
  colors = c("#716CAC", "white", "#c2473b")
)
ht <- Heatmap(
  matrix = data_plot_order %>% as.matrix(),
  col = col_fun_pro,
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = 'Protein',
  column_title_gp = gpar(fontsize = 14, fontface = "bold", hjust = 1),
  row_split = data_protein_use$group,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
pdf('./heatmap_onlypro_increase.pdf',width = 13,height = 10,bg = 'white')
draw(
  ht,
  heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))
dev.off()

scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte')
data_tran_use <- FetchData(object = scdata_tran_sub,layer = 'counts',vars = c('group_stage_cell',data_gene_use$Genes))
data_plot_trans_raw <- lapply(gene_order_name, function(gene_use){
  data_tran_use <- FetchData(object = scdata_tran_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_tran_use %>% 
    mutate(group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte",'Preovulatory_Oocyte'
    ))) %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes') %>%  
    tidyr::complete(Genes = gene_use, fill = list(cluster = NA)) %>%
    mutate(Genes = factor(Genes,levels = gene_use)) %>% 
    arrange(Genes)
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = names(gene_order_name))
  )
data_plot_trans_raw$pathway_name %>% unique()
data_plot_tran <- data_plot_trans_raw %>% dplyr::select(-c(pathway_name,Genes))

min_value <- min(data_plot_tran,na.rm = TRUE)
data_plot_tran <- data_plot_tran %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
data_tran_use <- data_gene_use
topanno_tran = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(min(data_plot_tran), 0, max(data_plot_tran)),
  colors = c("#0D7291", "white", "#FC7136")
)
ht_tran <- Heatmap(
  matrix = data_plot_tran %>% as.matrix(),
  col = col_fun,
  top_annotation = topanno_tran,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_tran_use$group,
  row_title = NULL,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht_tran,heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))



ht_legend <- Legend(
  col_fun = col_fun_pro, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_legend)
ht_grob1 <- grid.grabExpr(draw(
  ht,
  heatmap_legend_side = NULL, 
))
ht_grob2 <- grid.grabExpr(draw(
  ht_tran,
  annotation_legend_list = list(ht_legend), 
  merge_legends = TRUE, 
  column_title = "Transcript",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", just = "center")
))


library(gridExtra)
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(2, 1))
pdf('./heatmap_proincrease.pdf',width = 13,height = 10,bg = 'white')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(2, 1))
dev.off()

data_cluster <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% 
  stack() %>% 
  rename_all(~c('Genes','Cluster')) 

tmp <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(group !='protein') %>% 
  pull(Description)
pathway_select <- c(
  "nuclear division"#,"cytoplasmic translation"
)
data_pathway_raw <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(Description %in% pathway_select) %>% 
  arrange(group,desc(-log10(p.adjust)))
data_pathway_gene <- data_pathway_raw %>% 
  pull(geneID) %>% 
  str_split(pattern = '/')
names(data_pathway_gene) <- data_pathway_raw$group
data_gene_use <- stack(data_pathway_gene) %>% 
  rename_all(~c('Genes','group')) %>% 
  mutate(group = factor(group,levels = names(data_pathway_gene))) %>% 
  left_join(
    data_cluster %>% 
      dplyr::select(Cluster,Genes),
    by = 'Genes') %>% 
  arrange(Cluster)
data_gene_use$group %>% table()
data_gene_use$Cluster %>% table()
data_pathway_gene <- data_gene_use[,c('Genes','Cluster')] %>% 
  unstack() %>% 
  {Filter(function(x) length(x) > 0, .)}

scdata_pro$group_stage_cell
scdata_pro_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte')
data_plot_raw <- lapply(data_pathway_gene, function(gene_use){
  data_pro_use <- FetchData(object = scdata_pro_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_pro_use %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    dplyr::filter(rowSums(.) > 0) %>%
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes') %>% 
    tidyr::complete(Genes = gene_use, fill = list(cluster = NA))
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = names(data_pathway_gene))
  ) %>% 
  arrange(pathway_name,Genes)
data_plot_raw$pathway_name %>% unique()
data_plot <- data_plot_raw %>% dplyr::select(-c(pathway_name)) %>% 
  column_to_rownames('Genes')

min_value <- min(data_plot,na.rm = TRUE)
data_plot <- data_plot %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
gene_navalue <- data_plot %>% filter(is.na(rowSums(.))) %>% rownames()
data_plot[is.na(data_plot)] <- min(data_plot,na.rm = TRUE)
data_protein_use <- data_gene_use
library(ComplexHeatmap)
color_group_stage_cell <- c(
  '#F8B620','#FF800E','#ED444A','#D14B70'
)
group <- data.frame(
  group_stage_cell = colnames(data_plot) %>% str_remove(': Oocyte') %>% factor(.,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )),
  rowname = colnames(data_plot)
) %>% 
  column_to_rownames('rowname')
names(color_group_stage_cell) <- group$group_stage_cell
topanno = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  show_legend = FALSE,
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
draw(topanno)
group_row <- data.frame(
  group = data_protein_use$Cluster %>% str_wrap(width = 30,whitespace_only = TRUE),
  rowname = rownames(data_plot)
) %>% 
  column_to_rownames('rowname')
color_group_pathway <- my36colors[1:length(group_row$group %>% unique())]
names(color_group_pathway) <- group_row$group %>% unique()
leftanno = rowAnnotation(
  df = group_row, # Row annotation
  border = FALSE,
  show_annotation_name = FALSE,show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 24),
  simple_anno_size = unit(20, "points"),
  gap = unit(12, "points"),
  col = list(
    group = color_group_pathway
  )
)
draw(leftanno)
ht <- Heatmap(
  matrix = data_plot %>% as.matrix(),
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_protein_use$Cluster,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
ht = draw(ht)

gene_order <- row_order(ht) %>% 
  stack()#ht@row_names_param[["labels"]]
gene_order <- rownames(data_plot) %>% 
  .[gene_order$values]
gene_order_name <- data_plot_raw %>%
  dplyr::select(c(Genes,pathway_name)) %>% 
  filter(Genes %in% gene_order) %>% 
  unstack()
data_plot_order <- data_plot %>% 
  rownames_to_column('Genes') %>% 
  mutate(
    Genes = factor(Genes,levels = gene_order)
  ) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
data_plot_order[rownames(data_plot_order) %in% gene_navalue,] <- NA
col_fun_pro <- colorRamp2(
  breaks = c(min(data_plot_order,na.rm = TRUE), 0, max(data_plot_order,na.rm = TRUE)),
  colors = c("#716CAC", "white", "#c2473b")
)
ht <- Heatmap(
  matrix = data_plot_order %>% as.matrix(),
  col = col_fun_pro,
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = 'Protein',
  column_title_gp = gpar(fontsize = 14, fontface = "bold", hjust = 1),
  row_split = data_protein_use$Cluster,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht,
  heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))

scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte')
data_tran_use <- FetchData(object = scdata_tran_sub,layer = 'counts',vars = c('group_stage_cell',data_gene_use$Genes))
data_plot_trans_raw <- lapply(gene_order_name, function(gene_use){
  data_tran_use <- FetchData(object = scdata_tran_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_tran_use %>% 
    mutate(group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte",'Preovulatory_Oocyte'
    ))) %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes') %>%  
    tidyr::complete(Genes = gene_use, fill = list(cluster = NA)) %>%
    mutate(Genes = factor(Genes,levels = gene_use)) %>% 
    arrange(Genes)
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = names(gene_order_name))
  )
data_plot_trans_raw$pathway_name %>% unique()
data_plot_tran <- data_plot_trans_raw %>% 
  dplyr::select(-c(pathway_name)) %>% 
  column_to_rownames('Genes')

min_value <- min(data_plot_tran,na.rm = TRUE)
data_plot_tran <- data_plot_tran %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
data_tran_use <- data_gene_use
topanno_tran = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(min(data_plot_tran), 0, max(data_plot_tran)),
  colors = c("#0D7291", "white", "#FC7136")
)
ht_tran <- Heatmap(
  matrix = data_plot_tran %>% as.matrix(),
  col = col_fun,
  top_annotation = topanno_tran,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_tran_use$Cluster,
  row_title = NULL,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht_tran,heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))



ht_legend <- Legend(
  col_fun = col_fun_pro, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_legend)
ht_grob1 <- grid.grabExpr(draw(
  ht,
  heatmap_legend_side = NULL, 
))
ht_grob2 <- grid.grabExpr(draw(
  ht_tran,
  annotation_legend_list = list(ht_legend), 
  merge_legends = TRUE, 
  column_title = "Transcript",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", just = "center")
))


library(gridExtra)
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
pdf('./heatmap_Nuclear-division.pdf',width = 13,height = 10,bg = 'white')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()

data_cluster <- list(
  'pro_inc&tran_inc' = intersect(pro_inc, tran_inc),
  'pro_inc&tran_des' = intersect(pro_inc, tran_dec),
  'pro_dec&tran_ins' = intersect(pro_dec, tran_inc),
  'pro_dec&tran_des' = intersect(pro_dec, tran_dec),
  'pro_inc' = setdiff(pro_inc, c(tran_inc,tran_dec)),
  'pro_dec' = setdiff(pro_dec, c(tran_inc, tran_dec)),
  'tran_inc' = setdiff(tran_inc,c(pro_inc,pro_dec)),
  'tran_dec' = setdiff(tran_dec,c(pro_inc,pro_dec))
) %>% 
  stack() %>% 
  rename_all(~c('Genes','Cluster')) 

tmp <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(group !='protein') %>% 
  pull(Description)
pathway_select <- c(
  "cytoplasmic translation"
)
data_pathway_raw <- read.csv(
  './pathway_select.csv'
) %>% 
  filter(Description %in% pathway_select) %>% 
  arrange(group,desc(-log10(p.adjust)))
data_pathway_gene <- data_pathway_raw %>% 
  pull(geneID) %>% 
  str_split(pattern = '/')
names(data_pathway_gene) <- data_pathway_raw$group
data_gene_use <- stack(data_pathway_gene) %>% 
  rename_all(~c('Genes','group')) %>% 
  mutate(group = factor(group,levels = names(data_pathway_gene))) %>% 
  left_join(
    data_cluster %>% 
      dplyr::select(Cluster,Genes),
    by = 'Genes') %>% 
  arrange(Cluster)
data_gene_use$group %>% table()
data_gene_use$Cluster %>% table()
data_pathway_gene <- data_gene_use[,c('Genes','Cluster')] %>% 
  unstack() %>% 
  {Filter(function(x) length(x) > 0, .)}

scdata_pro$group_stage_cell
scdata_pro_sub <- subset(scdata_pro,subset = group_cell == 'Oocyte')
data_plot_raw <- lapply(data_pathway_gene, function(gene_use){
  data_pro_use <- FetchData(object = scdata_pro_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_pro_use %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    dplyr::filter(rowSums(.) > 0) %>%
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes') %>% 
    tidyr::complete(Genes = gene_use, fill = list(cluster = NA))
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = names(data_pathway_gene))
  ) %>% 
  arrange(pathway_name,Genes)
data_plot_raw$pathway_name %>% unique()
data_plot <- data_plot_raw %>% dplyr::select(-c(pathway_name)) %>% 
  column_to_rownames('Genes')

min_value <- min(data_plot,na.rm = TRUE)
data_plot <- data_plot %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
gene_navalue <- data_plot %>% filter(is.na(rowSums(.))) %>% rownames()
data_plot[is.na(data_plot)] <- min(data_plot,na.rm = TRUE)
data_protein_use <- data_gene_use
library(ComplexHeatmap)
color_group_stage_cell <- c(
  '#F8B620','#FF800E','#ED444A','#D14B70'
)
group <- data.frame(
  group_stage_cell = colnames(data_plot) %>% str_remove(': Oocyte') %>% factor(.,levels = c(
    "Secondary","Early antral","Antral","Preovulatory"
  )),
  rowname = colnames(data_plot)
) %>% 
  column_to_rownames('rowname')
names(color_group_stage_cell) <- group$group_stage_cell
topanno = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  show_legend = FALSE,
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
draw(topanno)
group_row <- data.frame(
  group = data_protein_use$Cluster %>% str_wrap(width = 30,whitespace_only = TRUE),
  rowname = rownames(data_plot)
) %>% 
  column_to_rownames('rowname')
color_group_pathway <- my36colors[1:length(group_row$group %>% unique())]
names(color_group_pathway) <- group_row$group %>% unique()
leftanno = rowAnnotation(
  df = group_row, # Row annotation
  border = FALSE,
  show_annotation_name = FALSE,show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = 24),
  simple_anno_size = unit(20, "points"),
  gap = unit(12, "points"),
  col = list(
    group = color_group_pathway
  )
)
draw(leftanno)
ht <- Heatmap(
  matrix = data_plot %>% as.matrix(),
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_protein_use$Cluster,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
ht = draw(ht)

gene_order <- row_order(ht) %>% 
  stack()#ht@row_names_param[["labels"]]
gene_order <- rownames(data_plot) %>% 
  .[gene_order$values]
gene_order_name <- data_plot_raw %>%
  dplyr::select(c(Genes,pathway_name)) %>% 
  filter(Genes %in% gene_order) %>% 
  unstack()
data_plot_order <- data_plot %>% 
  rownames_to_column('Genes') %>% 
  mutate(
    Genes = factor(Genes,levels = gene_order)
  ) %>% 
  arrange(Genes) %>% 
  column_to_rownames('Genes')
data_plot_order[rownames(data_plot_order) %in% gene_navalue,] <- NA
col_fun_pro <- colorRamp2(
  breaks = c(min(data_plot_order,na.rm = TRUE), 0, max(data_plot_order,na.rm = TRUE)),
  colors = c("#716CAC", "white", "#c2473b")
)
ht <- Heatmap(
  matrix = data_plot_order %>% as.matrix(),
  col = col_fun_pro,
  top_annotation = topanno,
  left_annotation = leftanno,
  column_split = group$group_stage_cell,
  column_title = 'Protein',
  column_title_gp = gpar(fontsize = 14, fontface = "bold", hjust = 1),
  row_split = data_protein_use$Cluster,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht,
  heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))

scdata_tran_sub <- subset(scdata_tran,subset = group_cell == 'Oocyte')
data_tran_use <- FetchData(object = scdata_tran_sub,layer = 'counts',vars = c('group_stage_cell',data_gene_use$Genes))
data_plot_trans_raw <- lapply(gene_order_name, function(gene_use){
  data_tran_use <- FetchData(object = scdata_tran_sub,vars = c('group_stage_cell',gene_use))
  data_plot <- data_tran_use %>% 
    mutate(group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary_Oocyte","Early antral_Oocyte","Antral_Oocyte",'Preovulatory_Oocyte'
    ))) %>% 
    group_by(group_stage_cell) %>% 
    summarise(across(everything(), mean, na.rm = TRUE)) %>% 
    column_to_rownames('group_stage_cell') %>% 
    {log10(.+1)} %>% {scale(.)} %>%
    t() %>% as.data.frame() %>%
    rownames_to_column('Genes') %>%  
    tidyr::complete(Genes = gene_use, fill = list(cluster = NA)) %>%
    mutate(Genes = factor(Genes,levels = gene_use)) %>% 
    arrange(Genes)
}) %>% 
  do.call(rbind,.) %>% 
  rownames_to_column('pathway_name') %>% 
  mutate(
    pathway_name = pathway_name %>% str_remove(pattern = "\\.\\d+$"),
    pathway_name = factor(pathway_name,levels = names(gene_order_name))
  )
data_plot_trans_raw$pathway_name %>% unique()
data_plot_tran <- data_plot_trans_raw %>% 
  dplyr::select(-c(pathway_name)) %>% 
  column_to_rownames('Genes')

min_value <- min(data_plot_tran,na.rm = TRUE)
data_plot_tran <- data_plot_tran %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), min_value, .)))
data_tran_use <- data_gene_use
topanno_tran = HeatmapAnnotation(
  df = group, # Column annotation
  border = FALSE,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 12),
  simple_anno_size = unit(10, "points"),
  gap = unit(10, "points"),
  col = list(
    group_stage_cell = color_group_stage_cell
  ),
  annotation_legend_param = list(
    group_stage_cell = list(title = 'Stage',title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 12),grid_width = unit(7, "mm"),gap = unit(5, "mm"))
  )
)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(min(data_plot_tran), 0, max(data_plot_tran)),
  colors = c("#0D7291", "white", "#FC7136")
)
ht_tran <- Heatmap(
  matrix = data_plot_tran %>% as.matrix(),
  col = col_fun,
  top_annotation = topanno_tran,
  column_split = group$group_stage_cell,
  column_title = NULL,
  row_split = data_tran_use$Cluster,
  row_title = NULL,
  row_title_rot = 0,
  row_gap = unit(3, "mm"),
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(
    title = 'Relative Expression',
    title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
    title_position = 'leftcenter-rot',
    legend_direction = 'vertical',
    legend_height = unit(40, units = "mm"),
    labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
  )
)
draw(
  ht_tran,heatmap_legend_side = "right", # Place heatmap legend at the bottom
  annotation_legend_side = "right",
  merge_legends = TRUE,padding = unit(c(10, 10, 10, 30), "mm"))



ht_legend <- Legend(
  col_fun = col_fun_pro, 
  title = 'Relative Expression',
  title_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5),
  title_position = 'leftcenter-rot',
  legend_height = unit(40, units = "mm"),
  labels_gp = grid::gpar(fontsize = 12, fontface = "bold", hjust = 0.5, vjust = 0.5)
)
draw(ht_legend)
ht_grob1 <- grid.grabExpr(draw(
  ht,
  heatmap_legend_side = NULL, 
))
ht_grob2 <- grid.grabExpr(draw(
  ht_tran,
  annotation_legend_list = list(ht_legend), 
  merge_legends = TRUE, 
  column_title = "Transcript",
  column_title_gp = gpar(fontsize = 14, fontface = "bold", just = "center")
))


library(gridExtra)
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
pdf('./heatmap_Cytoplasmic-translation.pdf',width = 13,height = 10,bg = 'white')
grid.arrange(ht_grob1, ht_grob2, ncol = 2, widths = c(1.1, 1))
dev.off()
