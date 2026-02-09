dir.create('./result_0813')
setwd('./result_0813')


library(Seurat)
library(rliger)
library(SingleCellExperiment)
library(SCENIC)
library(RcisTarget)
scdata.sce <- readRDS('./scdata_tpm.rds')
scdata.sce@assays$RNA$data_fpkm <- NULL
scdata.sce@assays$RNA$data_tpm <- NULL
scdata.sce[['RNA']] <- as(object = scdata.sce[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(scdata.sce)
saveRDS(sce,'./sce.rds')
exprMat <- counts(sce) %>% as.matrix()
cellInfo <- colData(sce)


data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
save.image('./SCENIC_step2.RData')
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
export2loom(scenicOptions, exprMat)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)



minGenes = 20
coexMethods = NULL
minJakkardInd = 0.8
signifGenesMethod = "aprox"
onlyPositiveCorr = TRUE
onlyBestGsPerMotif = TRUE
dbIndexCol = "features"
nCores <- getSettings(scenicOptions, "nCores")
tfModules_asDF <- tryCatch(
  loadInt(scenicOptions, "tfModules_asDF"),
  error = function(e) {
    if (getStatus(scenicOptions, asID = TRUE) < 2)
      e$message <-
        paste0(
          "It seems the co-expression modules have not been built yet. Please, run runSCENIC_1_coexNetwork2modules() first.\n",
          e$message
        )
    stop(e)
  }
)
if (!is.null(coexMethods))
  tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$method %in%
                                           coexMethods),]
if (!is.null(minJakkardInd))
  tfModules_asDF <- mergeOverlappingModules(tfModules_asDF,
                                            minJakkardInd = minJakkardInd)
if (nrow(tfModules_asDF) == 0)
  stop("The co-expression modules are empty.")
if ("BiocParallel" %in% installed.packages() && (nCores >
                                                 1)) {
  library(BiocParallel)
  register(MulticoreParam(nCores), default = TRUE)
}
msg <-
  paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
if (getSettings(scenicOptions, "verbose"))
  message(msg)
library(AUCell)
library(RcisTarget)
motifAnnot <- getDbAnnotations(scenicOptions)
if (is.null(names(getSettings(scenicOptions, "dbs")))) {
  names(scenicOptions@settings$dbs) <- scenicOptions@settings$dbs
  tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"),
                         "-", fixed = T), function(x)
                           x[grep("bp|kb", x)])
  if (all(lengths(tmp) > 0))
    names(scenicOptions@settings$dbs) <- tmp
}
loadAttempt <-
  sapply(getDatabases(scenicOptions), dbLoadingAttempt,
         indexCol = dbIndexCol)
if (any(!loadAttempt))
  stop("It is not possible to load the following databses: \n",
       paste(dbs[which(!loadAttempt)], collapse = "\n"))
genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions),
                                  function(dbFilePath) {
                                    rf <- arrow::ReadableFile$create(dbFilePath)
                                    fr <-
                                      arrow::FeatherReader$create(rf)
                                    genesInDb <- names(fr)
                                    rnktype <- "features"
                                    genesInDb <-
                                      genesInDb[genesInDb != rnktype]
                                  })))
featuresWithAnnot <- checkAnnots(scenicOptions, motifAnnot)
if (any(featuresWithAnnot == 0))
  warning("Missing annotations\n", names(which(rankingsInDb ==
                                                 0)))
tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
allTFs <- getDbTfs(scenicOptions)
tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in%
                                         allTFs),]
geneInDb <- tfModules_asDF$Target %in% genesInDb
missingGene <- sort(unique(tfModules_asDF[which(!geneInDb),
                                          "Target"]))
if (length(missingGene) > 0)
  warning(paste0(
    "Genes in co-expression modules not available in RcisTargetDatabases: ",
    paste(missingGene, collapse = ", ")
  ))
tfModules_asDF <- tfModules_asDF[which(geneInDb),]
if (all(is.na(tfModules_asDF$corr))) {
  warning("no correlation info available")
  tfModules_Selected <- tfModules_asDF
  tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF,
                                          tfModules_Selected$method, sep = "_")
}else {
  tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr ==
                                               1),]
  tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF,
                                          tfModules_Selected$method, sep = "_")
  if (!onlyPositiveCorr) {
    tfModules_IgnCorr <- tfModules_asDF[which(tfModules_asDF$corr !=
                                                1),]
    tfModules_IgnCorr$geneSetName <- paste0(tfModules_IgnCorr$TF,
                                            "_", tfModules_IgnCorr$method)
    posCorr <-
      tfModules_Selected[which(tfModules_Selected$geneSetName %in%
                                 unique(tfModules_IgnCorr$geneSetName)),]
    tfModules_IgnCorr <- rbind(tfModules_IgnCorr, posCorr)
    tfModules_IgnCorr$geneSetName <-
      paste0(tfModules_IgnCorr$geneSetName,
             "IgnCorr")
    tfModules_Selected <- rbind(tfModules_Selected,
                                tfModules_IgnCorr)
  }
}
tfModules_Selected$geneSetName <-
  factor(as.character(tfModules_Selected$geneSetName))
allGenes <- unique(tfModules_Selected$Target)
tfModules <-
  split(tfModules_Selected$Target, tfModules_Selected$geneSetName)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))
tfModules <- tfModules[which(lengths(tfModules) >= minGenes)]
saveRDS(tfModules, file = getIntName(scenicOptions, "tfModules_forEnrichment"))
if (getSettings(scenicOptions, "verbose")) {
  tfModulesSummary <- t(sapply(strsplit(names(tfModules),
                                        "_"), function(x)
                                          x[1:2]))
  message("tfModulesSummary:")
  print(cbind(sort(table(
    tfModulesSummary[, 2]
  ))))
}
msg <-
  paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
if (getSettings(scenicOptions, "verbose"))
  message(msg)
motifs_AUC <-
  lapply(getDatabases(scenicOptions), function(rnkName) {
    ranking <- importRankings(rnkName, columns = allGenes)
    message("Scoring database: ", ranking@description)
    RcisTarget::calcAUC(
      tfModules,
      ranking,
      aucMaxRank = 0.03 *
        getNumColsInDB(ranking),
      nCores = nCores,
      verbose = FALSE
    )
  })
saveRDS(motifs_AUC, file = getIntName(scenicOptions, "motifs_AUC"))
msg <-
  paste0(format(Sys.time(), "%H:%M"),
         "\tRcisTarget: Adding motif annotation")
message(msg)
scenicOptions@settings$nCores <- 1
motifEnrichment <- lapply(motifs_AUC, function(aucOutput) {
  tf <- sapply(setNames(strsplit(rownames(aucOutput),
                                 "_"), rownames(aucOutput)), function(x)
                                   x[[1]])
  addMotifAnnotation(
    aucOutput,
    nesThreshold = 3,
    digits = 3,
    motifAnnot = motifAnnot,
    motifAnnot_highConfCat = c("directAnnotation",
                               "inferredBy_Orthology"),
    motifAnnot_lowConfCat = c(
      "inferredBy_MotifSimilarity",
      "inferredBy_MotifSimilarity_n_Orthology"
    ),
    highlightTFs = tf
  )
})
motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment),
                                         function(dbName) {
                                           cbind(motifDb = dbName, motifEnrichment[[dbName]])
                                         }))
saveRDS(motifEnrichment,
        file = getIntName(scenicOptions,
                          "motifEnrichment_full"))
msg <- paste0("Number of motifs in the initial enrichment: ",
              nrow(motifEnrichment))
if (getSettings(scenicOptions, "verbose"))
  message(msg)
motifEnrichment_selfMotifs <-
  motifEnrichment[which(motifEnrichment$TFinDB !=
                          ""), , drop = FALSE]
msg <- paste0("Number of motifs annotated to the matching TF: ",
              nrow(motifEnrichment_selfMotifs))
if (getSettings(scenicOptions, "verbose"))
  message(msg)
rm(motifEnrichment)
if (nrow(motifEnrichment_selfMotifs) == 0)
  stop(
    "None of the co-expression modules present enrichment of the TF motif: There are no regulons."
  )
if (onlyBestGsPerMotif) {
  met_byDb <-
    split(motifEnrichment_selfMotifs,
          motifEnrichment_selfMotifs$motifDb)
  for (db in names(met_byDb)) {
    met <- met_byDb[[db]]
    met <- split(met, factor(met$highlightedTFs))
    met <- lapply(met, function(x) {
      rbindlist(lapply(split(x, x$motif), function(y)
        y[which.max(y$NES),]))
    })
    met_byDb[[db]] <- rbindlist(met)
  }
  motifEnrichment_selfMotifs <- rbindlist(met_byDb)
  rm(met_byDb)
  rm(met)
}
msg <-
  paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Pruning targets")
if (getSettings(scenicOptions, "verbose"))
  message(msg)
dbNames <- getDatabases(scenicOptions)
motifEnrichment_selfMotifs_wGenes <- lapply(names(dbNames),
                                            function(motifDbName) {
                                              ranking <- importRankings(dbNames[motifDbName],
                                                                        columns = allGenes)
                                              addSignificantGenes(
                                                resultsTable = motifEnrichment_selfMotifs[motifEnrichment_selfMotifs$motifDb ==
                                                                                            motifDbName,],
                                                geneSets = tfModules,
                                                rankings = ranking,
                                                plotCurve = FALSE,
                                                maxRank = 5000,
                                                method = signifGenesMethod,
                                                nMean = 100,
                                                nCores = nCores
                                              )
                                            })
suppressPackageStartupMessages(library(data.table))
motifEnrichment_selfMotifs_wGenes <-
  rbindlist(motifEnrichment_selfMotifs_wGenes)
saveRDS(
  motifEnrichment_selfMotifs_wGenes,
  file = getIntName(scenicOptions,
                    "motifEnrichment_selfMotifs_wGenes")
)
if (getSettings(scenicOptions, "verbose")) {
  message(
    format(Sys.time(), "%H:%M"),
    "\tNumber of motifs that support the regulons: ",
    nrow(motifEnrichment_selfMotifs_wGenes)
  )
  motifEnrichment_selfMotifs_wGenes[order(motifEnrichment_selfMotifs_wGenes$NES,
                                          decreasing = TRUE),][1:5, (1:ncol(motifEnrichment_selfMotifs_wGenes) -
                                                                       1), with = F]
}
if (!file.exists("output"))
  dir.create("output")
write.table(
  motifEnrichment_selfMotifs_wGenes,
  file = getOutName(scenicOptions,
                    "s2_motifEnrichment"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
if ("DT" %in% installed.packages() &&
    nrow(motifEnrichment_selfMotifs_wGenes) >
    0) {
  nvm <- tryCatch({
    colsToShow <- c("motifDb",
                    "logo",
                    "NES",
                    "geneSet",
                    "TF_highConf",
                    "TF_lowConf")
    motifEnrichment_2html <-
      viewMotifs(
        motifEnrichment_selfMotifs_wGenes,
        colsToShow = colsToShow,
        options = list(pageLength = 100)
      )
    fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
    dirName <- dirname(fileName)
    fileName <- basename(fileName)
    suppressWarnings(DT::saveWidget(motifEnrichment_2html,
                                    fileName))
    file.rename(fileName, file.path(dirName, fileName))
    if (getSettings(scenicOptions, "verbose"))
      message("\tPreview of motif enrichment saved as: ",
              file.path(dirName, fileName))
  }, error = function(e)
    print(e$message))
}
motifEnrichment.asIncidList <-
  apply(motifEnrichment_selfMotifs_wGenes,
        1, function(oneMotifRow) {
          genes <- strsplit(oneMotifRow["enrichedGenes"],
                            ";")[[1]]
          oneMotifRow <-
            data.frame(rbind(oneMotifRow), stringsAsFactors = FALSE)
          data.frame(oneMotifRow[rep(1, length(genes)), c("NES",
                                                          "motif",
                                                          "highlightedTFs",
                                                          "TFinDB",
                                                          "geneSet",
                                                          "motifDb")], genes, stringsAsFactors = FALSE)
        })
motifEnrichment.asIncidList <-
  rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) ==
                                              "highlightedTFs")] <-
  "TF"
colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) ==
                                              "TFinDB")] <- "annot"
colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList) ==
                                              "genes")] <- "gene"
motifEnrichment.asIncidList <-
  data.frame(motifEnrichment.asIncidList,
             stringsAsFactors = FALSE)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList,
                                   motifEnrichment.asIncidList$TF), function(tfTargets) {
                                     tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets,
                                                                                          tfTargets$gene), function(enrOneGene) {
                                                                                            highConfAnnot <- "**" %in% enrOneGene$annot
                                                                                            enrOneGeneByAnnot <-
                                                                                              enrOneGene
                                                                                            if (highConfAnnot)
                                                                                              enrOneGeneByAnnot <-
                                                                                              enrOneGeneByAnnot[which(enrOneGene$annot ==
                                                                                                                        "**"),]
                                                                                            bestMotif <-
                                                                                              which.max(enrOneGeneByAnnot$NES)
                                                                                            tf <-
                                                                                              unique(enrOneGene$TF)
                                                                                            cbind(
                                                                                              TF = tf,
                                                                                              gene = unique(enrOneGene$gene),
                                                                                              highConfAnnot = highConfAnnot,
                                                                                              nMotifs = nrow(enrOneGene),
                                                                                              bestMotif = as.character(enrOneGeneByAnnot[bestMotif,
                                                                                                                                         "motif"]),
                                                                                              NES = as.numeric(enrOneGeneByAnnot[bestMotif,
                                                                                                                                 "NES"]),
                                                                                              motifDb = as.character(enrOneGeneByAnnot[bestMotif,
                                                                                                                                       "motifDb"]),
                                                                                              coexModule = gsub(
                                                                                                paste0(tf,
                                                                                                       "_"),
                                                                                                "",
                                                                                                as.character(enrOneGeneByAnnot[bestMotif,
                                                                                                                               "geneSet"]),
                                                                                                fixed = TRUE
                                                                                              )
                                                                                            )
                                                                                          })), stringsAsFactors = FALSE)
                                     tfTable[order(tfTable$NES, decreasing = TRUE),]
                                   })
rm(motifEnrichment.asIncidList)
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
corrMat <- loadInt(scenicOptions, "corrMat", ifNotExists = "null")
if (!is.null(corrMat)) {
  regulonTargetsInfo$spearCor <- NA_real_
  for (tf in unique(regulonTargetsInfo$TF)) {
    regulonTargetsInfo[which(regulonTargetsInfo$TF ==
                               tf), "spearCor"] <-
      corrMat[tf, unlist(regulonTargetsInfo[which(regulonTargetsInfo$TF ==
                                                    tf), "gene"])]
  }
}
else
  warning("It was not possible to add the correlation to the regulonTargetsInfo table.")
linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists = "null")
if (!is.null(linkList) & ("weight" %in% colnames(linkList))) {
  if (is.data.table(linkList))
    linkList <- as.data.frame(linkList)
  uniquePairs <- nrow(unique(linkList[, c("TF", "Target")]))
  if (uniquePairs == nrow(linkList)) {
    linkList <-
      linkList[which(linkList$weight >= getSettings(scenicOptions,
                                                    "modules/weightThreshold")),]
    rownames(linkList) <- paste(linkList$TF, linkList$Target,
                                sep = "__")
    regulonTargetsInfo <- cbind(regulonTargetsInfo,
                                CoexWeight = linkList[paste(regulonTargetsInfo$TF,
                                                            regulonTargetsInfo$gene, sep = "__"), "weight"])
  }
  else {
    warning(
      "There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.",
      "\nThe co-expression weight was not added to the regulonTargetsInfo table."
    )
  }
}
else
  warning("It was not possible to add the weight to the regulonTargetsInfo table.")
saveRDS(regulonTargetsInfo,
        file = getIntName(scenicOptions,
                          "regulonTargetsInfo"))
write.table(
  regulonTargetsInfo,
  file = getOutName(scenicOptions,
                    "s2_regulonTargetsInfo"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)
rm(linkList)
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo,
                                         regulonTargetsInfo$highConfAnnot)
regulons <- NULL
if (!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]])) {
  regulons <- lapply(split(
    regulonTargetsInfo_splitByAnnot[["TRUE"]],
    regulonTargetsInfo_splitByAnnot[["TRUE"]][, "TF"]
  ),
  function(x)
    sort(as.character(unlist(x[, "gene"]))))
}
regulons_extended <- NULL
if (!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]])) {
  regulons_extended <-
    lapply(split(
      regulonTargetsInfo_splitByAnnot[["FALSE"]],
      regulonTargetsInfo_splitByAnnot[["FALSE"]][, "TF"]
    ),
    function(x)
      unname(unlist(x[, "gene"])))
  regulons_extended <- setNames(lapply(names(regulons_extended),
                                       function(tf)
                                         sort(unique(
                                           c(regulons[[tf]], unlist(regulons_extended[[tf]]))
                                         ))),
                                names(regulons_extended))
  names(regulons_extended) <- paste(names(regulons_extended),
                                    "_extended", sep = "")
}
regulons <- c(regulons, regulons_extended)
saveRDS(regulons, file = getIntName(scenicOptions, "regulons"))
incidList <- reshape2::melt(regulons)
incidMat <- table(incidList[, 2], incidList[, 1])
saveRDS(incidMat, file = getIntName(scenicOptions, "regulons_incidMat"))
rm(incidMat)
if (getSettings(scenicOptions, "verbose")) {
  length(regulons)
  summary(lengths(regulons))
}
scenicOptions@status$current <- 2
invisible(scenicOptions)

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(Seurat)
library(ggpubr)
library(paletteer)
rm(list = ls())
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
  rename(protein = Accessions)

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
data_tran <- scdata_tran@assays$RNA$counts %>% as.data.frame()
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


data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
res_tf <- data.frame(
  tf_name = c(
    data_raw_1$Symbol,data_raw_2$tf
  ) %>% unique()
)




vlnplot_pro <- function(plottype = 'first',scdata_pro,gene_select,celltype = 'Oocyte', num_col = NULL){
  mycolor_use <- if(celltype =='Oocyte'){my36colors[5:8]}else{my36colors[1:4]}
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
      expression = log10(expression+1),
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




gene_select <- c("Mideas","Zbed6")
vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Granulosa',
  gene_select = gene_select,
  num_col = 11
)
ggsave('./TF_barplot_second.png',width = 10,height = 6)
ggsave('./TF_barplot_second.pdf',width = 10,height = 6)


data_plot_tran <- data_tran %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rename_all(~c('abundance'))  %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% res_tf$tf_name) %>% 
  mutate(percentage = abundance/sum(abundance)) %>% 
  arrange(percentage) %>% 
  mutate(
    cumulative_intensity = cumsum(percentage),
    group_omits = 'transcription',
    rank = 1:nrow(.)
  )
data_plot_pro <- data_singlestage %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rename_all(~c('abundance')) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% res_tf$tf_name) %>% 
  mutate(percentage = abundance/sum(abundance)) %>% 
  arrange(percentage) %>% 
  mutate(
    cumulative_intensity = cumsum(percentage),
    group_omits = 'protein',
    rank = 1:nrow(.)
  )
data_plot <- data_plot_pro
data_plor_label <- data_plot %>% 
  arrange(group_omits,cumulative_intensity) %>% 
  group_by(group_omits) %>% 
  dplyr::filter(Genes %in% c("Mideas","Zbed6"))
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
    color = 'grey10',nudge_y = -0.3,
    show.legend = FALSE,
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0.8,
    segment.size  = 1,force_pull = 6,
    force = 900,
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
  './protein_abundance_tf.png',
  width = 6,height = 8
)
ggsave(
  './protein_abundance_tf.pdf',
  width = 5,height = 8
)


data_plot_pro <- data_singlestage %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rename_all(~c('abundance')) %>% 
  rownames_to_column('Genes') %>% 
  filter(Genes %in% res_tf$tf_name) %>% 
  mutate(
    percentage = abundance/sum(abundance),
    abundance_log = log10(abundance + 1)
  ) %>% 
  filter(abundance_log>0) %>% 
  arrange(abundance_log) %>% 
  mutate(
    rank = 1:nrow(.),
    group_omits = 'protein'
  )
data_plot <- data_plot_pro
data_plor_label <- data_plot %>% 
  arrange(group_omits,abundance_log) %>% 
  group_by(group_omits) %>% 
  dplyr::filter(Genes %in% c("Mideas","Zbed6"))
ggplot(data = data_plot, aes(x = rank, y = abundance_log, color = group_omits)) +
  geom_point() +
  geom_point(
    data = data_plor_label,
    aes(x = rank, y = abundance_log),
    color = 'grey10',
    size = 2
  ) +
  geom_text_repel(
    data = data_plor_label,
    aes(x = rank, y = abundance_log, label = Genes),
    color = 'grey10',nudge_y = -2,
    show.legend = FALSE,
    box.padding = 0.5,
    size = 6,
    segment.curvature = 0.6,
    segment.size  = 1,
    force = 900,
    fontface = "italic",
    seed = 42,
    arrow = arrow(length = unit(0.03, "npc")),
    max.overlaps = Inf
  ) +
  facet_wrap(~ group_omits, nrow = 1, scales = 'free_x') +
  scale_color_manual(values = '#C02B40') +
  labs(x = 'Abundance rank', y = 'Log10(Abundance)') +
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

gene_select <- c('Sub1','Gnrh1','Hmgb2')
vlnplot_pro(
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select,
  num_col = 11)
ggsave('./TF_barplot_Oocyte.pdf',width = 10,height = 6)
vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select,
  num_col = 11)
ggsave('./TF_barplot_Oocyte_second.pdf',width = 10,height = 6)
gene_select <- c('RHOT1','MGARP','UB','SATB1') %>% str_to_title()
vlnplot_pro(
  scdata_pro = scdata_pro,
  celltype = 'Oocyte',
  gene_select = gene_select,
  num_col = 11)
ggsave('./TF_barplot_Oocyte_IHC.pdf',width = 10,height = 6)

gene_select <- c('Cebpb','Foxo1','Stat3')
vlnplot_pro(
  scdata_pro = scdata_pro,
  celltype = 'Granulosa',
  gene_select = gene_select,
  num_col = 11)
ggsave('./TF_barplot_GC.pdf',width = 10,height = 6)
vlnplot_pro(
  plottype = 'second',
  scdata_pro = scdata_pro,
  celltype = 'Granulosa',
  gene_select = gene_select,
  num_col = 11)
ggsave('./TF_barplot_GC_second.pdf',width = 10,height = 6)

library(Seurat)
library(Biobase)
library(dplyr)
library(viper)
library(tidyr)
library(tibble)
setwd('./diann_copynumber/')
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
  rename(Accessions = ProteinGroups)
protein_list <- data_pro_clean[,c('Genes','Accessions')] %>% 
  dplyr::rename(protein = Accessions)
protein_list %>% 
  filter(Genes %in% c('Wt1'))

scdata <- readRDS('./scdata_use.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata$group_cell %>% unique()
scdata_use <- subset(scdata,subset = group_cell == 'Oocyte')
data_regulon <- read.csv('./network.txt',sep = '\t')
tf_list <- data_regulon$Regulator %>% unique()
protein_id <- protein_list$protein[protein_list$Genes %in% tf_list]
rownames(scdata_use)
data_exp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  {as.data.frame(t(apply(., 1, replace_outliers)))} %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  left_join(scdata@meta.data %>% select(group_stage,group_cell) %>% rownames_to_column('sample'),by = 'sample') %>% 
  filter(group_cell == 'Oocyte') %>% 
  select(-c('sample','group_cell')) %>% 
  group_by(group_stage) %>% 
  summarise(across(everything(),mean)) %>% 
  column_to_rownames('group_stage') %>% 
  t() %>% as.data.frame() %>% filter(rowSums(.)>0)
data_exp %>%
  rownames_to_column('Genes') %>%
  filter(Genes == 'P22561') %>%
  t()

standardize <- function(x) {
  (x - mean(x)) / sd(x)
}
data_standardized <- apply(data_exp, 1, standardize) %>% 
  t() %>% as.data.frame()
data_standardized[rownames(data_standardized) == 'O08609',]


library(cluster)


sil_width <- numeric(50)
for (i in 2:50) {
  km <- kmeans(data_standardized, centers = i, nstart = 20)
  ss <- silhouette(km$cluster, dist(data_standardized))
  sil_width[i] <- mean(ss[, 3])
}


silhouette_plot <- data.frame(k = 2:50, sil_width = sil_width[2:50])
ggplot(silhouette_plot, aes(x = k, y = sil_width)) +
  geom_point() +
  geom_line() +
  labs(title = "Silhouette Method for Determining Optimal Number of Clusters",
       x = "Number of Clusters",
       y = "Average Silhouette Width") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(angle = 30,length = unit(3,'mm'))),
    axis.text.x = element_text(size = 18,angle = -45,hjust = 0,vjust = 0.5),
    axis.text.y = element_text(),
    axis.title = element_text(size = 18),
    legend.position = 'none',
    plot.margin = margin(l = 2,r = 5,unit = 'mm')
  )


k = 28

kmeans_result <- kmeans(data_standardized, centers = k)


data_standardized$cluster <- kmeans_result$cluster


colnames(data_standardized)


data_long <- data_standardized %>%
  rownames_to_column('Regulon') %>% 
  pivot_longer(cols = names(.)[2:5], names_to = "Experiment", values_to = "Value") %>% 
  mutate(
    Experiment = factor(Experiment,levels = c("Secondary","Early antral","Antral","Preovulatory")),
    cluster = paste('Cluster',cluster,sep = ' '),
    cluster = factor(cluster,levels = paste('Cluster',1:k,sep = ' '))
  )

library(paletteer)
library(tidyr)
library(tibble)

ggplot(data_long, aes(x = Experiment, y = Value, color = factor(cluster), group = Regulon)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = 
                       c(
                         paletteer_d("ggsci::planetexpress_futurama"),
                         paletteer_d("ggsci::springfield_simpsons"))
  ) +
  facet_wrap(~cluster,ncol = 6,axes = 'all') +
  labs(title = "K-means Clustering of Regulons Based on scaled NES score",
       x = "Experiment",
       y = "Expression Level",
       color = "Cluster") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    axis.line = element_line(arrow = arrow(angle = 30,length = unit(3,'mm'))),
    axis.text.x = element_text(size = 18,angle = -45,hjust = 0,vjust = 0.5),
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    legend.position = 'none',
    plot.title = element_text(size = 22,hjust = 0.5),
    plot.margin = margin(l = 2,r = 15,unit = 'mm')
  )
ggsave('./data_res_tf_diann_cn_Oc.png',bg = 'white',width = 40,height = 20)
data_res <- data_exp %>% 
  mutate(cluster = kmeans_result$cluster) %>% 
  rownames_to_column('protein') %>% 
  left_join(protein_list,by = 'protein') %>% 
  select(names(.)[1],'Genes',names(.)[2:6]) %>% 
  rowwise() %>%
  mutate(
    express_sum = sum(c_across(3:6))
  ) %>% 
  ungroup() %>%
  mutate(
    quantile_group = ntile(express_sum, 4) %>% factor(.,levels = rev(1:4))
  ) %>% 
  group_by(quantile_group) %>% 
  arrange(desc(cluster), .by_group = TRUE)
write.csv(data_res,'./data_res_diann_cn.csv')
data_tmp <- data_exp %>% 
  rowwise() %>%
  mutate(
    express_sum = sum(c_across(1:4))
  )



data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
res_tf <- data.frame(
  Genes = c(
    data_raw_1$Symbol,data_raw_2$tf
  ) %>% unique(),
  group = 'tf'
)
colnames(data_res)
data_res_tf <- data_res %>% 
  left_join(res_tf,by = 'Genes') %>%
  filter(group == 'tf') %>%
  arrange(group,cluster)
write.csv(data_res_tf,'./data_res_tf_diann_cn_Oc.csv')




library(Seurat)
library(Biobase)
library(dplyr)
library(viper)
library(tidyr)
library(tibble)
setwd('./diann_copynumber/')
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
  rename(Accessions = ProteinGroups)
protein_list <- data_pro_clean[,c('Genes','Accessions')] %>% 
  dplyr::rename(protein = Accessions)
protein_list %>% 
  filter(Genes %in% c('Wt1'))

scdata <- readRDS('./scdata_use.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata$group_cell %>% unique()
data_regulon <- read.csv('./network.txt',sep = '\t')
tf_list <- data_regulon$Regulator %>% unique()
protein_id <- protein_list$protein[protein_list$Genes %in% tf_list]
rownames(scdata_use)
data_exp <- scdata@assays$RNA$counts %>% as.data.frame() %>% 
  {as.data.frame(t(apply(., 1, replace_outliers)))} %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  left_join(scdata@meta.data %>% select(group_stage,group_cell) %>% rownames_to_column('sample'),by = 'sample') %>% 
  filter(group_cell == 'Granulosa') %>% 
  select(-c('sample','group_cell')) %>% 
  group_by(group_stage) %>% 
  summarise(across(everything(),mean)) %>% 
  column_to_rownames('group_stage') %>% 
  t() %>% as.data.frame() %>% filter(rowSums(.)>0)
data_exp %>%
  rownames_to_column('Genes') %>%
  filter(Genes == 'P22561') %>%
  t()

standardize <- function(x) {
  (x - mean(x)) / sd(x)
}
data_standardized <- apply(data_exp, 1, standardize) %>% 
  t() %>% as.data.frame()
data_standardized[rownames(data_standardized) == 'O08609',]


library(cluster)


sil_width <- numeric(50)
for (i in 2:50) {
  km <- kmeans(data_standardized, centers = i, nstart = 20)
  ss <- silhouette(km$cluster, dist(data_standardized))
  sil_width[i] <- mean(ss[, 3])
}


silhouette_plot <- data.frame(k = 2:50, sil_width = sil_width[2:50])
ggplot(silhouette_plot, aes(x = k, y = sil_width)) +
  geom_point() +
  geom_line() +
  labs(title = "Silhouette Method for Determining Optimal Number of Clusters",
       x = "Number of Clusters",
       y = "Average Silhouette Width") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(angle = 30,length = unit(3,'mm'))),
    axis.text.x = element_text(size = 18,angle = -45,hjust = 0,vjust = 0.5),
    axis.text.y = element_text(),
    axis.title = element_text(size = 18),
    legend.position = 'none',
    plot.margin = margin(l = 2,r = 5,unit = 'mm')
  )


k = 28

kmeans_result <- kmeans(data_standardized, centers = k)


data_standardized$cluster <- kmeans_result$cluster


colnames(data_standardized)


data_long <- data_standardized %>%
  rownames_to_column('Regulon') %>% 
  pivot_longer(cols = names(.)[2:5], names_to = "Experiment", values_to = "Value") %>% 
  mutate(
    Experiment = factor(Experiment,levels = c("Secondary","Early antral","Antral","Preovulatory")),
    cluster = paste('Cluster',cluster,sep = ' '),
    cluster = factor(cluster,levels = paste('Cluster',1:k,sep = ' '))
  )

library(paletteer)
library(tidyr)
library(tibble)

ggplot(data_long, aes(x = Experiment, y = Value, color = factor(cluster), group = Regulon)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = 
                       c(
                         paletteer_d("ggsci::planetexpress_futurama"),
                         paletteer_d("ggsci::springfield_simpsons"))
  ) +
  facet_wrap(~cluster,ncol = 6,axes = 'all') +
  labs(title = "K-means Clustering of Regulons Based on scaled NES score",
       x = "Experiment",
       y = "Expression Level",
       color = "Cluster") +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 20),
    axis.line = element_line(arrow = arrow(angle = 30,length = unit(3,'mm'))),
    axis.text.x = element_text(size = 18,angle = -45,hjust = 0,vjust = 0.5),
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    legend.position = 'none',
    plot.title = element_text(size = 22,hjust = 0.5),
    plot.margin = margin(l = 2,r = 15,unit = 'mm')
  )
ggsave('./data_res_tf_diann_cn_Gc.png',bg = 'white',width = 40,height = 20)
data_res <- data_exp %>% 
  mutate(cluster = kmeans_result$cluster) %>% 
  rownames_to_column('protein') %>% 
  left_join(protein_list,by = 'protein') %>% 
  select(names(.)[1],'Genes',names(.)[2:6]) %>% 
  rowwise() %>%
  mutate(
    express_sum = sum(c_across(3:6))
  ) %>% 
  ungroup() %>%
  mutate(
    quantile_group = ntile(express_sum, 4) %>% factor(.,levels = rev(1:4))
  ) %>% 
  group_by(quantile_group) %>% 
  arrange(desc(cluster), .by_group = TRUE)
write.csv(data_res,'./data_res_diann_cn_Gc.csv')
data_tmp <- data_exp %>% 
  rowwise() %>%
  mutate(
    express_sum = sum(c_across(1:4))
  )



data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
res_tf <- data.frame(
  Genes = c(
    data_raw_1$Symbol,data_raw_2$tf
  ) %>% unique(),
  group = 'tf'
)
colnames(data_res)
data_res_tf <- data_res %>% 
  left_join(res_tf,by = 'Genes') %>%
  filter(group == 'tf') %>%
  arrange(group,cluster)
write.csv(data_res_tf,'./data_res_tf_diann_cn_Gc.csv')




library(dplyr)
library(tibble)
library(Seurat)
library(scales)
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
)

my_other36colors <- sapply(my36colors, function(col) {
  col <- col2rgb(col) / 255
  col <- pmin(col * 1.2, 1)
  col <- rgb(col[1], col[2], col[3])
  return(col)
})


my72colors <- c(my36colors, my_other36colors)
show_col(my72colors[37:72])
show_col(my72colors[1:36])
data_raw <- read.csv('./Mus_musculus_TF.txt',sep = '\t')
Family_tf <- data_raw$Family %>% unique()
my36colors <- my72colors[1:length(Family_tf)]
names(my36colors) <- Family_tf

tmp <- data_raw$Family %>% table() %>%
  as.data.frame() %>% arrange(desc(Freq)) %>% 
  rename_all(~c('Family','Number')) %>% 
  mutate(Family = factor(Family,levels = Family))
ggplot(data = tmp,aes(Family,Number,fill = Family)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = 'black'),
    legend.position = 'none',
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
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
  rename(protein = Accessions)

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


data_raw %>% colnames()
data_raw$Family %>% table()
scdata_pro_sub <- subset(scdata_pro,subset = group_cell =='Oocyte')
gene_useful <- scdata_pro_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot %>% colnames()

ggplot(data_plot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = Family)) +
  geom_text(aes(x = 3.5, y = ylab, label = lab, angle = angle), hjust = 0.5, size = 5) +
  scale_fill_manual(values = my36colors) +
  xlim(2, 4) +
  coord_polar(theta = "y") +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )
file_path <- './protein_Oocyte.png'
ggsave(file_path,width = 9,height = 9)
file_path <- './protein_Oocyte.pdf'
ggsave(file_path,width = 9,height = 9)

scdata_pro_sub <- subset(scdata_pro,subset = group_cell =='Granulosa')
gene_useful <- scdata_pro_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot %>% colnames()

ggplot(data_plot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = Family)) +
  geom_text(aes(x = 3.5, y = ylab, label = lab, angle = angle), hjust = 0.5, size = 5) +
  scale_fill_manual(values = my36colors) +
  xlim(2, 4) +
  coord_polar(theta = "y") +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )
file_path <- './protein_Granulosa.png'
ggsave(file_path,width = 9,height = 9)
file_path <- './protein_Granulosa.pdf'
ggsave(file_path,width = 9,height = 9)

scdata_trans_sub <- subset(scdata_tran,subset = group_cell =='Oocyte')
gene_useful <- scdata_trans_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot %>% colnames()

ggplot(data_plot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = Family)) +
  geom_text(aes(x = 3.5, y = ylab, label = lab, angle = angle), hjust = 0.5, size = 5) +
  scale_fill_manual(values = my36colors) +
  xlim(2, 4) +
  coord_polar(theta = "y") +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )
file_path <- './transcript_Oocyte.png'
ggsave(file_path,width = 9,height = 9)
file_path <- './transcript_Oocyte.pdf'
ggsave(file_path,width = 9,height = 9)

library(dplyr)
data_raw %>% colnames()
data_raw$Family %>% table()
scdata_pro_sub <- subset(scdata_pro,subset = group_cell =='Oocyte')
gene_useful <- scdata_pro_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot_protein_O <- data_plot

scdata_trans_sub <- subset(scdata_tran,subset = group_cell =='Oocyte')
gene_useful <- scdata_trans_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot_trans_O <- data_plot
df_counts <- bind_rows(
  data_plot_protein_O %>% mutate(omic = "protein"),
  data_plot_trans_O   %>% mutate(omic = "trans")
) %>% select(Family, Counts, omic)

tab <- xtabs(Counts ~ Family + omic, data = df_counts)

chisq <- chisq.test(tab)
chisq <- chisq.test(tab, simulate.p.value = TRUE, B = 100000)
enrich_df <- as.data.frame.matrix(chisq$observed / chisq$expected)
enrich_df
std_resid <- chisq$stdres %>% as.data.frame() %>% 
  pivot_wider(id_cols = 'Family',names_from = 'omic',values_from = 'Freq') %>% 
  column_to_rownames('Family')# standardized residuals matrix (Family x omic)
std_resid
pearson_resid <- chisq$residuals
pearson_resid

fisher_2x2_safe <- function(m, B = 100000, seed = 1) {
  stopifnot(is.matrix(m), all(dim(m) == c(2, 2)))
  if (any(m < 0) || any(is.na(m))) stop("Invalid counts in 2x2 table.")
  if (sum(m) == 0) {
    return(list(or = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  }
  
  ft <- tryCatch(
    fisher.test(m),
    error = function(e) {
      set.seed(seed)
      fisher.test(m, simulate.p.value = TRUE, B = B)
    }
  )
  
  a <- m[1, 1]; b <- m[1, 2]; c <- m[2, 1]; d <- m[2, 2]
  a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
  or <- (a2 * d2) / (b2 * c2)
  se_log_or <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
  ci_low  <- exp(log(or) - 1.96 * se_log_or)
  ci_high <- exp(log(or) + 1.96 * se_log_or)
  
  list(or = or, ci_low = ci_low, ci_high = ci_high, p = unname(ft$p.value))
}

per_family_or_fdr <- function(tab, B = 100000, seed = 1) {
  stopifnot(is.matrix(tab), ncol(tab) == 2)
  if (is.null(colnames(tab))) colnames(tab) <- c("protein", "trans")
  
  col1 <- colnames(tab)[1]  # e.g. protein
  col2 <- colnames(tab)[2]  # e.g. trans
  
  tot1 <- sum(tab[, col1])
  tot2 <- sum(tab[, col2])
  
  res <- map_dfr(rownames(tab), function(fam) {
    a <- tab[fam, col1]                 # family in protein
    c <- tab[fam, col2]                 # family in trans
    b <- tot1 - a                       # others in protein
    d <- tot2 - c                       # others in trans
    
    m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c(col1, col2), c("family", "others")))
    
    ft <- fisher_2x2_safe(m, B = B, seed = seed)
    
    tibble(
      Family = fam,
      a_family_in_protein = as.integer(a),
      b_others_in_protein = as.integer(b),
      c_family_in_trans   = as.integer(c),
      d_others_in_trans   = as.integer(d),
      OR = ft$or,
      CI95_low = ft$ci_low,
      CI95_high = ft$ci_high,
      p_value = ft$p
    )
  }) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value)
  
  res
}


or_fdr_tbl <- per_family_or_fdr(tab, B = 100000, seed = 1)

or_fdr_tbl <- or_fdr_tbl %>%
  mutate(
    log2OR = log2(OR),
    direction = case_when(
      is.na(OR) ~ NA_character_,
      OR > 1 ~ "Enriched_in_protein",
      OR < 1 ~ "Enriched_in_trans",
      TRUE   ~ "No_difference"
    )
  )

or_fdr_tbl
write.csv(or_fdr_tbl,'./tf_Oocyte_fdr.csv')


family_levels <- enrich_df %>% 
  arrange(desc(protein)) %>% 
  rownames()

data_ggplot <- enrich_df %>%
  rownames_to_column("Family") %>%
  pivot_longer(cols = -Family, names_to = "omic", values_to = "RoE") %>%
  left_join(or_fdr_tbl %>% select(Family, FDR), by = "Family") %>%
  mutate(
    omic = factor(omic, levels = c("protein", "trans")),
    Family = factor(Family, levels = family_levels),
    

    signif = case_when(
      omic != "protein" ~ "",
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    

    label = if_else(signif == "", sprintf("%.2f", RoE), paste0(sprintf("%.2f", RoE), signif))
  )

limit_min <- min(data_ggplot$RoE, na.rm = TRUE)
limit_max <- max(data_ggplot$RoE, na.rm = TRUE)

p <- ggplot(data_ggplot, aes(x = Family, y = omic)) +
  geom_tile(aes(fill = RoE), color = "white", linewidth = 0.2) +
  geom_text(aes(label = signif), size = 3.2) +
  guides(fill = guide_colorbar(title = "Ro/e", title.vjust = 1)) +
  scale_fill_gradient2(
    low = "#418EA5", mid = "#FCFADD", high = "#CE584B",midpoint = 1,
    limits = c(limit_min, limit_max),
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "right"
  )
p
file_path <- './tf_Oocyte_fisher.pdf'
cairo_pdf(file_path,width = 16,height = 3)
p
dev.off()

library(patchwork)
colnames(or_fdr_tbl)
dfp <- or_fdr_tbl %>%

filter(!is.na(FDR) & FDR < 0.05) %>%

  mutate(
    log2OR = log2(OR),
    log2CI_low  = log2(CI95_low),
    log2CI_high = log2(CI95_high),
    direction = if_else(log2OR >= 0, "Enriched_in_protein", "Enriched_in_trans"),
    P_txt = formatC(FDR, format = "e", digits = 2),
    OR_CI_txt = sprintf("%.2f (%.2f–%.2f)", OR, CI95_low, CI95_high)
  ) %>%
  arrange(desc(log2OR),FDR) %>%
  mutate(Family = factor(Family, levels = rev(Family)))


n_row <- nlevels(dfp$Family)

x_family <- 0
x_p      <- 1.2
x_orci   <- 2.7

p_tbl <- ggplot(dfp, aes(y = Family)) +

  geom_hline(yintercept = n_row + 0.5, linewidth = 0.6, color = "black") +
  

  geom_text(aes(x = x_family, label = as.character(Family)), hjust = 0, size = 3.6) +
  geom_text(aes(x = x_p,      label = P_txt),               hjust = 0, size = 3.6) +
  geom_text(aes(x = x_orci,   label = OR_CI_txt),           hjust = 0, size = 3.6) +
  

  annotate("text", x = x_p,    y = n_row + 0.9, label = "FDR (BH)",
           hjust = 0, fontface = "bold", size = 3.8) +
  annotate("text", x = x_orci, y = n_row + 0.9, label = "Odds ratio (95% CI)",
           hjust = 0, fontface = "bold", size = 3.8) +
  
  scale_x_continuous(limits = c(-0.1, 4.6), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.6))) +
  coord_cartesian(clip = "off") +
  

  theme_classic() +
  theme(

    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_line(linewidth = 0.6, color = "black"),
    

    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )


p_for <- ggplot(dfp, aes(y = Family)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "black") +
  geom_errorbarh(aes(xmin = log2CI_low, xmax = log2CI_high),
                 height = 0.18, linewidth = 0.6, color = "black") +
  geom_point(aes(x = log2OR, fill = direction),
             shape = 21, size = 3.2, stroke = 0.4, color = "black") +
  scale_fill_manual(values = c(
    Enriched_in_protein = "#D55E00",
    Enriched_in_trans   = "#009E73"
  ), guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.6))) +
  scale_x_continuous(
    limits = c(-4, 8),
    breaks = c(-4, 0, 4, 8),
    labels = c(-4, 0, 4, 8)
  ) +#
  labs(x = "log2(Odds ratio)", y = NULL) +
  theme_classic() +
  theme(

    

    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    

    axis.line.x  = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )


(p_tbl + p_for) + plot_layout(widths = c(2, 1.52))
file_path <- './tf_Oocyte_forset.pdf'
cairo_pdf(file_path,width = 10,height = 5)
(p_tbl + p_for) + plot_layout(widths = c(2, 1.52))
dev.off()


scdata_trans_sub <- subset(scdata_tran,subset = group_cell =='Granulosa')
gene_useful <- scdata_trans_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot %>% colnames()

ggplot(data_plot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(aes(fill = Family)) +
  geom_text(aes(x = 3.5, y = ylab, label = lab, angle = angle), hjust = 0.5, size = 5) +
  scale_fill_manual(values = my36colors) +
  xlim(2, 4) +
  coord_polar(theta = "y") +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = 'none'
  )


file_path <- './transcript_Granulosa.png'
ggsave(file_path,width = 9,height = 9)

file_path <- './transcript_Granulosa.pdf'
ggsave(file_path,width = 9,height = 9)

library(dplyr)
data_raw %>% colnames()
data_raw$Family %>% table()
scdata_pro_sub <- subset(scdata_pro,subset = group_cell =='Granulosa')
gene_useful <- scdata_pro_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  filter(Family != 'Others') %>% 
  slice(1:10) %>%
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot_protein_O <- data_plot

scdata_trans_sub <- subset(scdata_tran,subset = group_cell =='Granulosa')
gene_useful <- scdata_trans_sub@assays$RNA$counts %>% 
  as.data.frame() %>% 
  rownames_to_column('Genes') %>% 
  rowwise() %>%
  mutate(count_gt_0 = sum(c_across(where(is.numeric)) > 0)) %>%
  ungroup() %>%  
  filter(count_gt_0>0)
data_plot <- data_raw %>% 
  filter(Symbol %in% gene_useful$Genes) %>% 
  group_by(Family) %>% 
  summarise(Counts = n(),.groups = 'drop') %>% 
  arrange(desc(Counts)) %>% 
  mutate(
    percent = Counts/sum(Counts),
    ymax = cumsum(percent),
    ymin = c(0,ymax[-nrow(.)]),
    ylab = (ymax+ymin)/2,
    lab = paste(Family,'(',round(percent * 100,digits = 2), '%)',sep = ''),
    midpoint = (ymin + ymax) / 2,
    angle = 90 - (midpoint / max(midpoint)) * 360,  # Adjust angle to point towards center
    angle = ifelse(angle < -90, angle + 180, angle)
  )
data_plot_trans_O <- data_plot
df_counts <- bind_rows(
  data_plot_protein_O %>% mutate(omic = "protein"),
  data_plot_trans_O   %>% mutate(omic = "trans")
) %>% select(Family, Counts, omic)

tab <- xtabs(Counts ~ Family + omic, data = df_counts)

chisq <- chisq.test(tab)
chisq <- chisq.test(tab, simulate.p.value = TRUE, B = 100000)
enrich_df <- as.data.frame.matrix(chisq$observed / chisq$expected)
enrich_df
std_resid <- chisq$stdres %>% as.data.frame() %>% 
  pivot_wider(id_cols = 'Family',names_from = 'omic',values_from = 'Freq') %>% 
  column_to_rownames('Family')# standardized residuals matrix (Family x omic)
std_resid
pearson_resid <- chisq$residuals
pearson_resid

fisher_2x2_safe <- function(m, B = 100000, seed = 1) {
  stopifnot(is.matrix(m), all(dim(m) == c(2, 2)))
  if (any(m < 0) || any(is.na(m))) stop("Invalid counts in 2x2 table.")
  if (sum(m) == 0) {
    return(list(or = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  }
  
  ft <- tryCatch(
    fisher.test(m),
    error = function(e) {
      set.seed(seed)
      fisher.test(m, simulate.p.value = TRUE, B = B)
    }
  )
  
  a <- m[1, 1]; b <- m[1, 2]; c <- m[2, 1]; d <- m[2, 2]
  a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
  or <- (a2 * d2) / (b2 * c2)
  se_log_or <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
  ci_low  <- exp(log(or) - 1.96 * se_log_or)
  ci_high <- exp(log(or) + 1.96 * se_log_or)
  
  list(or = or, ci_low = ci_low, ci_high = ci_high, p = unname(ft$p.value))
}

per_family_or_fdr <- function(tab, B = 100000, seed = 1) {
  stopifnot(is.matrix(tab), ncol(tab) == 2)
  if (is.null(colnames(tab))) colnames(tab) <- c("protein", "trans")
  
  col1 <- colnames(tab)[1]  # e.g. protein
  col2 <- colnames(tab)[2]  # e.g. trans
  
  tot1 <- sum(tab[, col1])
  tot2 <- sum(tab[, col2])
  
  res <- map_dfr(rownames(tab), function(fam) {
    a <- tab[fam, col1]                 # family in protein
    c <- tab[fam, col2]                 # family in trans
    b <- tot1 - a                       # others in protein
    d <- tot2 - c                       # others in trans
    
    m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c(col1, col2), c("family", "others")))
    
    ft <- fisher_2x2_safe(m, B = B, seed = seed)
    
    tibble(
      Family = fam,
      a_family_in_protein = as.integer(a),
      b_others_in_protein = as.integer(b),
      c_family_in_trans   = as.integer(c),
      d_others_in_trans   = as.integer(d),
      OR = ft$or,
      CI95_low = ft$ci_low,
      CI95_high = ft$ci_high,
      p_value = ft$p
    )
  }) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value)
  
  res
}


or_fdr_tbl <- per_family_or_fdr(tab, B = 100000, seed = 1)

or_fdr_tbl <- or_fdr_tbl %>%
  mutate(
    log2OR = log2(OR),
    direction = case_when(
      is.na(OR) ~ NA_character_,
      OR > 1 ~ "Enriched_in_protein",
      OR < 1 ~ "Enriched_in_trans",
      TRUE   ~ "No_difference"
    )
  )

or_fdr_tbl
write.csv(or_fdr_tbl,'./tf_Granulosa_fdr.csv')


family_levels <- enrich_df %>% 
  arrange(desc(protein)) %>% 
  rownames()

data_ggplot <- enrich_df %>%
  rownames_to_column("Family") %>%
  pivot_longer(cols = -Family, names_to = "omic", values_to = "RoE") %>%
  left_join(or_fdr_tbl %>% select(Family, FDR), by = "Family") %>%
  mutate(
    omic = factor(omic, levels = c("protein", "trans")),
    Family = factor(Family, levels = family_levels),
    

    signif = case_when(
      omic != "protein" ~ "",
      !is.na(FDR) & FDR < 0.001 ~ "***",
      !is.na(FDR) & FDR < 0.01  ~ "**",
      !is.na(FDR) & FDR < 0.05  ~ "*",
      TRUE ~ ""
    ),
    

    label = if_else(signif == "", sprintf("%.2f", RoE), paste0(sprintf("%.2f", RoE), signif))
  )

limit_min <- min(data_ggplot$RoE, na.rm = TRUE)
limit_max <- max(data_ggplot$RoE, na.rm = TRUE)

p <- ggplot(data_ggplot, aes(x = Family, y = omic)) +
  geom_tile(aes(fill = RoE), color = "white", linewidth = 0.2) +
  geom_text(aes(label = signif), size = 3.2) +
  guides(fill = guide_colorbar(title = "Ro/e", title.vjust = 1)) +
  scale_fill_gradient2(
    low = "#418EA5", mid = "#FCFADD", high = "#CE584B",midpoint = 1,
    limits = c(limit_min, limit_max),
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "right"
  )
p
file_path <- './tf_Granulosa_fisher.pdf'
cairo_pdf(file_path,width = 16,height = 3)
p
dev.off()

library(patchwork)
colnames(or_fdr_tbl)
dfp <- or_fdr_tbl %>%

  filter(!is.na(FDR) & FDR < 0.05) %>%

  mutate(
    log2OR = log2(OR),
    log2CI_low  = log2(CI95_low),
    log2CI_high = log2(CI95_high),
    direction = if_else(log2OR >= 0, "Enriched_in_protein", "Enriched_in_trans"),
    P_txt = formatC(FDR, format = "e", digits = 2),
    OR_CI_txt = sprintf("%.2f (%.2f–%.2f)", OR, CI95_low, CI95_high)
  ) %>%
  arrange(desc(log2OR),FDR) %>%
  mutate(Family = factor(Family, levels = rev(Family)))


n_row <- nlevels(dfp$Family)

x_family <- 0
x_p      <- 1.2
x_orci   <- 2.7

p_tbl <- ggplot(dfp, aes(y = Family)) +

  geom_hline(yintercept = n_row + 0.5, linewidth = 0.6, color = "black") +
  

  geom_text(aes(x = x_family, label = as.character(Family)), hjust = 0, size = 3.6) +
  geom_text(aes(x = x_p,      label = P_txt),               hjust = 0, size = 3.6) +
  geom_text(aes(x = x_orci,   label = OR_CI_txt),           hjust = 0, size = 3.6) +
  

  annotate("text", x = x_p,    y = n_row + 0.9, label = "FDR (BH)",
           hjust = 0, fontface = "bold", size = 3.8) +
  annotate("text", x = x_orci, y = n_row + 0.9, label = "Odds ratio (95% CI)",
           hjust = 0, fontface = "bold", size = 3.8) +
  
  scale_x_continuous(limits = c(-0.1, 4.6), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.6))) +
  coord_cartesian(clip = "off") +
  

  theme_classic() +
  theme(

    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_line(linewidth = 0.6, color = "black"),
    

    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )


p_for <- ggplot(dfp, aes(y = Family)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "black") +
  geom_errorbarh(aes(xmin = log2CI_low, xmax = log2CI_high),
                 height = 0.18, linewidth = 0.6, color = "black") +
  geom_point(aes(x = log2OR, fill = direction),
             shape = 21, size = 3.2, stroke = 0.4, color = "black") +
  scale_fill_manual(values = c(
    Enriched_in_protein = "#D55E00",
    Enriched_in_trans   = "#009E73"
  ), guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.6))) +
  scale_x_continuous(
    limits = c(-10, 5),
    breaks = c(-10,-5, 0, 5),
    labels = c(-10,-5, 0, 5)
  ) +#
  labs(x = "log2(Odds ratio)", y = NULL) +
  theme_classic() +
  theme(

    

    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    

    axis.line.x  = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )


(p_tbl + p_for) + plot_layout(widths = c(2, 1.52))
file_path <- './tf_Granulosa_forset.pdf'
cairo_pdf(file_path,width = 10,height = 5)
(p_tbl + p_for) + plot_layout(widths = c(2, 1.52))
dev.off()




library(dplyr)
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
data_raw_1 <- read.csv('./animalTFDB4_Mus_musculus_TF.txt',sep = '\t')
data_raw_2 <- read.csv('./trrust_rawdata.mouse.tsv',sep = '\t',col.names = c(
  'tf','tg','impact','protein_id'
))
res_tf <- data.frame(
  tf_name = c(
    data_raw_1$Symbol,data_raw_2$tf
  ) %>% unique()
)


data_plot <- list(
  'animalTFDB4' = data_raw_1$Symbol %>% unique(),
  'trrust' = data_raw_2$tf %>% unique()
  )

library(eulerr)
plot(
  euler(data_plot),
  labels = list(col="black",font=3,cex=2.5),
  edges = list(col="darkgreen",lwd=5,
               lty=1:3),
  quantities = list(cex=2)
)


res_tf$tf_name %>% unique() %>% length()
setwd('./')




data_pro_clean <- read.csv('./data_copy_num_sorted.csv',row.names = 1) %>% 
  dplyr::select(-MolecularWeight) %>% 
  rename(Accessions = ProteinGroups)
sum((data_pro_clean$Genes %>% nchar())==0)
data_gene_protein <- data_pro_clean[,c('Genes','Accessions')]
data_pro_clean <- data_pro_clean %>% 
  filter(nchar(Genes)>0)
protein_list <- data_pro_clean %>% 
  dplyr::select(Genes,Accessions) %>% 
  rename(protein = Accessions)

scdata_pro <- readRDS('./scdata_use_filteroutliers.rds')
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
  column_to_rownames(var = 'Genes')
scdata_pro <- CreateSeuratObject(
  counts = data_singlestage %>% as.matrix(),
  meta.data = scdata_pro@meta.data
)

scdata_tran <- readRDS('./scdata_tpm_filteroutliters.rds')
data_tran <- scdata_tran@assays$RNA$counts %>% as.data.frame()
meta_data_tran <- scdata_tran@meta.data %>% 
  mutate(group_stage_cell = paste(group_stage,group_cell,sep = ': ')) %>% 
  mutate(group_stage_cell = factor(group_stage_cell,levels = group_stage_cell %>% unique())) %>% 
  mutate(group_omics = 'transcription')
gene_select <- rownames(scdata_tran)[grepl('H2',rownames(scdata_tran))]
data_tmp <- FetchData(scdata_tran,vars = c(gene_select,'sample'))
scdata_use <- subset(scdata_tran,features = res_tf$tf_name)
scdata_use@assays$RNA$counts %>% as.data.frame() %>% 
  filter(rowSums(.)>0) %>% nrow()
data_plot <- list(
  exp_pro = data_pro %>% 
    rownames_to_column('Genes') %>% 
    filter(Genes %in% res_tf$tf_name) %>% 
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes),
  exp_tran = data_tran %>% 
    rownames_to_column('Genes') %>% 
    filter(Genes %in% res_tf$tf_name) %>% 
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes)
)
library(eulerr)
library(ggplot2)

pdf('./venn_pro-trans_20250515.pdf',width = 9,height = 6)
par(mar = c(30, 5, 5, 5))
plot(
  euler(data_plot),
  labels = list(col="black",font=3,cex=2.5),
  edges = list(col=my36colors[1:4],lwd=5,lty=1),
  fills = list(fill = 'white',font = 3,cex=2),
  quantities = list(cex=2)
)
dev.off()
data_plot[[1]][!(data_plot[[1]] %in% data_plot[[2]])]


data_plot <- list(
  TF_Oo_Trans = subset(scdata_tran,features = res_tf$tf_name,subset = group_cell == 'Oocyte')@assays$RNA$counts %>% 
    as.data.frame() %>% 
    rownames_to_column('Genes') %>%
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes),
  TF_GC_Trans = subset(scdata_tran,features = res_tf$tf_name,subset = group_cell == 'Granulosa')@assays$RNA$counts %>% as.data.frame() %>% 
    rownames_to_column('Genes') %>%
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes),
  TF_Oo_Pro = subset(scdata_pro,features = res_tf$tf_name,subset = group_cell == 'Oocyte')@assays$RNA$counts %>% as.data.frame() %>%
    rownames_to_column('Genes') %>%
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes),
  TF_GC_Pro = subset(scdata_pro,features = res_tf$tf_name,subset = group_cell == 'Granulosa')@assays$RNA$counts %>% as.data.frame() %>%
    rownames_to_column('Genes') %>%
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes)
)
plot(
  euler(data_plot),
  labels = list(col="black",font=3,cex=2.5),
  edges = list(col=my36colors[1:4],lwd=5,lty=1),
  fills = list(fill = 'white',font = 3,cex=2),
  quantities = list(cex=2)
)
library(UpSetR)
pdf('./upset.pdf',width = 12,height = 5.6)
upset(data = fromList(data_plot),nintersects = 100)
dev.off()


upset(data = fromList(data_plot),      
      

      

      

      
      number.angles = 0,
      
      point.size=4,
      
      line.size=1,
      
      mainbar.y.label="Intersection size",
      
      main.bar.color = 'black',
      
      matrix.color="black",
      
      sets.x.label="Set size",
      

      
      mb.ratio = c(0.7, 0.3),
      
      order.by = "freq",
      
      decreasing = c(T,F),
      
      text.scale=c(1.5,1.5,1.5,1.5,1.5,1),
      
      shade.color="red",
      

      
      queries=list(

        list(query=intersects,params=list("TF_Oo_Trans","TF_GC_Trans"),color="yellow",active=T),
        list(query=intersects,params=list("TF_Oo_Pro","TF_GC_Pro","TF_Oo_Trans","TF_GC_Trans"),color="purple",active=T)
      )
      
)

max_length <- max(sapply(data_plot, length))
data_plot_filled <- lapply(data_plot, function(x) {
  length(x) <- max_length
  return(x)
})
df_res <- data_plot_filled %>% as.data.frame()
write.csv(
  df_res,
  './suppl_table_tf.csv',
  quote = FALSE,row.names = FALSE)

df <- as.data.frame(my_list_filled)


scdata_use <- subset(scdata_pro,features = res_tf$tf_name)
scdata_use <- subset(
  scdata_use,subset = group_cell == 'Oocyte'
)



scdata_use<-FindVariableFeatures(scdata_use,nfeatures = nrow(scdata_use))
high_variable <- VariableFeatures(scdata_use)

library(tidyverse)
library(broom)
library(palmerpenguins)
library(BiocSingular)
library(ggplot2)
library(ggrepel)
library(stats)
data_use <- scdata_use@assays$RNA$counts %>% as.data.frame()
PCA.sample <-
  data.frame(t(data_use[rownames(data_use) %in% high_variable, ]))

PCA.sample<-PCA.sample[ , which(apply(PCA.sample, 2, var) != 0)]
pca.result<-prcomp(PCA.sample)
pca<-data.frame(pca.result$x)
pca.label<-cbind(sample=rownames(pca),pca) %>% 
  left_join(scdata_use@meta.data %>% dplyr::select(sample,group_cell,group_stage_cell),by = 'sample') %>% 
  mutate(
    group_cell = factor(group_cell,levels = c(
      'Oocyte','Granulosa'
    )),
    group_stage_cell = factor(group_stage_cell,levels = c(
      "Secondary: Oocyte","Early antral: Oocyte","Antral: Oocyte","Preovulatory: Oocyte",
      "Secondary: Granulosa","Early antral: Granulosa","Antral: Granulosa","Preovulatory: Granulosa"
    ))
  )

summ = summary(pca.result)
xlab = paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab = paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

p_pca <- ggplot(data = pca.label,aes(x = PC1,y = PC2,color = group_stage_cell))+ 
  geom_point(size = 8) +
  guides(fill = "none") +
  scale_color_manual(values = c(
    '#FFAD72','#F76D5E','#D82632','#A50021'
  )) +
  labs(x=xlab,y=ylab,color = "Groups")+
  theme_bw()+
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 24),
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm')
  )
p_pca
scdata_use <- NormalizeData(scdata_use)
scdata_use <- ScaleData(scdata_use,features = rownames(scdata_use))
scdata_use <- RunPCA(
  scdata_use,
  features = VariableFeatures(object = scdata_use),
  npcs = 30)
ElbowPlot(scdata_use,ndims = 30)
pc.num = 1:20
scdata_use <- RunUMAP(
  scdata_use,n.neighbors = 10,
  dims = pc.num)
DimPlot(scdata_use,group.by = 'group_stage_cell',pt.size = 3)
data_plot <-
  DimPlot(
    object = scdata_use,
    pt.size = 1,
    reduction = 'umap',
    group.by = 'group_stage'
  )$data %>%
  rownames_to_column('sample')
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage,
         label = sample
       )) +
  geom_point(size = 6) +
  scale_color_manual(values = c('#FFAD72','#F76D5E','#D82632','#A50021')) +
  guides(color = guide_legend(
    title = 'Stage',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )
ggsave(
  './umap_tf_pro.pdf',
  width = 9,height = 6
)
ggsave(
  './umap_tf_pro.png',
  width = 9,height = 6
)


scdata_use <- subset(scdata_tran,features = res_tf$tf_name)
scdata_use <- subset(
  scdata_use,subset = group_cell == 'Oocyte'
)
tmp_oocyte <- scdata_use@assays$RNA$counts %>% as.data.frame() %>% filter(rowSums(.)>0)



number_high <- rownames(scdata_use) %>% length()
scdata_use<-FindVariableFeatures(scdata_use,nfeatures = 800)
high_variable <- VariableFeatures(scdata_use)

library(tidyverse)
library(broom)
library(palmerpenguins)
library(BiocSingular)
library(ggplot2)
library(ggrepel)
library(stats)
data_use <- scdata_use@assays$RNA$counts %>% as.data.frame()
PCA.sample <-
  data.frame(t(data_use[rownames(data_use) %in% high_variable, ]))

PCA.sample<-PCA.sample[ , which(apply(PCA.sample, 2, var) != 0)]
pca.result<-prcomp(PCA.sample)
pca<-data.frame(pca.result$x)
pca.label<-cbind(sample=rownames(pca),pca) %>% 
  left_join(scdata_use@meta.data %>% dplyr::select(sample,group_cell,group_stage_cell),by = 'sample') %>% 
  mutate(
    group_cell = factor(group_cell,levels = c(
      'Oocyte','Granulosa'
    )),
  )

summ = summary(pca.result)
xlab = paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab = paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

p_pca <- ggplot(data = pca.label,aes(x = PC1,y = PC2,color = group_stage_cell))+ 
  geom_point(size = 8) +
  guides(fill = "none") +
  scale_color_manual(values = c(
    '#FFAD72','#F76D5E','#D82632','#A50021'
  )) +
  labs(x=xlab,y=ylab,color = "Groups")+
  theme_bw()+
  theme(
    plot.title = element_blank(),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 24),
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm')
  )
p_pca



scdata_use <- NormalizeData(scdata_use)
scdata_use <- ScaleData(scdata_use,features = rownames(scdata_use))
scdata_use<-FindVariableFeatures(scdata_use,nfeatures = 500)
high_variable <- VariableFeatures(scdata_use)
scdata_use <- RunPCA(
  scdata_use,
  features = VariableFeatures(object = scdata_use),
  npcs = 30)
ElbowPlot(scdata_use,ndims = 30)
pc.num = 1:30
scdata_use <- RunUMAP(
  scdata_use,n.neighbors = 18,
  dims = pc.num)
DimPlot(scdata_use,group.by = 'group_stage_cell',pt.size = 3)
data_plot <-
  DimPlot(
    object = scdata_use,
    pt.size = 1,
    reduction = 'umap',
    group.by = 'group_stage'
  )$data %>%
  rownames_to_column('sample')
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage,
         label = sample
       )) +
  geom_point(size = 6) +
  scale_color_manual(values = c('#FFAD72','#F76D5E','#D82632','#A50021')) +
  guides(color = guide_legend(
    title = 'Stage',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )
ggsave(
  './umap_tf_trans.pdf',
  width = 9,height = 6
  )
ggsave(
  './umap_tf_trans.png',
  width = 9,height = 6
)


scdata_use <- subset(scdata_pro,features = res_tf$tf_name)
scdata_use <- subset(
  scdata_use,subset = group_cell == 'Granulosa'
)
tmp <- scdata_use@assays$RNA$counts %>% as.data.frame() %>% filter(rowSums(.)>0)
scdata_use <- NormalizeData(scdata_use)
scdata_use <- ScaleData(scdata_use,features = rownames(scdata_use))
scdata_use<-FindVariableFeatures(scdata_use,nfeatures = nrow(scdata_use))
scdata_use <- RunPCA(
  scdata_use,
  features = VariableFeatures(object = scdata_use),
  npcs = 30)
ElbowPlot(scdata_use,ndims = 30)
pc.num = 1:30
scdata_use <- RunUMAP(
  scdata_use,n.neighbors = 8,
  dims = pc.num)
DimPlot(scdata_use,group.by = 'group_stage_cell',pt.size = 3)
data_plot <-
  DimPlot(
    object = scdata_use,
    pt.size = 1,
    reduction = 'umap',
    group.by = 'group_stage'
  )$data %>%
  rownames_to_column('sample')
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage,
         label = sample
       )) +
  geom_point(size = 6) +
  scale_color_manual(values = c('#B9DDF1FF','#7EAED3FF','#5081AEFF','#2A5783FF')) +
  guides(color = guide_legend(
    title = 'Stage',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )

ggsave(
  './umap_tf_pro_Granulosa.pdf',
  width = 9,height = 6
)
ggsave(
  './umap_tf_pro_Granulosa.png',
  width = 9,height = 6
)


scdata_use <- subset(scdata_tran,features = res_tf$tf_name)
scdata_use <- subset(
  scdata_use,subset = group_cell == 'Granulosa'
)
scdata_use <- NormalizeData(scdata_use)
scdata_use <- ScaleData(scdata_use,features = rownames(scdata_use))
scdata_use<-FindVariableFeatures(scdata_use,nfeatures = 500)
high_variable <- VariableFeatures(scdata_use)
scdata_use <- RunPCA(
  scdata_use,
  features = VariableFeatures(object = scdata_use),
  npcs = 30)
ElbowPlot(scdata_use,ndims = 30)
pc.num = 1:30
scdata_use <- RunUMAP(
  scdata_use,n.neighbors = 8,
  dims = pc.num)
DimPlot(scdata_use,group.by = 'group_stage_cell',pt.size = 3)
data_plot <-
  DimPlot(
    object = scdata_use,
    pt.size = 1,
    reduction = 'umap',
    group.by = 'group_stage'
  )$data %>%
  rownames_to_column('sample')
ggplot(data_plot,
       mapping = aes(
         x = umap_1,
         y = umap_2,
         color = group_stage,
         label = sample
       )) +
  geom_point(size = 6) +
  scale_color_manual(values = c('#B9DDF1FF','#7EAED3FF','#5081AEFF','#2A5783FF')) +
  guides(color = guide_legend(
    title = 'Stage',
    title.theme = element_text(size = 20),
    override.aes = list(size = 4)
  )) +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 20
    ),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 18)
  )
ggsave(
  './umap_tf_trans_Granulosa.pdf',
  width = 9,height = 6
)
ggsave(
  './umap_tf_trans_Granulosa.png',
  width = 9,height = 6
)



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

data_tran_regulon_Oocyte <- read.csv('./data_standardized.csv',row.names = 1)
data_tran_regulon_GC <- read.csv('./data_standardized.csv',row.names = 1)
data_tf_use <- data_tf_all %>% 
  filter(
    TF %in% unique(rownames(data_tran_regulon_Oocyte),rownames(data_tran_regulon_GC)),
    Group %in% 'trrust'
  )
data_tf_use$Group %>% unique()


data_plot <- list(
  exp_pro = data_pro %>% 
    rownames_to_column('Genes') %>% 
    filter(Genes %in% res_tf$tf_name) %>% 
    rowwise() %>%
    summarise(
      positive_count = sum(c_across(-Genes) > 0),
      at_least_two_positive = positive_count >= 3,
      Genes = Genes
    ) %>% 
    filter(at_least_two_positive) %>% 
    pull(Genes),
  exp_Regulon = data_tf_use$TF
)
library(eulerr)
library(ggplot2)
pdf('./venn_pro-regulon.pdf',width = 9,height = 6)
par(mar = c(30, 5, 5, 5))
plot(
  euler(data_plot),
  labels = list(col="black",font=3,cex=2.5),
  edges = list(col=my36colors[1:4],lwd=5,lty=1),
  fills = list(fill = 'white',font = 3,cex=2),
  quantities = list(cex=2)
)
dev.off()

library(Seurat)
library(Biobase)
library(dplyr)
library(viper)
setwd('./data_tf/')

scdata <- readRDS('./scdata_tpm_filteroutliters.rds')
scdata@assays[["RNA"]]@layers[["data_tpm"]] <- NULL
scdata@assays[["RNA"]]@layers[["data_fpkm"]] <- NULL
scdata <- subset(scdata,subset = group_cell == 'Oocyte')
scdata$group_stage %>% unique()
group_stage_use <- 'Preovulatory'
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
tmp %>% 
  filter(Regulon %in% c('Smarca4','Gata4'))


par(mai = c(0,0,0.5,1))
tmp <- data.frame(
  nes = mrs$es$nes,
  p.value = mrs$es$p.value,
  activity = (mrs$es$nes/max(abs(mrs$es$nes)))
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
mrs_select <- c('Smarca4','Gata4')
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
