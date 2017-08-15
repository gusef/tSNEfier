setwd("C:/Users/Daniel Gusenleitner/Dropbox (Partners HealthCare)/github_repos/tSNEfier/example")
rm(list = ls ())
gc()
require(Biobase)

#load the gene expression data
counts_cpm <- readRDS('c:/work/projects/melanoma_datasets/processed_datasets/RNAseq_nodedup_cpm.RDS')

#filter out samples that do not pass IHC QC - after tr
counts_cpm <- counts_cpm[,is.na(pData(counts_cpm)$HnE.Quality) | pData(counts_cpm)$HnE.Quality != 'F']

#drop the NE sample in BORI
counts_cpm$BORI[counts_cpm$BORI=='NE'] <- NA
counts_cpm$BORI <- droplevels(counts_cpm$BORI)
counts_cpm <- counts_cpm[,!is.na(counts_cpm$BORI)]
counts_cpm$notPD <- as.factor(c('notPD','PD')[(counts_cpm$BORI=='PD')+1])

#bring the data into a log2 space and drop empty genes
counts_cpm <- counts_cpm[rowSums(exprs(counts_cpm))>1,]
exprs(counts_cpm) <- log2(exprs(counts_cpm)+1)

#get only the genes that have a gene symbol and remove all duplicates
genes <- fData(counts_cpm)$hgnc_symbol
counts_cpm <- counts_cpm[!is.na(genes) & !duplicated(genes),]
featureNames(counts_cpm) <- fData(counts_cpm)$hgnc_symbol

#remove all rnaseq data that was not pre treatment
counts_cpm <- counts_cpm[,counts_cpm$Visit.Code == 'SCREEN']
saveRDS(counts_cpm,file='gene_expression_data.rds')

#make one comprehensive .gmt file
file_dir <- 'gene_sets'
files <- dir(file_dir)
gene_sets <- lapply(files,function(x)readLines(file.path(file_dir,x)))
gene_sets <- unlist(gene_sets)
names(gene_sets) <- sapply(strsplit(gene_sets,'\t'),function(x)x[1])
gene_sets <- gene_sets[!duplicated(names(gene_sets))]

#add the Nanostring gene immune gene list
genes <- as.character(read.csv('C:/work/projects/2017_03_064_figures_and_tables/nivo_biomarker/Nanostring_genelist.txt')[,2])
genes <- sort(unique(genes))
gene_sets <- c(paste(c('Nanostring_Immuneset','',genes),collapse = '\t'),gene_sets)

writeLines(gene_sets,con = file('gene_sets.gmt'))



