#Install required packages
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("GEOmetadb", quietly = TRUE))
  BiocManager::install("GEOmetadb")
if(!requireNamespace("kableExtra", quietly = TRUE))
  install.packages("kableExtra")
if(!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
if(!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if(!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("plyr")
library(BiocManager)
library(GEOmetadb)
library(knitr)
library(edgeR)
library(biomaRt)
library(plyr)
#Get the meta data get SQLiteFile()
if(!file.exists('GEOmetadb.sqlite'))
  getSQLiteFile()

#Connect to the meta data database
con <- dbConnect(SQLite(), 'GEOmetadb.sqlite')

#Get the tables
geo_tables <- dbListTables(con)



#Query the tables for RNASeq data sets
sql <- paste("SELECT DISTINCT gse.title, gse.gse, gpl.title, gse.submission_date, gse.supplementary_file",
             "FROM",
             " gse JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             " JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             " gse.submission_date > '2015-01-01' AND",
             " gse.title LIKE '%MS%' AND",
             " gpl.title LIKE '%HiSeq%' AND",
             " gpl.organism = 'Homo sapiens'",
             sep=" ")
rs <- dbGetQuery(con, sql)

#Disconnect from the database
dbDisconnect(con)

unlist(lapply(rs$supplementary_file,
              FUN = function(x){x <- unlist(strsplit(x, ";")) ;
              x <- x[grep(x,pattern="txt",ignore.case = TRUE)];
              tail(unlist(strsplit(x,"/")),n=1)}))[1:10]
counts_files <- rs$supplementary_file[grep(rs$supplementary_file, pattern = "count", ignore.case = TRUE)]

#Get the files for GSE81475
sfiles = getGEOSuppFiles('GSE81475')
fnames = rownames(sfiles)
b2 = read.delim(fnames[2],header=TRUE)

#Get the GEO description
gse <- getGEO("GSE81475", GSEMatrix=FALSE)
current_gpl <- names(GPLList(gse))[1]
current_gpl_info <- Meta(getGEO(current_gpl))


exp = read.delim(fnames[2],header = TRUE,check.names = FALSE)



samples <-data.frame(t(exp))

#Get the gene counts
summarized_gene_counts <- sort(table(exp$Geneid), decreasing = TRUE)
length(summarized_gene_counts[which(summarized_gene_counts == 1)])
kable(summarized_gene_counts[which(summarized_gene_counts > 1)], format = "html")

dim(exp)
#Translate out counts into counts per million
cpms = cpm(exp[,2:1609])
rownames(cpms) <- exp[,1]
keep = rowSums(cpms > 1) >= 3

exp_filtered = exp[keep,]

exp_filtered$Geneid

#Distribution of our data
par("mar")
par(mar=c(1,1,1,1))
data2plot <- log2(cpm(exp_filtered[,2:1609]))
#The following line crashes my computer.  
#boxplot(data2plot, xlab = "Samples", ylab = "log2 CPM", las = 2, cex = 0.5, cex.lab = 0.5, cex.axis = 0.5, main = "	Zika Virus Disrupts Phospho-TBK1 Localization and Mitosis in Human Neural Stem Cell Model Systems")
filtered_data_matrix <- as.matrix(exp_filtered[,2:1609])



#Identifier Mapping
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

exp_filtered$Geneid
typeof(exp_filtered$Geneid)
exp_filtered[1]
exp_filtered_gene_id <- list()
exp_filtered_gene_id$ensembl_gene <- c(unlist(lapply(exp_filtered$Geneid, FUN=function(x){x <- unlist(strsplit(as.character(x), split = '\\|'))[[1]]})))
typeof(exp_filtered_gene_id$ensembl_gene[1])
rownames(filtered_data_matrix) <- exp_filtered_gene_id$ensembl_gene
#Get normalized counts
d = DGEList(counts=filtered_data_matrix)
d = calcNormFactors(d)
normalized_counts <- cpm(d)

conversion_stash <- "id_conversion.rds"
if(file.exists(conversion_stash)){
  id_conversion <- readRDS(conversion_stash)
} else{
  id_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = c("ensembl_gene_id"),
                         values = unlist(exp_filtered_gene_id$ensembl_gene),
                         mart = ensembl)
  saveRDS(id_conversion, conversion_stash)
}
nrow(normalized_counts) - nrow(id_conversion)
typeof(normalized_counts)
normalized_counts_annot <- merge(id_conversion, normalized_counts, by.x = 1, by.y = 0, all.y = TRUE)
kable(normalized_counts_annot[1:5,1:5],type = "html")
ensembl_id_missing_gene <- normalized_counts_annot$ensembl_gene_id[which(is.na(normalized_counts_annot$hgnc_symbol))]
length(ensembl_id_missing_gene)
typeof(normalized_counts_annot)
exp_dataframe <- data.frame(normalized_counts_annot)
