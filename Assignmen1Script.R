#Install required packages
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("GEOmetadb", quietly = TRUE))
  BiocManager::install("GEOmetadb")
if(!requireNamespace("kableExtra", quietly = TRUE))
  install.packages("kableExtra")
if(!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
if(!requireNamespace("biomaRt"), quietly = TRUE)
  BiocManager::install("biomaRt")
library(BiocManager)
library(GEOmetadb)
library(knitr)
library(kableExtra)
library(edgeR)
library(biomaRt)
#Get the meta data get SQLiteFile()
if(!file.exists('GEOmetadb.sqlite'))
  getSQLiteFile()

#Connect to the meta data database
con <- dbConnect(SQLite(), 'GEOmetadb.sqlite')

#Get the tables
geo_tables <- dbListTables(con)
geo_tables
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
unlist(lapply(rs$supplementary_file,
              FUN = function(x){x <- unlist(strsplit(x, ";")) ;
              x <- x[grep(x,pattern="txt",ignore.case = TRUE)];
              tail(unlist(strsplit(x,"/")),n=1)}))[1:10]
counts_files <- rs$supplementary_file[grep(rs$supplementary_file, pattern = "count", ignore.case = TRUE)]

#Get the files for GSE81475
sfiles = getGEOSuppFiles('GSE81475')
fnames = rownames(sfiles)
b2 = read.delim(fnames[1],header=TRUE)
head(b2)

#Get the GEO description
gse <- getGEO("GSE81475", GSEMatrix=FALSE)
kable(data.frame(head(Meta(gse))), format = "html")
current_gpl <- names(GPLList(gse))[1]
current_gpl_info <- Meta(getGEO(current_gpl))
current_gpl_info$title
current_gpl_info$last_update_date
current_gpl_info$organism

exp = read.delim(fnames[2],header = TRUE,check.names = FALSE)
exp
kable(exp[1:25,1:10], format = "html")
dim(exp)
colnames(exp)
rownames(exp)
samples <-data.frame(t(exp))
table(exp$Geneid)
#Get the gene counts
summarized_gene_counts <- sort(table(exp$Geneid), decreasing = TRUE)
length(summarized_gene_counts[which(summarized_gene_counts == 1)])
kable(summarized_gene_counts[which(summarized_gene_counts > 1)], format = "html")

#Translate out counts into counts per million
cpms = cpm(exp[,2:1609])
rownames(cpms) <- exp[,1]
rownames(cpms)
keep = rowSums(cpms > 1) >= 3
exp_filtered = exp[keep,]
length(exp_filtered)
dim(exp_filtered)

#Distribution of our data
par("mar")
par(mar=c(1,1,1,1))
data2plot <- log2(cpm(exp_filtered[,2:1609]))
#The following line crashes my computer.  
#boxplot(data2plot, xlab = "Samples", ylab = "log2 CPM", las = 2, cex = 0.5, cex.lab = 0.5, cex.axis = 0.5, main = "	Zika Virus Disrupts Phospho-TBK1 Localization and Mitosis in Human Neural Stem Cell Model Systems")
filtered_data_matrix <- as.matrix(exp_filtered[,2:1609])
rownames(samples)
d = DGEList(counts=filtered_data_matrix)
d = calcNormFactors(d)
normalized_counts <- cpm(d)

#Identifier Mapping
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

exp_filtered$Geneid
typeof(exp_filtered$Geneid)
exp_filtered[1]
exp_filtered_gene_id <- data.frame(apply(exp_filtered$Geneid, FUN=function(x){x <- unlist(strsplit(as.character(x), split = '\\|'))}))
exp_filtered_gene_id[1]
conversion_stash <- "id_conversion.rds"
if(file.exists(conversion_stash)){
  id_conversion <- readRDS(conversion_stash)
} else{
  id_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = c("ensemble_gene_id"),
                         values = exp_filtered$ensemble_id,
                         mart = ensembl)
  saveRDS(id_conversion, conversation_stash)
}


