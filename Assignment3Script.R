if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")
if(!requireNamespace("edgeR", quietly = TRUE))
  install.packages("edgeR")
if(!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")
if(!requireNamespace("GEOquery", quietly = TRUE))
  install.packages("GEOquery")
if(!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")
BiocManager::install("RCy3")

library(GEOquery)
library(edgeR)
library(biomaRt)
library(tidyr)
library(RCy3)

#The following code is copied from previous assignments done by me(Emily Ayala)
#for the course BCB420.  This code was also copied from lecture slides written by Ruth Isserlin
#Load the data
sfiles = getGEOSuppFiles('GSE98922')

fnames = rownames(sfiles)

#Read the file
exp = read.delim(fnames[1],header = TRUE,check.names = FALSE)

#Get rid of low counts
cpms <- edgeR:: cpm(exp[,2:10])
rownames(cpms) <- exp[,1]
keep = rowSums(cpms > 1) >= 3
exp_filtered = exp[keep,]
rownames(exp_filtered) <- rownames(cpms[keep,])

filtered_data_matrix <- as.matrix(exp_filtered[,2:10])
rownames(filtered_data_matrix) <- rownames(exp_filtered)

#Separate into groups
samples <- data.frame(c('control', '1'), c('control', '2'), c('control', '3'), c('treatmentA', '1'), c('treatmentA', '2'), c('treatmentA', '3'), c('treatmentD', '1'), c('treatmentD', '2'), c('treatmentD', '3'))
rownames(samples) <- c("treatment", "trial")
colnames(samples) <- colnames(exp_filtered)[2:10]
samples <- data.frame(t(samples))

d <- DGEList(counts = filtered_data_matrix, group = samples$treatment)

#Create model
model_design_pat <- model.matrix(~ samples$treatment + samples$trial)

#Estimate dispersion
d <- estimateDisp(d, model_design_pat)

#Calculate normalization factors
d <- calcNormFactors(d)

#fit model
fit <- glmQLFit(d, model_design_pat)

#calculate differential expression
qlf.pos_vs_neg <- glmQLFTest(fit)

qlf_output_hits <- topTags(qlf.pos_vs_neg, sort.by = "PValue", n = nrow(filtered_data_matrix))

#Create thresholded lists of gene
qlf_output_hits_withgn <- merge(exp[,1:2], qlf_output_hits, by.x = 1, by.y = 0)
qlf_output_hits_withgn[,"rank"] <- - log(qlf_output_hits_withgn$PValue, base = 10) * sign(qlf_output_hits_withgn$logFC)
qlf_output_hits_withgn <- qlf_output_hits_withgn[order(qlf_output_hits_withgn$rank),]

#Make a ranked genelist
ranked_list <- data.frame(symbol = qlf_output_hits_withgn$Var.1, rank = qlf_output_hits_withgn$rank)

#Write that ranked list to a file
write.table(x=ranked_list, file=file.path("data", "ranked_gene_list.rnk"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#The following is copied from Enrichment Map Analysis Pipeline by Ruth Isserlin 2019-11-18 
#https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/Protocol2_createEM.html#24_run_gsea
#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  install.packages("RCurl")}, 
         finally = library("RCurl"))

#use library
tryCatch(expr = { library("limma")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("limma")}, 
         finally = library("limma"))
tryCatch(expr = { library("Biobase")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("Biobase")}, 
         finally = library("Biobase"))
tryCatch(expr = { library("ggplot2")}, 
         error = function(e) { install.packages("ggplot2")}, 
         finally = library("ggplot2"))

#For creating json and communicating with cytoscape
tryCatch(expr = { library("httr")}, 
         error = function(e) { install.packages("httr")}, 
         finally = library("httr"))
tryCatch(expr = { library("RJSONIO")}, 
         error = function(e) { install.packages("RJSONIO")}, 
         finally = library("RJSONIO"))

#Configurable Parameters
parameters <- read.table("parameters.txt", header = FALSE, sep = "", dec = ".", stringsAsFactors = FALSE)
parameters[1, 1]
gsea_jar <- parameters$V1[1] 
java_version <- 11
working_dir <- parameters$V1[2]

analysis_name <- "Assignment 3"
rnk_file <- "ranked_gene_list.rnk"

#Download the latest pathway definition file
gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"

filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(working_dir,paste("Supplementary_Table3_",gmt_file,sep="") )

download.file(
  paste(gmt_url,gmt_file,sep=""),
  destfile=dest_gmt_file
)

#Create a command to run
command <- paste("", gsea_jar, "GSEAPreRanked -gmx", dest_gmt_file, "-rnk", file.path(working_dir, rnk_file), "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ", analysis_name, "-plot_top_x 20 -rnd_seed 54321 -set_max 200 -set_min 15 -zip_report false -out" ,working_dir, "> gsea_output", sep = " ")
#system2(gsea_jar, args = c("GSEAPreRanked", "-gmx", dest_gmt_file, "-rnk", file.path(working_dir, rnk_file), "-collapse", "false", "-nperm", "1000", "-scoring_scheme", "weighted", "-rpt_label", analysis_name, "-plot_top_x", "20", "-rnd_seet", "54321", "-set_max", "200", "-set_min", "15", "-zip_report", "false", "-out", working_dir), stdout = "C:/Users/Emily/Desktop/Desktop/Y5S2/BCB420/Assignment3BCB420/data/gsea_output.txt", stderr = "")
system(command)

gsea_directories <- list.files(path = working_dir, pattern = "\\.GseaPreranked")

details = file.info(file.path(getwd(), working_dir, gsea_directories))

details = details[with(details, order(as.POSIXct(mtime), decreaseing = TRUE)), ]

gsea_output_dir <- row.names(details)[1]

#Set up connection from R to cytoscape
port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

version.url = paste(base.url, "version", sep="/")
cytoscape.open = TRUE

tryCatch(expr = { GET(version.url)},
         error = function(e) {return (cytoscape.open = FALSE)}, finally = function(r){ return(cytoscape.open = TRUE)})
if(!cytoscape.open){
  #try and launch cytoscape
  print("Cytoscape is not open.  Please launch cytoscape.")
} else{
  cytoscape.version =  GET(version.url)
  cy.version = fromJSON(rawToChar(cytoscape.version$content))
  print(cy.version)
  
}

#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library("RCy3")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("RCy3")}, finally = library("RCy3"))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_gsea_threshold <- "0.05"
qvalue_gsea_threshold <- "0.01"

similarity_threshold <- "0.375"
similarity_metric = "COMBINED"

cur_model_name <- analysis_name

gsea_results_path <- file.path(gsea_output_dir,"edb")
gsea_results_filename <- file.path(gsea_results_path,"results.edb")

#although there is a gmt file in the gsea edb results directory it have been filtered to 
#contain only genes represented in the expression set.  If you use this fltered file you 
#will get different pathway connectivity depending on the dataset being used.  We recommend 
#using original gmt file used for the gsea analysis and not the filtered one in the results directory.
gmt_gsea_file <- file.path(getwd(),dest_gmt_file)
gsea_ranks_file <- file.path(gsea_results_path,list.files(gsea_results_path,pattern=".rnk"))

#######################################
#create EM
current_network_name <- paste(cur_model_name,pvalue_gsea_threshold,qvalue_gsea_threshold,sep="_")

em_command = paste('enrichmentmap build analysisType="gsea" gmtFile=',gmt_gsea_file,
                   'pvalue=',pvalue_gsea_threshold, 'qvalue=',qvalue_gsea_threshold,
                   'similaritycutoff=',similarity_threshold,
                   'coefficients=',similarity_metric,'ranksDataset1=', 
                   gsea_ranks_file,'enrichmentsDataset1=',gsea_results_filename, 
                   'filterByExpressions=false',
                   'expressionDataset1=',file.path(getwd(),working_dir,expression_file),
                   sep=" ")

#enrichment map command will return the suid of newly created network.
response <- commandsGET(em_command)

current_network_suid <- 0
#enrichment map command will return the suid of newly created network unless it Failed.  
#If it failed it will contain the word failed
if(grepl(pattern="Failed", response)){
  paste(response)
} else {
  current_network_suid <- response
}

#check to see if the network name is unique
current_names <- getNetworkList()
if(current_network_name %in% current_names){
  #if the name already exists in the network names then put the SUID in front
  # of the name (this does not work if you put the suid at the end of the name)
  current_network_name <- paste(current_network_suid,current_network_name,  sep="_")
}
response <- renameNetwork(title=current_network_name, 
                          network = as.numeric(current_network_suid),base.url)