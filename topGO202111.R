#!/usr/bin/env Rscript

# Author: Xiaodong Li
# Date:   Nov. 29th, 2021
# R version: R 4.0.2 (2020-06-22)

if (!requireNamespace("topGO", quietly = TRUE))
    BiocManager::install("topGO")
if (!requireNamespace("Rgraphviz", quietly = TRUE))
    BiocManager::install("Rgraphviz")

if (!requireNamespace("getopt", quietly = TRUE)){
    install.packages("getopt");
    q()
}

# library packages
suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))
suppressMessages(library(getopt))

#############################
# getting parameters
#############################
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
                'help'          , 'h', 0, "logical"  , "for help",
                'input'         , 'i', 1, "character", "input the gene names file [required, csv format, only one column, no header]",
                'outdir'        , 'o', 1, "character", "output file outdir [optional, default: cwd]",
                'map'           , 'm', 1, "character", "a GO-to-genes mappings file separated by tab delimiter [required]",
                'name'          , 'n', 1, "character", "out file name prefix [optional, default: output]",
                'firstSigNodes' , 'N', 1, "integer"  , "the number of top scoring GO terms which .... [optional, default: 5]",
                'ontology'      , 'O', 1, "character", "character string specifying the ontology of interest (BP, MF or CC) [optional, default: BP]"
        ), byrow=TRUE, ncol=5);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
    cat("")
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$input) )    { print_usage(spec) }
if ( is.null(opt$outdir) )    { opt$outdir=getwd()}
if ( is.null(opt$map) )    { print_usage(spec) }
if ( is.null(opt$name) )    { opt$name="output" }
if ( is.null(opt$firstSigNodes) )    { opt$firstSigNodes=5 }
if ( is.null(opt$ontology) )    { opt$ontology="BP" }

#set some reasonable defaults for the options that are needed,
#but were not specified.
if( !file.exists(opt$outdir) ){
    if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
        stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
    }
}


#############################
# Data preparation
#############################
# 1 A mapping between gene identifiers and GO terms
## Annotations need to be provided either as gene-to-GOs or as GO-to-genes mappings.
## Gene-to-GOs text format: gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3, ....
geneID2GO <- readMappings(file = opt$map)
#str(head(geneID2GO))
geneNames <- names(geneID2GO)
#head(geneNames)
#length(geneNames)

# 2 A list of gene identifiers and optionally the gene-wise scores
genes <- read.csv(opt$input, header = FALSE, colClasses = c("character"))
myInterestingGenes <- genes[, 1]
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
#str(geneList)

# 3 The GO hierarchical structure.This structure is obtained from the GO.db package.

# create a topGOdata object
GOdata <- new("topGOdata", ontology = opt$ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
print(GOdata)

# all GO terms are available for analysis
ug <- usedGO(GOdata)

#############################
# Running the enrichment tests
#############################

# runTest functions return an object of type topGOresult
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
pvalFis <- score(resultFis)

#############################
# Save the results
#############################

# set your output data directory
setwd(opt$outdir)

# Summarising the results
allRes <- GenTable(GOdata, classic = resultFis, orderBy = "classic", ranksOf = "classic", topNodes = length(ug))
write.csv(allRes, paste(opt$name, "topGO.csv", sep = "_"), row.names = FALSE)

# Visualising the GO structure
## plot the induced subgraph to the current graphic device
# showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')

## save the resulting graph into a PDF or PS file.
printGraph(GOdata, resultFis, firstSigNodes = opt$firstSigNodes, fn.prefix = opt$name, useInfo = "all", pdfSW = TRUE)

