#===============================================================================
# author: Benedicte Rhone
# date: 03/2016 
#-------------------------------------------------------------------------------
# Description: GOterms enrichment analysis on a set of genomic region under 
# selection using the topGO R package
#
# Citation: Alexa A and Rahnenfuhrer J (2016). topGO: Enrichment Analysis for 
# Gene Ontology. R package version 2.22.0

#-------------------------------------------------------------------------------
# Inmput files description: 2 files
# GOterms-anotation-file.txt: A table containing at least the gene ID and the 
# corresponding GO terms annotation
# datafile.txt: A file containing at least a list of the name of the genes under 
# selection
#
# Steps of the script:
# 1) Read Input data
# 2) Data preparation for topGO analysis
# 3) topGO Analysis
#===============================================================================

### Install packages
source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
biocLite("Rgraphviz")

### Set the working directory
setwd("..........")

### Load packages
library(topGO)
library(Rgraphviz)

#######################
# 1) Read Input data  #
#######################

#read annotation file
Annot<-read.table("GOterms-anotation-file.txt", header=FALSE, sep="\t", fill=TRUE)
Annot<-Annot[,c(1,2)]
names(Annot)<-c("geneid","GOterm")

#return number of annotated genes
length(unique(Annot$geneid)) 	

###datafile : list of genes identified as under selection
data<-read.table("datafile.txt", header=FALSE, sep="\t", fill=TRUE)
names(data)<-c("geneid")
dataDF<-as.data.frame(data)

#return no. genes found as selected
length(unique(dataDF$geneid))

###########################################
# 2) Data preparation for topGO analysis  #
###########################################
# list of the annotated genes  - eliminate redundancy
ListAnnotUnique<-as.vector(unique(Annot$geneid))
geneID2GO <- list(NULL)

for (i in 1 : length(ListAnnotUnique)) {
  temp<-Annot[Annot$geneid==ListAnnotUnique[i],]
  geneID2GO[[i]]<-as.character(temp$GOterm)
}
names(geneID2GO)<-as.character(ListAnnotUnique)

#Building the list of genes of interest
geneNames2GO <- names(geneID2GO)
geneListdataGO <- factor(as.integer(geneNames2GO %in% dataDF$geneid))
names(geneListdataGO) <- geneNames2GO
length(which(geneListdataGO==1))  #returns the number of genes under selection and annotated thus usable for the analysis
str(geneListdataGO)


######################
# 3) topGO Analysis  #
######################
# creation of a topGO data object #### 
GOdata <- new("topGOdata",
			  description = "GO analysis on genomic regions under selection, BP",
              ontology = "BP", # or CC or MF
              allGenes = geneListdataGO,
              nodeSize =5,   # delete categories that have too few genes : in our case, a high number of genes gives similar results with a nodeSize of 5 or of 10 (in the tutorial: values between 5 and 10 give more stable results) 
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO
)

#returns the description of the topGO object
GOdata

#### 5 types of statistical tests and 6 algorithms dealing with the GO graph structure are available in the topGO R package (See description in the tutorial)
#### Here we used a Fisher exact test based on genes counting combine with 2 algorithms 

#### Fisher test with Classic algorithm not taking into account the hierarchical link between GOterms
Fisherclassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

#returns a table listing the top 20 GO terms found as significantly enriched
Fisherclassic.table<-GenTable(GOdata, classicFisher = Fisherclassic, topNodes=20)
Fisherclassic.table
write.table(Fisherclassic.table, "resultFisherclassic.table.txt" , sep=";", quote=FALSE)

#returns a subgraph induced by the top 4 GO terms found as significantly enriched 
showSigOfNodes(GOdata, score(Fisherclassic), firstSigNodes = 4, useInfo ='all')
#generates a pdf file of the subgraph induced by the top 5 GO terms found as significantly enriched 
printGraph(GOdata, Fisherclassic, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE)


#### Fisher test with weight01 algorithm taking into account the hierarchical link between GOterms
FisherWeight01<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#returns a table listing the top 50 GO terms found as significantly enriched
FisherWeight01.table<-GenTable(GOdata, Weight01 = FisherWeight01, topNodes=50)
FisherWeight01.table
write.table(FisherWeight01.table, "FisherWeight01.table.txt" , sep=";", quote=FALSE)

#returns a subgraph induced by the top 4 GO terms found as significantly enriched 
showSigOfNodes(GOdata, score(FisherWeight01), firstSigNodes = 4, useInfo ='all')
#generates a pdf file of the subgraph induced by the top 5 GO terms found as significantly enriched 
printGraph(GOdata, FisherWeight01, firstSigNodes = 5, useInfo = "all", pdfSW = TRUE)

#returns a summary table of the top 10 GO terms found as significantly enriched with the FisherWeight01 compared to the Fisherclassic tests of enrichment
allRes<-GenTable(GOdata, classicFisher = Fisherclassic, weight01=FisherWeight01, classic=Fisherclassic, orderBy="weight01", ranksOf="classic", topNodes=10)
allRes   
write.table(allRes  , "allRes.pcadapt.txt" , sep=";", quote=FALSE)

