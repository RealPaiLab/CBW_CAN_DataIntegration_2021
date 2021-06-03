##  1. In practice I recommend running a predictor design with 2-3 different sets of pathway definitions, and comparing the predictive pathway themes. For instance, it could be useful to compare results from using all curated pathways, to just domain-specific ones.

##  Feature design is something of an art, and the choice of pathways depends on what your goals with building the predictor are. Are you looking to prioritize a known set of biological processes or interested in general discovery? These are tradeoffs.

##  At the very least I would recommend running with all curated pathways as a baseline, because you may generate novel hypotheses.

## 
##  2. Note that while in this design we group gene expression measures into pathways, the same design can be used to group other types of data based on prior knowledge. For instance, measures from imaging data could be grouped by prior knowledge of correlated networks of regions of interest subserving specific functions.

## 


## ---- class.source="codeblock",echo=FALSE, fig.cap="Lab 2 design: We will integrate clinical and gene expression data. Each layer will be converted into a single patient similarity network using Pearson correlation for pairwise similarity.", echo=FALSE----
knitr::include_graphics("images/Lab2_design.jpg")


## ---- class.source="codeblock", class.source="codeblock",eval=TRUE---------------------------------------
suppressMessages(library(curatedTCGAData))
brca <- suppressMessages(
   curatedTCGAData(
	   "BRCA",c("mRNAArray"),
	   dry.run=FALSE)
	)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
brca


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
source("prepare_data.R")
brca <- prepareDataForCBW(brca, setBinary=TRUE)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))
groupList <- list()


## This design changes model-building time to several hours, so avoid large gene sets (e.g. the full set of ~44,000 Gene Ontology terms, or even ~29,000 GO Biological Process terms). A reasonable start is a compilation of pathways from all curated pathway databases, as in below.

## Whichever list you use can be pruned by constraining the min/max number of genes in a set, but the size is something to keep in mind.


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
x <- fetchPathwayDefinitions("March",2021)
x


## ---- class.source="codeblock",echo=FALSE, fig.cap="Lab 2: Example of GMT file format."------------------
knitr::include_graphics("images/GMT_screenshot.png")


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
pathList <- readPathways(fetchPathwayDefinitions("March",2021))
head(pathList)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
gmtFile <- sprintf("%s2/cancer_pathways.gmt",tempdir())
if (!file.exists(sprintf("%s2",tempdir()))) {
	dir.create(sprintf("%s2",tempdir()))
}
download.file("https://raw.githubusercontent.com/RealPaiLab/CBW_CAN_DataIntegration_2021/master/supporting_files/c6.all.v7.4.symbols.gmt",gmtFile)
x <- readPathways(gmtFile)
x[1:3]


## The pathway definition file should use the same identifier type as your patient data. For instance, if the genes in your transcriptomic data are represented using [HGNC symbols](https://www.genenames.org/tools/search/), then your pathway definition file must also use HGNC symbols (e.g. *ID2S*), and not a different type of identifier, such as Ensembl IDs (which look like this: *ENSG00000010404*).


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
groupList[["BRCA_mRNAArray-20160128"]] <- pathList 


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
pheno <- colData(brca)
head(pheno[,c("patient.age_at_initial_pathologic_diagnosis","STAGE")])


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
groupList[["clinical"]] <- list(
      age="patient.age_at_initial_pathologic_diagnosis",
	   stage="STAGE"
)


## ---- class.source="codeblock",eval=FALSE----------------------------------------------------------------
## makePSN_NamedMatrix(..., writeProfiles=TRUE,...)`


## ---- class.source="codeblock",eval=FALSE----------------------------------------------------------------
## makePSN_NamedMatrix(,...,
## 	simMetric="custom", customFunc=normDiff,
## 	writeProfiles=FALSE)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() 

	# make RNA nets (Pearson correlation)
	rna <- "BRCA_mRNAArray-20160128"
	if (!is.null(groupList[[rna]])) { ## REMEMBER TO CHECK FOR NULL
		netList <- makePSN_NamedMatrix(
			dataList[[rna]],
			rownames(dataList[[rna]]),
			groupList[[rna]],
			netDir,
			verbose=FALSE,
			writeProfiles=TRUE,			## define Pearson similarity as before
			...) 
	}
	
	# make clinical nets (normalized difference)
	netList2 <- c()
	if (!is.null(groupList[["clinical"]])) {
	netList2 <- makePSN_NamedMatrix(
		dataList$clinical, 
		rownames(dataList$clinical),
		groupList[["clinical"]],netDir,
		simMetric="custom",customFunc=normDiff, ### Notice simMetric & customFunc
		writeProfiles=FALSE,
		sparsify=TRUE,
		verbose=FALSE,
		...)
	}
	netList <- c(unlist(netList),unlist(netList2))
	return(netList)
}


## ----lab2-buildpredictor,class.source="codeblock",eval=TRUE----------------------------------------------
t0 <- Sys.time()
set.seed(42) # make results reproducible
outDir <- paste(sprintf("%s2",tempdir()),"pred_output",sep=getFileSep()) # use absolute path
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
numSplits <- 2L
model <- suppressMessages(
   buildPredictor(
	   dataList=brca,
	   groupList=groupList,
	   makeNetFunc=makeNets,
	   outDir=outDir, 
	   numSplits=numSplits, 
	   featScoreMax=2L, 
	   featSelCutoff=1L,
	   numCores=4L
	   )
)
t1 <- Sys.time()
print(t1-t0) # time taken


## ----lab2-getresults, class.source="codeblock", eval=TRUE------------------------------------------------
outFile <- sprintf("%s2/CBW_Lab2_full.rda",tempdir())
download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/Lab2_files/brca_binary_pathways.rda",
	destfile=outFile)
lnames <- load(outFile)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
source("helper.R")
results <- getResults(brca,model_full,
	featureSelCutoff=9L,
	featureSelPct=0.9)


## ----class.source="codeblock",eval=TRUE------------------------------------------------------------------
perf <- results$performance
round(mean(perf$splitAUROC),2)*100
round(mean(perf$splitAUPR),2)*100
round(mean(perf$splitAccuracy),2)*100


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
confMat <- confusionMatrix(model_full)


## Particularly when there is class imbalance, i.e. one class has several-fold the number of samples than the other, predictors can achieve a baseline high accuracy by "cheating" and simply predicting all samples as being of the dominating label. e.g. if you have 99:1 class imbalance of two classes A and B, the predictor can achieve 99% accuracy simply by calling all samples of type "A"!

## 
## In practice, class imbalance should be handled using suitable performance evaluation metrics and sampling proportionally for training/test sets.


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
table(colData(brca)$STATUS,useNA="always")


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
summary(confMat)


## ----lab2-makeemapinput ,eval=TRUE-----------------------------------------------------------------------
emap <- makeInputForEnrichmentMap (
	model=model_full,
	results=results,
	pathwayList=pathList,
	EMapMinScore=7L, 
	EMapMaxSore=10L,
	EMapPctPass=0.7,
	outDir="/home/spai/data" ### SP: CHANGE TO WORKSPACE ON AWS
)


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
emap


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
###plotEmap(gmtFiles[[1]],nodeAttrFiles[[1]],
###         groupClusters=TRUE, hideNodeLabels=TRUE)


## ----echo=FALSE, out.width="140%",  fig.cap="Lab 2: EnrichmentMap for top-scoring pathways predictive of Luminal A status. Nodes show pathways scoring 7+ out of 10 in over 70% of train/test splits. Edges connect pathways with shared genes. Yellow bubbles show pathway clusters, labelled by a AutoAnnotate. Node fills indicate top score in 70%+ of the splits. Speech bubbles show pathways in example clusters."----
knitr::include_graphics("images/Lab2_EMap.jpg")


##  (1) In practice, it can take a good portion of an hour to adjust the layout of the EnrichmentMap, and often longer to explore the contents. The automatically-generated pathway theme labels are often revised upon looking closer inspection of the pathways within.


## ---- class.source="codeblock",eval=TRUE-----------------------------------------------------------------
sessionInfo()

