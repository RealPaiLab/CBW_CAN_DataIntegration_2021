## ----echo=FALSE, fig.cap="Lab 2: Predictor design"---------------
knitr::include_graphics("images/Lab2_design.jpg")


## ----load-data,eval=TRUE-----------------------------------------
suppressMessages(library(curatedTCGAData))
brca <- suppressMessages(
   curatedTCGAData(
	   "BRCA",c("mRNAArray"),
	   dry.run=FALSE)
	)


## ----eval=TRUE---------------------------------------------------
brca


## ----eval=TRUE---------------------------------------------------
source("prepare_data.R")
brca <- prepareDataForCBW(brca, setBinary=TRUE)


## ----eval=TRUE---------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))
groupList <- list()


## ----eval=TRUE---------------------------------------------------
x <- fetchPathwayDefinitions("March",2021)
x


## ----echo=FALSE, fig.cap="Lab 2: Predictor design"---------------
knitr::include_graphics("images/GMT_screenshot.png")


## ----eval=TRUE---------------------------------------------------
pathList <- readPathways(fetchPathwayDefinitions("March",2021))
head(pathList)


## ----eval=TRUE---------------------------------------------------
gmtFile <- sprintf("%s/cancer_pathways.gmt",tempdir())
download.file("https://raw.githubusercontent.com/RealPaiLab/CBW_CAN_DataIntegration_2021/master/supporting_files/c6.all.v7.4.symbols.gmt",gmtFile)
x <- readPathways(gmtFile)
x[1:3]


## ----eval=TRUE---------------------------------------------------
groupList[["BRCA_mRNAArray-20160128"]] <- pathList 


## ----eval=TRUE---------------------------------------------------
pheno <- colData(brca)
head(pheno[,c("patient.age_at_initial_pathologic_diagnosis","STAGE")])


## ----eval=TRUE---------------------------------------------------
groupList[["clinical"]] <- list(
      age="patient.age_at_initial_pathologic_diagnosis",
	   stage="STAGE"
)


## ----eval=FALSE--------------------------------------------------
## makePSN_NamedMatrix(..., writeProfiles=TRUE,...)`


## ----eval=FALSE--------------------------------------------------
## makePSN_NamedMatrix(,...,
## 	simMetric="custom", customFunc=normDiff,
## 	writeProfiles=FALSE)


## ----eval=TRUE---------------------------------------------------
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


## ----eval=TRUE---------------------------------------------------
t0 <- Sys.time()
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),"pred_output",sep=getFileSep()) # use absolute path
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


## ----eval=TRUE---------------------------------------------------
outFile <- sprintf("%s/CBW_Lab2_full.rda",tempdir())
download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/gmt_pathway_2021/brca_binary_pathways.rda",
	destfile=outFile)
lnames <- load(outFile)


## ----eval=TRUE---------------------------------------------------
source("helper.R")
results <- getResults(brca,model_full,
	featureSelCutoff=9L,
	featureSelPct=0.9)


## ----eval=TRUE---------------------------------------------------
confMat <- confusionMatrix(model_full)


## ----eval=TRUE---------------------------------------------------
summary(confMat)


## ----eval=TRUE---------------------------------------------------
emap <- makeInputForEnrichmentMap (
	model=model_full,
	results=results,
	pathwayList=pathList_full,
	EMapMinScore=7L, 
	EMapMaxSore=10L,
	EMapPctPass=0.7,
	outDir="/home/spai/data" ### SP: CHANGE TO WORKSPACE ON AWS
)


## ----eval=TRUE---------------------------------------------------
emap


## ----eval=TRUE---------------------------------------------------
###plotEmap(gmtFiles[[1]],nodeAttrFiles[[1]],
###         groupClusters=TRUE, hideNodeLabels=TRUE)


## ----eval=TRUE---------------------------------------------------
knitr::include_graphics("images/Lab2_EMap.jpg")


## ----eval=TRUE---------------------------------------------------
sessionInfo()

