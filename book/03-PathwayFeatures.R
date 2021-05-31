## ----eval=TRUE----------------------------------------------------------------
### suppressMessages(library(curatedTCGAData))
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### brca <- suppressMessages(
###    curatedTCGAData(
### 	   "BRCA",c("mRNAArray"),
### 	   dry.run=FALSE)
### 	)
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### source("prepare_data.R")
### brca <- prepareDataForCBW(brca)
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### suppressWarnings(suppressMessages(require(netDx)))
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### groupList <- list()
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### x <- fetchPathwayDefinitions("March",2021)
### x
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### ###system2("head",args=c("-n","3", x))
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### message("here")
### pathList <- readPathways(fetchPathwayDefinitions("March",2021))
### head(pathList)
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### #####TBA
### ###myGMTfile <- "supporting_files/c6.all.c7.4.symbols.gmt"
### ###file.exists(myGMTfile)
### ####x <- readPathways(myGMTfile)
### ###head(x)
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### message("grouplist")
### groupList[["BRCA_mRNAArray-20160128"]] <- pathList
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### message("pheno")
### pheno <- colData(brca)
### head(pheno[,c("patient.age_at_initial_pathologic_diagnosis","STAGE")])
### 
### 
### ## ----eval=TRUE----------------------------------------------------------------
### message("clinical")
### groupList[["clinical"]] <- list(
###       age="patient.age_at_initial_pathologic_diagnosis",
### 	   stage="STAGE"
### )
### 
### 
### ## ----eval=FALSE---------------------------------------------------------------
### ## makePSN_NamedMatrix(..., writeProfiles=TRUE,...)`
### 
### 
### ## ----eval=FALSE---------------------------------------------------------------
### ##    makePSN_NamedMatrix(,...,
### ##                        simMetric="custom", customFunc=normDiff,
### ##                        writeProfiles=FALSE)
### 
### 
### ## ----eval=FALSE---------------------------------------------------------------
### makeNets <- function(dataList, groupList, netDir,...) {
### 	netList <- c()
### 	for (nm in setdiff(names(groupList),"clinical")) {
### 		if (!is.null(groupList[[nm]])) { ## REMEMBER TO CHECK FOR NULL
### 		netList <- makePSN_NamedMatrix(
### 			dataList[[nm]],
### 			rownames(dataList[[nm]]),
### 			groupList[[nm]],
### 			netDir,
### 			verbose=FALSE,
### 			writeProfiles=TRUE,
### 			...)
### 		}
### 	}
### 	
### 	# make clinical nets (normalized difference)
### 	netList2 <- c()
### 	if (!is.null(groupList[["clinical"]])) {
### 	netList2 <- makePSN_NamedMatrix(
### 		dataList$clinical,
### 		rownames(dataList$clinical),
### 		groupList[["clinical"]],netDir,
### 		simMetric="custom",customFunc=normDiff, ### Notice simMetric & customFunc
### 		writeProfiles=FALSE,
### 		sparsify=TRUE,verbose=TRUE,...)
### 	}
### 	netList <- c(unlist(netList),unlist(netList2))
### 	return(netList)
### }
### 
### 
### 
### 
### set.seed(42) # make results reproducible
### outDir <- paste(tempdir(),"pred_output",sep=getFileSep()) # use absolute path
### if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
### numSplits <- 10L
### t0 <- Sys.time()
### model <- suppressMessages(
###    buildPredictor(
### 	   dataList=brca,
### 	   groupList=groupList,
### 	   makeNetFunc=makeNets,
### 	   outDir=outDir,
### 	   numSplits=numSplits,
### 	   featScoreMax=10L,
### 	   featSelCutoff=9L,
### 	   numCores=12L
### 	   )
### )
### print(Sys.time()-t0)
### save(brca,model,file="/home/spai/data/brca_pathways.rda")

load("/home/spai/data/brca_pathways.rda")

## ----eval=FALSE---------------------------------------------------------------
 source("helper.R")
 results <- getResults(brca,model,featureSelCutoff=8L,
 	featureSelPct=0.7)


## ----eval=FALSE---------------------------------------------------------------
## print("get confusion matrix here")


## ----eval=TRUE----------------------------------------------------------------
###plotEmap(gmtFiles[[1]],nodeAttrFiles[[1]],
###         groupClusters=TRUE, hideNodeLabels=TRUE)


## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

