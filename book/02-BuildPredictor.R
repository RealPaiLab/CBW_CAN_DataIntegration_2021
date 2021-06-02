## ----echo=TRUE, fig.cap="Lab 2: Predictor design", out.width="50%"----
knitr::include_graphics("images/Lab1_design.jpg")


## ----eval=TRUE---------------------------------------------------
suppressMessages(library(curatedTCGAData))


## ----eval=TRUE---------------------------------------------------
curatedTCGAData(diseaseCode="BRCA", assays="*",dry.run=TRUE)


## ----eval=TRUE---------------------------------------------------
brca <- suppressMessages(curatedTCGAData("BRCA",
                                         c("mRNAArray","Methylation_methyl27", 
										 "RPPAArray","miRNASeqGene"),
                                         dry.run=FALSE))


## ----eval=TRUE---------------------------------------------------
brca


## ----eval=TRUE---------------------------------------------------
summary(assays(brca))


## ----eval=TRUE---------------------------------------------------
names(assays(brca))


## ----eval=TRUE---------------------------------------------------
mir <- assays(brca)[["BRCA_miRNASeqGene-20160128"]]
head(mir[,1:5])


## ----eval=TRUE---------------------------------------------------
pheno <- colData(brca)
colnames(pheno)[1:20]
head(pheno[,1:5])


## ----eval=TRUE---------------------------------------------------
source("prepare_data.R")
brca <- prepareDataForCBW(brca)


## ----eval=TRUE---------------------------------------------------
pheno <- colData(brca)
head(pheno[,c("ID","STATUS")])
table(pheno$STATUS,useNA="always")  # good practice: useNA="always" shows missing values


## ----eval=TRUE---------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))


## ----eval=TRUE---------------------------------------------------
expr <- assays(brca)
groupList <- list()
for (k in 1:length(expr)) {	# loop over all layers
	cur <- expr[[k]]; nm <- names(expr)[k]

	# all measures from this layer go into our single PSN
	groupList[[nm]] <- list(nm=rownames(cur)) 

	# assign same layer name as in input data
	names(groupList[[nm]])[1] <- nm;
}


## ----eval=TRUE---------------------------------------------------
summary(groupList)


## ----eval=TRUE---------------------------------------------------
names(groupList[["BRCA_mRNAArray-20160128"]])
length(groupList[["BRCA_mRNAArray-20160128"]][[1]])
head(groupList[["BRCA_mRNAArray-20160128"]][[1]])


## ---- eval=TRUE--------------------------------------------------
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() # initialize before is.null() check
	
	layerNames <- c("BRCA_miRNASeqGene-20160128",
		"BRCA_mRNAArray-20160128",
		"BRCA_RPPAArray-20160128",
		"BRCA_Methylation_methyl27-20160128")
	
	for (nm in layerNames){  			## for each layer
		if (!is.null(groupList[[nm]])){ ## must check for null for each layer
			netList_cur <- makePSN_NamedMatrix(
				dataList[[nm]],
				rownames(dataList[[nm]]),	## names of measures (e.g. genes, CpGs)
				groupList[[nm]],			## how to group measures in that layer
				netDir,						## leave this as-is, netDx will figure out where this is.
				verbose=FALSE, 			
				writeProfiles=TRUE,   		## use Pearson correlation-based similarity
				...
				)

			netList <- c(netList,netList_cur)	## just leave this in
		}
	}
	return(unlist(netList))	## just leave this in 
}



## ----eval=TRUE---------------------------------------------------
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
	"pred_output",sep=getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
t0 <- Sys.time()
model <- suppressMessages(buildPredictor(
	dataList=brca,			## your data
	groupList=groupList,	## grouping strategy
	makeNetFunc=makeNets,	## function to build PSNs
	outDir=outDir, 			## output directory
	trainProp=0.8,			## pct of samples to use to train model in
							## each split
	numSplits=2L,			## number of train/test splits
  	featSelCutoff=1L,		## threshold for calling something
							## feature-selected
  	featScoreMax=2L,		## max score for feature selection
    numCores=8L,			## set higher for parallelizing
  	debugMode=FALSE,
  	keepAllData=FALSE,	## set to TRUE for debugging or low-level files used by the predictor
    logging="none"
  ))
t1 <- Sys.time()
print(t1-t0)


## ----eval=TRUE---------------------------------------------------
summary(model)
summary(model$Split1)


## ----eval=TRUE---------------------------------------------------
outFile <- sprintf("%s/CBW_Lab1_full.rda",tempdir())
download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/Lab1_files/Lab1_10splits.rda",
	destfile=outFile)
lnames <- load(outFile)


## ----eval=TRUE---------------------------------------------------
source("helper.R")
results <- getResults(brca,model=model_full,featureSelCutoff=9L,
	featureSelPct=0.9)


## ----eval=TRUE---------------------------------------------------
summary(results)


## ----eval=TRUE---------------------------------------------------
results$performance


## ---- eval=TRUE--------------------------------------------------
results$featureScores


## ----eval=TRUE---------------------------------------------------
results$selectedFeatures


## ----fig.width=8,fig.height=8, eval=TRUE-------------------------
psnFile <- sprintf("%s/psn.rda",tempdir())
download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/Lab1_files/Lab1_PSN.rda",
	destfile=outFile)
load(outFile)
#psn <- getPSN(brca,groupList_full,makeNets_full,results$selectedFeatures)

tsne <- tSNEPlotter(
	psn$patientSimNetwork_unpruned, 
	colData(brca)
	)


## ----------------------------------------------------------------
sessionInfo()

