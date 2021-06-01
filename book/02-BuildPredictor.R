## ----eval=TRUE--------------------------------------------------------------
suppressMessages(library(curatedTCGAData))
## ----eval=TRUE--------------------------------------------------------------
brca <- suppressMessages(curatedTCGAData("BRCA",
                                         c("mRNAArray","Methylation_methyl27", 
										 "RPPAArray","miRNASeqGene"),
                                         dry.run=FALSE))
## ----eval=TRUE--------------------------------------------------------------

## ----eval=TRUE--------------------------------------------------------------
source("prepare_data.R")
brca <- prepareDataForCBW(brca)
pheno <- colData(brca)

## ----eval=TRUE--------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))


## ----eval=TRUE--------------------------------------------------------------
expr <- assays(brca)
groupList <- list()
for (k in 1:length(expr)) {	# loop over all layers
	cur <- expr[[k]]; nm <- names(expr)[k]

	# all measures from this layer go into our single PSN
	groupList[[nm]] <- list(nm=rownames(cur)) 

	# assign same layer name as in input data
	names(groupList[[nm]])[1] <- nm;
}

## ---- eval=TRUE-------------------------------------------------------------
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
				verbose=TRUE, 			
				writeProfiles=TRUE,   		## use Pearson correlation-based similarity
				...
				)

			netList <- c(netList,netList_cur)	## just leave this in
		}
	}
	return(unlist(netList))	## just leave this in 
}



## ----eval=TRUE--------------------------------------------------------------
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
	"pred_output",sep=getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
t0 <- Sys.time()
model <- buildPredictor(
	dataList=brca,			## your data
	groupList=groupList,	## grouping strategy
	makeNetFunc=makeNets,	## function to build PSNs
	outDir=outDir, 			## output directory
	trainProp=0.8,			## pct of samples to use to train model in
							## each split
	numSplits=10L,			## number of train/test splits
  	featSelCutoff=9L,		## threshold for calling something
							## feature-selected
  	featScoreMax=10L,		## max score for feature selection
    numCores=10L,			## set higher for parallelizing
  	debugMode=FALSE,
  	keepAllData=FALSE,	## set to TRUE for debugging or low-level files used by the predictor
    logging="none"
  )
t1 <- Sys.time()
print(t1-t0)


## ----eval=TRUE--------------------------------------------------------------
source("helper.R")
results <- getResults(brca,model,featureSelCutoff=2L,
	featureSelPct=0.5)

## ----fig.width=8,fig.height=8, eval=TRUE------------------------------------
psn <- getPSN(brca,groupList,makeNets,results$selectedFeatures)

tsne <- tSNEPlotter(
	psn$patientSimNetwork_unpruned, 
	colData(brca)
	)


## ---------------------------------------------------------------------------
sessionInfo()

