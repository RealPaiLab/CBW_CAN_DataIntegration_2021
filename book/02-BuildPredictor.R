rm(list=ls())
suppressWarnings(suppressMessages(require(netDx)))
suppressMessages(library(curatedTCGAData))
brca <- suppressMessages(curatedTCGAData("BRCA",
                                         c("mRNAArray","Methylation_methyl27", 
										 "RPPAArray","miRNASeqGene"),
                                         dry.run=FALSE))
source("prepare_data.R")
brca <- prepareDataForCBW(brca)

expr <- assays(brca)
groupList <- list()
for (k in 1:length(expr)) {
	cur <- expr[[k]]; nm <- names(expr)[k]
	# all measure names should be in rownames column
	groupList[[nm]] <- list(nm=rownames(cur)) 
	names(groupList[[nm]])[1] <- nm;
}

makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() # initialize before is.null() check
	
	layerNames <- c("BRCA_miRNASeqGene-20160128",
		"BRCA_mRNAArray-20160128",
		"BRCA_RPPAArray-20160128",
		"BRCA_Methylation_methyl27-20160128")

	# Similarity defined as Pearson correlation
	for (nm in layerNames){
		if (!is.null(groupList[[nm]])){ # IMPORANT check
			netList_cur <- makePSN_NamedMatrix(dataList[[nm]],
				rownames(dataList[[nm]]),
				groupList[[nm]],
				netDir,verbose=TRUE, 
				writeProfiles=TRUE,...)
			netList <- c(netList,netList_cur)
		}

	}
	return(unlist(netList))
}


set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
	"pred_output",sep=getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
##set keepAllData=TRUE to not delete at the end of the predictor run.
##This can be useful for debugging.
t0 <- Sys.time()
model <- buildPredictor(
      dataList=brca,		## your data
	  groupList=groupList,	## grouping strategy
      makeNetFunc=makeNets,	## function to build PSNs
      outDir=outDir, 		## output directory
      numSplits=2L,			## number of train/test splits
	  featScoreMax=2L,		## max score for feature selection
      featSelCutoff=1L,		## threshold for calling something feature-selected
      numCores=4L,			## set higher for parallelizing
	  debugMode=FALSE,
	  keepAllData=FALSE,
      logging="none"
	  )
print(Sys.time()-t0)

save(brca, model,file=sprintf("%s/output.rda",outDir))

source("helper.R")
results <- getResults(brca,model,featureSelCutoff=7L,
	featureSelPct=0.7)

psn <- getPSN(brca,groupList,makeNets,results$selectedFeatures)
#tsne <- plot_tSNE(psn$patientSimNetwork_unpruned,
#	colData(brca))
tsne <- tSNEPlotter(
	psn$patientSimNetwork_unpruned, 
	colData(brca)
	)



sessionInfo()

