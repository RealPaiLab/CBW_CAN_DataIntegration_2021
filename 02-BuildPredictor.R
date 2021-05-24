## ----eval=FALSE---------------------------------------------------------------
## numSplits <- 100     # num times to split data into train/blind test samples
## featScoreMax <- 10   # num folds for cross-validation, also max score for a network
## featSelCutoff <- 9
## netScores <- list()  # collect <numSplits> set of netScores
## perf <- list()       # collect <numSplits> set of test evaluations
## 
## for k in 1:numSplits
##  [train, test] <- splitData(80:20) # split data using RNG seed
##   featScores[[k]] <- scoreFeatures(train, featScoreMax)
##  topFeat[[k]] <- applyFeatCutoff(featScores[[k]])
##  perf[[k]] <- collectPerformance(topFeat[[k]], test)
## end


## ----eval=TRUE----------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))


## ----eval=TRUE----------------------------------------------------------------
suppressMessages(library(curatedTCGAData))


## ----eval=TRUE----------------------------------------------------------------
curatedTCGAData(diseaseCode="BRCA", assays="*",dry.run=TRUE)


## ----eval=TRUE----------------------------------------------------------------
brca <- suppressMessages(curatedTCGAData("BRCA",
                                         c("mRNAArray","Methylation_methyl27", "RPPAArray","miRNASeqGene"),
                                         dry.run=FALSE))


## ----eval=TRUE----------------------------------------------------------------
brca


## ----eval=TRUE----------------------------------------------------------------
pheno <- colData(brca)
colnames(pheno)[1:20]
head(pheno[,1:5])


## ----eval=TRUE----------------------------------------------------------------
staget <- sub("[abcd]","",sub("t","",colData(brca)$pathology_T_stage))
staget <- suppressWarnings(as.integer(staget))
colData(brca)$STAGE <- staget

pam50 <- colData(brca)$PAM50.mRNA
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
colData(brca)$pam_mod <- pam50

tmp <- colData(brca)$PAM50.mRNA
idx <- union(which(tmp %in% c("Normal-like","Luminal B","HER2-enriched")),
             		which(is.na(staget)))
pID <- colData(brca)$patientID
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

## remove duplicate assays mapped to the same sample
smp <- sampleMap(brca)
expr <- assays(brca)
for (k in 1:length(expr)) {
	samps <- smp[which(smp$assay==names(expr)[k]),]
	notdup <- samps[which(!duplicated(samps$primary)),"colname"]
	brca[[k]] <- suppressMessages(brca[[k]][,notdup])
}


## ----eval=TRUE----------------------------------------------------------------
pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- colData(brca)$pam_mod


## ----eval=FALSE---------------------------------------------------------------
## summary(brca)


## ----eval=TRUE----------------------------------------------------------------
groupList <- list()

for (k in 1:length(expr)) {
	cur <- expr[[k]]; nm <- names(expr)[k]
	# all measure names should be in rownames column
	groupList[[nm]] <- list(nm=rownames(cur)) 
	names(groupList[[nm]])[1] <- nm;
}


## ----eval=TRUE----------------------------------------------------------------
summary(groupList)


## ----eval=FALSE---------------------------------------------------------------
## names(groupList[["BRCA_mRNAArray-20160128"]])
## head(groupList[["BRCA_mRNAArray-20160128"]][[1]])


## ---- eval=TRUE---------------------------------------------------------------
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



## ----eval=TRUE----------------------------------------------------------------
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
	"pred_output",sep=getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
## set keepAllData=TRUE to not delete at the end of the predictor run.
## This can be useful for debugging.
out <- buildPredictor(
      dataList=brca,		## your data
	  groupList=groupList,	## grouping strategy
      makeNetFunc=makeNets,	## function to build PSNs
      outDir=outDir, 		## output directory
      numSplits=2L,			## number of train/test splits
	  featScoreMax=2L,		## max score for feature selection
      featSelCutoff=1L,		## threshold for calling something feature-selected
      numCores=1L,			## set higher for parallelizing
	  debugMode=FALSE,
      logging="default"
	  )


## ----eval=FALSE---------------------------------------------------------------
## summary(out)
## summary(out$Split1)


## ----eval=FALSE---------------------------------------------------------------
## numSplits <- 2
## st <- unique(colData(brca)$STATUS)
## acc <- c()         # accuracy
## predList <- list() # prediction tables
## 
## featScores <- list() # feature scores per class
## for (cur in unique(st)) featScores[[cur]] <- list()
## 
## for (k in 1:numSplits) {
## 	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
## 	# predictions table
## 	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
## 	                 sprintf("%s_SCORE",st))]
## 	predList[[k]] <- tmp
## 	# accuracy
## 	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
## 	# feature scores
## 	for (cur in unique(st)) {
## 	   tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
## 	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
## 	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
## 	}
## }


## ----eval=FALSE---------------------------------------------------------------
## print(acc)


## ----eval=FALSE,fig.width=8,fig.height=10-------------------------------------
## predPerf <- plotPerf(predList, predClasses=st)


## ----eval=FALSE---------------------------------------------------------------
## featScores2 <- lapply(featScores, getNetConsensus)
## summary(featScores2)
## head(featScores2[["LumA"]])


## ----eval=FALSE---------------------------------------------------------------
## featSelNet <- lapply(featScores2, function(x) {
##     callFeatSel(x, fsCutoff=1, fsPctPass=0)
## })
## print(head(featScores2[["LumA"]]))


## ----eval=FALSE---------------------------------------------------------------
## Emap_res <- getEMapInput_many(featScores2,pathList,
##     minScore=1,maxScore=2,pctPass=0,out$inputNets,verbose=FALSE)


## ----eval=FALSE---------------------------------------------------------------
## gmtFiles <- list()
## nodeAttrFiles <- list()
## 
## for (g in names(Emap_res)) {
##     outFile <- paste(outDir,sprintf("%s_nodeAttrs.txt",g),sep=getFileSep())
##     write.table(Emap_res[[g]][["nodeAttrs"]],file=outFile,
##         sep="\t",col=TRUE,row=FALSE,quote=FALSE)
##     nodeAttrFiles[[g]] <- outFile
## 
##     outFile <- paste(outDir,sprintf("%s.gmt",g),sep=getFileSep())
##     conn <- suppressWarnings(
##          suppressMessages(base::file(outFile,"w")))
##     tmp <- Emap_res[[g]][["featureSets"]]
##     gmtFiles[[g]] <- outFile
## 
##     for (cur in names(tmp)) {
##         curr <- sprintf("%s\t%s\t%s", cur,cur,
##             paste(tmp[[cur]],collapse="\t"))
##         writeLines(curr,con=conn)
##     }
## close(conn)
## }


## ----eval=FALSE---------------------------------------------------------------
## ###plotEmap(gmtFiles[[1]],nodeAttrFiles[[1]],
## ###         groupClusters=TRUE, hideNodeLabels=TRUE)
## 


## ----eval=FALSE---------------------------------------------------------------
## featScores2 <- lapply(featScores, getNetConsensus)
## featSelNet <- lapply(featScores2, function(x) {
##     callFeatSel(x, fsCutoff=2, fsPctPass=1)
## })


## ----eval=FALSE---------------------------------------------------------------
## print(featSelNet)


## ----eval=FALSE---------------------------------------------------------------
## topPath <- gsub(".profile","",
## 		unique(unlist(featSelNet)))
## topPath <- gsub("_cont.txt","",topPath)
## ## create groupList limited to top features
## g2 <- list();
## for (nm in names(groupList)) {
## 	cur <- groupList[[nm]]
## 	idx <- which(names(cur) %in% topPath)
## 	message(sprintf("%s: %i pathways", nm, length(idx)))
## 	if (length(idx)>0) g2[[nm]] <- cur[idx]
## }


## ----eval=FALSE---------------------------------------------------------------
## psn <- suppressMessages(
##    plotIntegratedPatientNetwork(brca,
##   groupList=g2, makeNetFunc=makeNets,
##   aggFun="MEAN",prune_pctX=0.30,prune_useTop=TRUE,
##   numCores=1L,calcShortestPath=TRUE,
##   showStats=FALSE,
##   verbose=FALSE, plotCytoscape=FALSE)
## )


## ----fig.width=8,fig.height=8, eval=FALSE-------------------------------------
## tsne <- plot_tSNE(psn$patientSimNetwork_unpruned,colData(brca))
## summary(tsne)
## class(tsne)


## -----------------------------------------------------------------------------
sessionInfo()

