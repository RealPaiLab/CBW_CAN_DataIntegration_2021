# CBW helper functions


#' prints performance, returns list of output.
#'
#' @param dat (MultiAssayExperiment) input data
#' @param res (list) output of buildPredictor() function
#' @param featureSelCutoff (integer) cutoff score for feature selection.
#' A feature must have minimum of this score for specified fraction of splits 
#' (see featureSelPct) to pass.
#' @param featureSelPct (numeric between 0 and 1) cutoff percent for feature selection.
#' A feature must have minimum score of featureSelCutoff for featureSelPct of 
#' train/test splits to pass.
#' @returns list of results.
#' - selectedFeatures (list of character vectors): list, one per class
#' - performance (list of mixed datatypes) accuracy, confusion matrix, AUROC
#' Side effect of plotting ROC curve if binary classifier
#' @export
getResults <- function(dat, res, featureSelCutoff=1L, 
    featureSelPct=0, pathwayList=NULL){

numSplits <- length(grep("^Split",names(res)))
st <- unique(colData(dat)$STATUS)
message(sprintf("Detected %i splits and %i classes", numSplits, length(st)))

acc <- c()         # accuracy
predList <- list() # prediction tables
featScores <- list() # feature scores per class
for (cur in unique(st)) featScores[[cur]] <- list()

# collect accuracy and feature scores
for (k in 1:numSplits) {
	pred <- res[[sprintf("Split%i",k)]][["predictions"]];
	# predictions table
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
	                 sprintf("%s_SCORE",st))]
	predList[[k]] <- tmp
	# accuracy
	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
	# feature scores
	for (cur in unique(st)) {
	   tmp <- res[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
	}
}
message("* Plotting performance")
predPerf <- plotPerf(predList, predClasses=st)

message("* Compiling feature scores and calling selected features")
feats <- runOverallFeatureSelection(featScores, 
    featureSelCutoff = featureSelCutoff,
    featureSelPct = featureSelPct,
    cleanNames = TRUE
)

# Enrichment map
if (!is.null(pathwayList)){
    message("* Pathway List detected - creating input for EnrichmentMap")
    browser()
}

return(list(
    selectedFeatures=feats,
    featureScores=feats$featScores,
    performance=list(meanAccuracy=mean(acc),
                    splitAccuracy=acc)
))

}

#' wrapper to call selected features
#'
#' @param featureSelCutoff (integer) cutoff score for feature selection.
#' A feature must have minimum of this score for specified fraction of splits 
#' (see featureSelPct) to pass.
#' @param featureSelPct (numeric between 0 and 1) cutoff percent for feature selection.
#' A feature must have minimum score of featureSelCutoff for featureSelPct of 
#' train/test splits to pass.
#' @param cleanNames (logical) remove internal suffixes for human readability
#' @return (list) Feature scores for all splits, plus those passing selection for overall predictor
#' featScores: (matrix) feature scores for each split
#' selectedFeatures: (list) features passing selection for each class; one key per class
runOverallFeatureSelection <- function(featScores, featureSelCutoff, 
    featureSelPct, printRes=TRUE, cleanNames=TRUE){
featScores2 <- lapply(featScores, getNetConsensus)
if (cleanNames) {
    featScores2 <- lapply(featScores2,function(x){
        x$PATHWAY_NAME <- sub(".profile","",x$PATHWAY_NAME)
        x$PATHWAY_NAME <- sub("_cont.txt","",x$PATHWAY_NAME)
        colnames(x)[1] <- "Feature"
        x
    })
}
featSelNet <- lapply(featScores2, function(x) {
    x <- callFeatSel(x, fsCutoff=1, fsPctPass=0)
})


return(list(
    featScores=featScores2,
    selectedFeatures=featSelNet
))
}

#' wrapper to create input files for Enrichment Map
#'
#' @param res (list) output of buildPredictor() function
#' @param featScores (list) keys are classes, values are matrices with feature scores for all 
#' train/test splits. Output of runOverallFeatureSelection()
#' @param pathwayList (list) output of readPathwayFile() used to make pathway-level features for predictor
#' @param EMapMinScore (integer) minimum score for Enrichment Map
#' @param EMapMaxScore (integer) maximum score for Enrichment Map
#' @param featureSelPct (numeric between 0 and 1) percent of splits for which feature must have score in range
#'  [EMapMinScore,EMapMaxScore] to be included for EnrichmentMap visualization
#' @return 
makeInputForEnrichmentMap <- function(res,featScores,pathwayList,EMapMinScore=0L, EMapMaxSore=1L,
    featureSelPct)
{
message("* Getting EnrichmentMap input")
Emap_res <- getEMapInput_many(featScores,
    pathwayList,
    minScore=EMapMinScore,
    maxScore=EMapMaxSore,
    pctPass=featureSelPct,
    res$inputNets,
    verbose=FALSE
)
gmtFiles <- list()
nodeAttrFiles <- list()

message("* Writing files for network visualization")
for (g in names(Emap_res)) {
    outFile <- paste(outDir,sprintf("%s_nodeAttrs.txt",g),sep=getFileSep())
    write.table(Emap_res[[g]][["nodeAttrs"]],file=outFile,
        sep="\t",col=TRUE,row=FALSE,quote=FALSE)
    nodeAttrFiles[[g]] <- outFile

    outFile <- paste(outDir,sprintf("%s.gmt",g),sep=getFileSep())
    conn <- suppressWarnings(
         suppressMessages(base::file(outFile,"w")))
    tmp <- Emap_res[[g]][["featureSets"]]
    gmtFiles[[g]] <- outFile

    for (cur in names(tmp)) {
        curr <- sprintf("%s\t%s\t%s", cur,cur,
            paste(tmp[[cur]],collapse="\t"))
        writeLines(curr,con=conn)
    }
close(conn)
}
}

#' get the integrated patient similarity network made of selected features
#' @param dat (MultiAssayExperiment) input data
#' @param res (list) output of buildPredictor() function
#' @param aggFun - see plotIntegratedPatientNetwork()
#' @param prune_pctX - see plotIntegratedPatientNetwork()
#' @param prune_useTop - see plotIntegratedPatientNetwork()
#' @param numCores - see plotIntegratedPatientNetwork()
#' @param calcShortestPath - see plotIntegratedPatientNetwork()
#' @return (list) output of plotIntegratedPatientNetwork()
getPSN <- function(dat, groupList, makeNets, selectedFeatures, plotCytoscape=FALSE,
    aggFun="MEAN", prune_pctX=0.30, prune_useTop=TRUE,numCores=1L,calcShortestPath=FALSE
    ){
    
topPath <- gsub(".profile","", unique(unlist(selectedFeatures)))
topPath <- gsub("_cont.txt","",topPath)

## create groupList limited to top features
g2 <- list();
for (nm in names(groupList)) {
	cur <- groupList[[nm]]
	idx <- which(names(cur) %in% topPath)
	message(sprintf("%s: %i features", nm, length(idx)))
	if (length(idx)>0) g2[[nm]] <- cur[idx]
}

message("* Making integrated PSN")
psn <- suppressMessages(
   plotIntegratedPatientNetwork(brca,
  groupList=g2, makeNetFunc=makeNets,
  aggFun=aggFun,
  prune_pctX=prune_pctX,
  prune_useTop=prune_useTop,
  numCores=numCores,
  calcShortestPath=calcShortestPath,
  showStats=FALSE,
  verbose=TRUE, 
  plotCytoscape=plotCytoscape)
)

return(psn)
}

#' creates a MultiAssayExperiment object from user-provided files
#' 
#' @param dataTables (list of matrices/data.frames) one key per data type, 
#' with table in corresponding value. Keys are user-defined names for the data layers
#' @param pheno (data.frame) sample metadata. Rows are patients and columns are variables.
#' * Patient IDs must be in rownames
createMultiAssayExperiment <- function(dataTables,pheno){
    summList <- lapply(dataTables, function(x){
        SummarizedExperiment(x)
    })

    return(MultiAssayExperiment(experiments=summList,
                        colData=pheno)
    )
}
