# setup brca data

prepareDataForCBW <- function(dat) {
### clean up stage variable
staget <- sub("[abcd]","",sub("t","",colData(dat)$pathology_T_stage))
staget <- suppressWarnings(as.integer(staget))
colData(dat)$STAGE <- staget

### remove NA PAM50 calls, remove normal samples
tmp <- colData(dat)$PAM50.mRNA
idx <- union(which(tmp %in% c("Normal-like","HER2-enriched")),
             		which(is.na(staget)))
idx <- union(idx, which(is.na(tmp)))				
pID <- colData(dat)$patientID
tokeep <- setdiff(pID, pID[idx])
dat <- dat[,tokeep,]
pam50 <- colData(dat)$PAM50.mRNA

### where a patient has multiple instances of the same assay
### just keep the first instance encountered
smp <- sampleMap(dat)
expr <- assays(dat)
for (k in 1:length(expr)) {
	samps <- smp[which(smp$assay==names(expr)[k]),]
	notdup <- samps[which(!duplicated(samps$primary)),"colname"]
	#message(sprintf("%s: %i notdup", names(expr)[k], length(notdup)))
	dat[[k]] <- suppressMessages(dat[[k]][,notdup])
}

### create ID, STATUS columns, remove spaces/hyphens from patient labels
pID <- colData(dat)$patientID
colData(dat)$ID <- pID
colData(dat)$STATUS <- pam50
colData(dat)$STATUS <- gsub(" ",".",colData(dat)$STATUS)
colData(dat)$STATUS <- gsub("-",".",colData(dat)$STATUS)

return(dat)
}

