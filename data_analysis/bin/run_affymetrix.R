##----------------------------
## Read and parse arguments
##----------------------------
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

args <- commandArgs(TRUE)

## test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least 2 arguments must be supplied: --targets_file=targets.txt --contrast=A-B", call.=FALSE)    
} else if (length(args)>1) {
    argsDF <- as.data.frame(do.call(rbind, parseArgs(args)))
    args <- as.list(as.character(argsDF$V2))
    names(args) <- argsDF$V1
}



##----------------------------
## Libraries
##----------------------------
##library(affycoretools) ##Needed only for QC plots
.libPaths("/nfs/users/nfs_r/rb11/local/lib/lib-r/3.2.2")
.libPaths()
##(.packages(all.available=TRUE))

##----------------------------
## Read data
##----------------------------
library(limma, verbose=FALSE, lib.loc="/nfs/users/nfs_r/rb11/local/lib/lib-r/3.2.2")

cat("\nReading data ...\n\n")

wd <- getwd()

if (file.exists(args$targets_file)) {
    targetsPath <- dirname(args$targets_file)
    targetsFile <- basename(args$targets_file)
    targets <- readTargets(file=targetsFile, path=targetsPath)
} else {
    cat(paste("Read targets file: file '", args$targets_file, "' could not be found, wd=", wd, "\n", sep=""))
    quit(save = "no", status = 1, runLast = FALSE)
}

cat(paste("Running contrast: ", args$contrast, "\n", sep=""))

if (exists('contrast', where=args)) {
    expContrast <- args$contrast
} else {
    cat("Contrasts: a string representing a contrast is required, please provide one using the option --contrast")
    quit(save = "no", status = 1, runLast = FALSE)
}

cat(paste("Using column: ", args$group_column, ", as grouping column\n", sep=""))

if (exists('group_column', where=args)) {
    groupVariable <- args$group_column
} else {
    cat("Group column: a name for a grouping column in targets is required, please provide one using the option e.g. --group_column=Sample\n")
    quit(save = "no", status = 1, runLast = FALSE)
}


inputDir <- ifelse(exists('input_dir', where=args), args$input_dir, "input")


library(affy, verbose=FALSE, lib.loc="/nfs/users/nfs_r/rb11/local/lib/lib-r/3.2.2")


if (file.exists(file.path(wd, inputDir))){    
    rawAffyData <- ReadAffy(filenames=targets$filename, sampleNames=targets$filename, celfile.path=inputDir, verbose=T)
} else {
    cat(paste("Read CEL files: directory '", inputDir, "' could not be found, wd=", wd, "\n", sep=""))
    quit(save = "no", status = 1, runLast = FALSE)
}

eset <- rma(rawAffyData)


##----------------------------
## EBAYES FIT
##----------------------------
cat("\nComparing and ranking...\n\n")

fac <- factor(targets[[groupVariable]])

lev <- levels(fac)

design <- model.matrix(~0+fac)

colnames(design) <- lev

fit <- lmFit(eset, design=design)

contrast.matrix <- makeContrasts(contrast=expContrast, levels=design)

contrastFit <- contrasts.fit(fit, contrast.matrix)

expFit <- eBayes(contrastFit)


##----------------------------
### ANNOTATION
##----------------------------

cat("\nObtaining annotation information...\n\n")

library(annotate, verbose=FALSE, lib.loc="/nfs/users/nfs_r/rb11/local/lib/lib-r/3.2.2")
library("hgu133plus2.db", verbose=FALSE, lib.loc="/nfs/users/nfs_r/rb11/local/lib/lib-r/3.2.2")

##probeLocations <- as.matrix(lapply(lookUp(expFeatureNames, "hgu133plus2", what="CHRLOC"), function(x) unique(paste(names(x),x, sep=":"))))
##refSeq <- as.matrix(lapply(lookUp(expFeatureNames, "hgu133plus2", what="REFSEQ"), function(x) paste0(x, collapse=",")))

expFeatureNames <- featureNames(eset)

geneSymbols <- as.matrix(getSYMBOL(expFeatureNames, "hgu133plus2")) 

ensemblID <- as.matrix(lookUp(expFeatureNames, "hgu133plus2", what="ENSEMBL"))

addedAnnotation <- cbind(geneSymbols, ensemblID)

colnames(addedAnnotation) <- c("Symbol", "Ensembl")


##-------------------------------------
## Generate Toptable and write to file
##-------------------------------------

cat("\nGenerating Toptable ")

expFit$genes <- data.frame(ID=expFeatureNames, addedAnnotation)

expToptable <- topTable(expFit, coef=colnames(contrast.matrix), n=length(expFit$F), adjust="fdr")

outputDir <- ifelse(exists('output_dir', where=args), args$output_dir, ".")

fileName <- file.path(outputDir, paste(expContrast, "_TopGenes.txt", sep = ""))

cat(paste("and saving it into file: ", fileName, "\n\n", sep=""))

write.table(as.matrix(expToptable[order(abs(expToptable$adj.P.Val), decreasing = FALSE), TRUE]),file = fileName, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAll Done.\n\n")

##-------------------------------
## QC (just performing DE is OK)
##-------------------------------
##celNames <- sampleNames(rawAffyData)
##rd <- function (x) { paste(targets$sample[match(x, targets$filename)], targets$bioRep[match(x, targets$filename)], sep="") }
##sampleNames <- unlist(lapply(celNames, rd))
##pdf(file="qc/affy_data_qc_1.1.pdf", onefile = TRUE)
##boxplot(rawAffyData, main="Raw Data", names=sampleNames)
##boxplot(as.data.frame(exprs(exp.eset)), main="RMA Normalised Data", names=sampleNames) ##for the normalised data
##plotPCA(exp.eset, addtext=sampleNames)
##dev.off()


