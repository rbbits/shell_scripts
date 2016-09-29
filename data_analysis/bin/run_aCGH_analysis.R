#=====================================================================#
# Script can be run in non-interective mode using Rscript:
#
# $ /software/R-3.0.0/bin/Rscript ./bin/run_aCGH_analysis.R --species=mouse --targets=targets.txt
#
# Or it can be run in interactive mode like this to include the arguments:
#
# $ /software/R-3.0.0/bin/R-3.0.0 --args --targets=targets.txt --species=mouse
#=====================================================================#

##------------------------------------
## load libraries
##------------------------------------
require("limma")
require("DNAcopy")


##------------------------------------
## functions
##------------------------------------
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")


##------------------------------------
## check basic working directory
##------------------------------------
workingDir <- getwd()
inputDir <- "input"
outputDir <- "output"
resourcesDir <- "resources"

if ( ! file.exists(file.path(workingDir, inputDir)) ){
    stop(paste('cannot access input directory', file.path(workingDir, inputDir),': No such directory', sep=''), call. = FALSE)
}


##-----------------------------------
## read arguments 
##-----------------------------------
targetsFile <- ""
species <- ""

if (! interactive()){
    cat('\nINFO: Running in non-interactive mode\n\n')    
    args <- commandArgs(TRUE)
    argsDF <- as.data.frame(do.call(rbind, parseArgs(args)))
    args <- as.list(as.character(argsDF$V2))
    names(args) <- argsDF$V1

    if(exists('targets', where = args)){
        targetsFile <- args$targets
    }else{
        targetsFile <- ""
    }
    
    if(! exists('species', where = args)){
        stop('--species: either mouse or human is required for argument', call.=FALSE)
    } else {
        species <- args$species
    }
} else {
    cat('\nINFO: Running in interactive mode\n\n')    
    if (species == "") {
        stop('Running in interactive mode: make sure the following variable exists if you dont like this error: species <- human | mouse\n\n', call. = FALSE)
    }    
}


##-----------------------------------
## analysis
##-----------------------------------

if (species == "mouse") {
    chroms <- 19
    xchrom <- 20
    ychrom <- 21
    geneInfo <- "resources/geneInfoMouse.rda"
} else if (species == "human") {
    chroms <- 22
    xchrom <- 23
    ychrom <- 24
    geneInfo <- "resources/geneInfoHuman.rda"
} else {
    stop('--species: either mouse or human is required for this variable', call.=FALSE)
}

if (targetsFile == ""){
    targetsFile <- file.path(workingDir, resourcesDir, targetsFile)
}

if ( ! file.exists(targetsFile) ){
    stop(paste('cannot access targets file', file.path(workingDir, targetsFile),': No such file or directory', sep=''), call. = FALSE)
}

dir.create(file.path(workingDir, outputDir), recursive = TRUE)
cat('directory', file.path(workingDir, outputDir), 'has been created\n')

## read targets
targets = readTargets(targetsFile)

cat('\nReading input files...\n\n')

## create list of files
arrayFiles <-  file.path(workingDir, inputDir, targets$FileName)

## define columns to read
listCols <- list(R="rMedianSignal", G="gMedianSignal", Rb="rBGMeanSignal", Gb="gBGMeanSignal")
listAnno <- c("Row", "Col", "ControlType", "ProbeName", "GeneName", "SystematicName", "PositionX", "PositionY", "gIsFeatNonUnifOL", "rIsFeatNonUnifOL", "gIsBGNonUnifOL", "rIsBGNonUnifOL", "gIsFeatPopnOL", "rIsFeatPopnOL", "gIsBGPopnOL", "rIsBGPopnOL", "rIsSaturated", "gIsSaturated")

## read raw data
RG <- read.maimages(arrayFiles, source="agilent", columns=listCols, annotation=listAnno, names=basename(arrayFiles))

## prepare to add new columns
## block size: 456 rows
block <- c(floor(1:nrow(RG$genes)/456)+1)
chrom <- gsub("chr([0-9XY]+):.*", "\\1", RG$genes[, "SystematicName"])
start <- gsub("chr([0-9XY]+):([0-9]+)-[0-9]+", "\\2", RG$genes[, "SystematicName"])
end <- gsub("chr([0-9XY]+):([0-9]+)-([0-9]+)", "\\3", RG$genes[, "SystematicName"])
posit <- round((as.numeric(start) + as.numeric(end))/2)
stat <- paste("Chrom", chrom, sep="")
naMe <- c(which(!chrom %in% c(1:chroms, "X", "Y")), which(RG$genes[, "ControlType"] != 0))
numSamples <- length(targets$FileName)

stat[naMe] <- ""
chrom[naMe] <- NA
posit[naMe] <- NA
chrom[chrom == "X"] <- xchrom
chrom[chrom == "Y"] <- ychrom

# A value of 1 indicates that the Cy3 channel is the reference,
# whilst a value of -1 equates to Cy5 being the reference
# NB: WTSI traditionally uses Cy3 for reference
# (As many 1/-1's as number of samples)
RG$design <- c(rep(1, times=numSamples))

## add extra columns
RG$genes$Chr = as.numeric(chrom)
RG$genes$Position = as.numeric(posit)
RG$genes$Block = as.numeric(block)
RG$genes$Column = as.numeric(RG$genes$Col)
RG$genes$Status = stat

## get layout
RG$printer = getLayout(RG$genes)
types <- readSpotTypes("resources/SpotTypes.txt")
RG$genes$Status <- controlStatus(types, RG)

cat('\nData normalisation ... \n')
## obtain weights for normalization
# First set all to 1
weights.for.norm <- rep(1,times=length(RG$genes$Chr))
# Then change empty spots to 0
weights.for.norm[is.na(RG$genes$Chr)] <- 0
# Also chromosome X and Y spots to 0
weights.for.norm[RG$genes$Chr == xchrom | RG$genes$Chr == ychrom] <- 0

## Normalise Within array
MA.norm <- normalizeWithinArrays(backgroundCorrect(RG, method = "normexp"), RG$printer, weights=weights.for.norm, method = "loess")

## Normalise Between array (necessary? seems not to hurt)
MA.norm.ba <- normalizeBetweenArrays(MA.norm, method="Aquantile")


## Prepare for CNA

# drop dark spots, controls, NA, etc (similar o the same? to naMe)
dropMe <- c(which(!chrom %in% c(1:ychrom)), which(MA.norm.ba$genes[, "ControlType"] != 0))
# map chromosome names
# 1:24 human or 1:21 mouse
map<-1:ychrom
#names(map)<-paste(c(1:chroms,"X","Y"), sep="")
names(map)<-paste(c(1:ychrom), sep="")
# genome data (M)
genomdat<-MA.norm.ba$M[-dropMe,]
# chromosome (a bit repetitive)
genomdat.chrom <- chrom[-dropMe] 
genomdat.maploc <- as.numeric(posit[-dropMe])
genomdat.chrom_order<-map[(genomdat.chrom)]

## Run CNA
cat('\n\nData segmentation ...\n\n')
cna <- CNA(genomdat, genomdat.chrom_order, genomdat.maploc, data.type = "logratio", sampleid = colnames(MA.norm.ba$M))

## smooth data
cna.smooth <- smooth.CNA(cna)

## segment data into regions of estimated equal copy
## number using CBS (circular binary segmentation) on smoothed data.Undo
## method is needed to get rid of unnecessary change-points. Below all splits
## that are not at least three SDs apart are removed. Apply a less restrictive
## standard deviation parameter for better visualisation --On a Case by Case basis-- (default = 3)
cna.smooth.segment <- segment(cna.smooth, undo.splits="sdundo", undo.SD=3, verbose=1)

# Plot whole genome
#jpeg(file=paste(workingDir, outputDir,"cna_3sd_wholegenome.jpg", sep="/"), width=1440, height=1440, pointsize=10, quality=75, bg="white", res=NA)
pdf(file=paste(workingDir, outputDir, "cna_3sd_wholegenome.pdf", sep="/"), width=30, height=30, pointsize=36)
plot(cna.smooth.segment, plot.type = "whole", pt.cols=c("green","blue"),ylim=range(-2:2), cex.lab=1)
dev.off()

# Plot by chromosome within a study
pdf(file=paste(workingDir, outputDir, "cna_3sd_samplebychrom.pdf", sep="/"), width=30, height=30, pointsize=36)
plot(cna.smooth.segment, plot.type = "samplebychrom", pt.cols=c("green","blue"), ylim=range(-2:2))
dev.off()

# Plot by chromosome across studies
pdf(file=paste(workingDir, outputDir, "cna_3sd_plateau.pdf", sep="/"), width=30, height=30, pointsize=36)
plot(cna.smooth.segment, plot.type = "plateau", pt.cols=c("green","blue"), ylim=range(-2:2))
dev.off()

# This plot focuses on the chromosome across the study's samples
pdf(file=paste(workingDir, outputDir, "cna_3sd_chrombysample.pdf", sep="/"), width=30, height=30, pointsize=36)
plot(cna.smooth.segment, plot.type = "chrombysample", pt.cols=c("green","blue"), ylim=range(-2:2))
dev.off()

# subset, for example for this study chrom by chrom
# and plot together with the w and s plots to compare
#pdf(file=paste(workingDir, outputDir, "cna_chromosome.pdf", sep = "/"), width=30, height=30, pointsize=36)
#for (c in 1:22) {
#  cna.smooth.segment.chr <- subset(cna.smooth.segment, c)
#  plot(cna.smooth.segment.chr, plot.type="w", pt.cols=c("green","blue"))
#  plot(cna.smooth.segment.chr, plot.type="s")
#}
#dev.off()

## b-Identifying Segment Gain Or Loss (SGOL)

## First, the segment list we created  is a data frame but can not be used
## directly for computation each row contains the segment data for a given segment
## within a sample. As samples usually differ in the number of altered regions and
## even for the same region where multiple samples show alterations , the margins
## of the altered matrix rarely align across samples. Use at this point, CNTools
## packages. object = "CNSeg" a reduced segment can be generated in three ways; by
## chromosomal regions that overlap across sample (by = region), by genes (by =
## gene), or by pair of samples with chromosome regions aligned (by = pair). User
## may choose to imput cells (by region or gene only) where a value can not be
## assigned by setting imput = TRUE. The X and Y chromosomes can dropped by stting
## XY = FALSE.

## For computation, the segment list needs to be converted into a matrix format with
## segments as rows and samples CNTools packages.
require(CNTools)
seg.list <- cna.smooth.segment[["output"]]
cnseg <- CNSeg(seg.list)

## Segment means can be assigned to genes  within the segments for each samples.
## I have  produced a genome file "geneInfoHuman.rda" or "geneInfoMouse.rda"
load(geneInfo)

## by gene:
## by assigning segment means to genes within the segments for each sample a reduced
## segment matrix (with genes as rows and samples as columns) can be obtained. In this case
## the number of rows of the matrix does not vary with the samples of concern
convertedData <- getRS(cnseg, by ="gene", imput=FALSE, XY=FALSE, geneMap = geneInfoHuman, what="median")
reduced_bygene <- rs(convertedData)
write.table(reduced_bygene, file=paste(workingDir, outputDir, "reduced_bygene_3sd.xls", sep="/"), sep="\t")

## by region:
## by aligning chromosome segments across samples and then assigning segment means
## to overlapping chromosomal fragments  a reduced segment matrix (with the
## fragments as rows and samples as columns) can be obtained. This approach is straight-
## forward but the drawback is the number of chromosomal fragments (or number of rows
## of the matrix) varies depending on the samples used for the calculation.
# NOTE: getRS by region fails if run on a single sample
convertedData<- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what="median")
reduced_byregion=rs(convertedData)
write.table(reduced_byregion, file=paste(workingDir, outputDir, "reduced_region_3sd.xls", sep="/"), sep="\t")

## Working on the data converted from segment list, we can compute the SGOL scores
## for genes (or chromosomal fragment if samples are aligned by regions) by calculating the
## summations (parameter method) for all the positive values over a set threshold and all
## the negative values below a set threshold (threshold below).
require(cghMCR, quietly = TRUE)

SGOLScores <- SGOL(convertedData, threshold = c(-0.3, 0.3), method = sum)

pdf(file=paste(workingDir, outputDir, "sgol_scores.pdf", sep="/"), width=30, height=30, pointsize=36)
plot(SGOLScores)
dev.off()

## Based on the SGOL scores, genes in regions of gains or losses can be obtained by a
## set of thresholds, say -1.5 and 1.5
GOIGains <- SGOLScores[which(as.numeric(unlist(gol(SGOLScores[, "gains"]))) > 1.5), "gains"]
GOILosses <- SGOLScores[which(as.numeric(unlist(gol(SGOLScores[, "losses"]))) < -1.5), "losses"]

write.table(gol(GOIGains), file=paste(workingDir, outputDir, "list_of_gains.txt", sep="/"), sep="\t")
write.table(gol(GOILosses), file=paste(workingDir, outputDir, "list_of_losses.txt", sep="/"), sep="\t")

## c. Identifying Minimum Common Regions (MCR)

## Only segements with means less or greater than the lower or upper threshold values
## will be considered as altered regions and included in the subsequent analysis
cghmcr <- cghMCR(cna.smooth.segment, gapAllowed = 500, alteredLow = 0.9, alteredHigh = 0.9, recurrence = 75)
mcrs <- MCR(cghmcr)

## adding the number of samples in the region
list.samples <- strsplit(mcrs[,9],",")
num.samples <- summary(list.samples)
num.samples <- as.numeric(num.samples[,"Length"])
mcrs.bind <-cbind (mcrs[, c ("chromosome", "status", "loc.start", "loc.end","mcr.start", "mcr.end","samples")], as.numeric(num.samples))
write.table(mcrs.bind, file=paste(workingDir, outputDir, "mcrs.xls", sep="/"), sep="\t", row.names=FALSE)
