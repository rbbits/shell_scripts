#=====================================================================#
# Script can be run in non-interective mode using Rscript:
#
# $ /software/R-3.0.0/bin/Rscript ./bin/run_snapCGH_analysis.R --species=mouse --targets=resources/targets.txt
#
# Via LSF:
#
# $ bsub -J "snapCGH_analysis" -o log/snapcgh_analysis.%J.o -e log/snapcgh_analysis.%J.e -R 'select[mem>12000] rusage[mem=12000]' -M 12000 "/software/R-3.0.0/bin/Rscript ./bin/run_snapCGH_analysis.R --targets=resources/targets.txt --species=mouse"
#
# Or it can be run in interactive mode like this to include the arguments:
#
# $ /software/R-3.0.0/bin/R-3.0.0 --args --targets=targets.txt --species=mouse
#
# Then run inside R using source("bin/run_snapCGH_analysis.R")
#=====================================================================#

##------------------------------------
## load libraries
##------------------------------------
library("snapCGH")


##------------------------------------
## functions
##------------------------------------
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")


##------------------------------------
## check basic working directory
##------------------------------------
workingDir <- getwd()
inputDir <- "input"
outputDir <- "output/snapcgh"
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

if (species == "mouse") {
    chroms <- 19
    xchrom <- 20
    ychrom <- 21
} else if (species == "human") {
    chroms <- 22
    xchrom <- 23
    ychrom <- 24
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


##-----------------------------------
## analysis
##-----------------------------------

## read targets
targets = readTargets(targetsFile)

cat('\nReading input files...\n\n')
RG = read.maimages( file.path(workingDir, inputDir, targets$FileName), source = "agilent")

block <- c(floor(1:nrow(RG$genes)/456)+1)
chrom <- gsub("chr([0-9XY]+):.*", "\\1", RG$genes[, "SystematicName"])
start <- gsub("chr([0-9XY]+):([0-9]+)-[0-9]+", "\\2", RG$genes[, "SystematicName"])
end <- gsub("chr([0-9XY]+):([0-9]+)-([0-9]+)", "\\3", RG$genes[, "SystematicName"])
posit <- round((as.numeric(start) + as.numeric(end))/2000000) # middle point in Mb (values in Mb required by snapCGH library)
stat <- paste("Chrom", chrom, sep="")
naMe <- c(which(!chrom %in% c(1:chroms, "X", "Y")), which(RG$genes[, "ControlType"] != 0))
numSamples <- length(targets$FileName)

stat[naMe] <- ""
chrom[naMe] <- NA
posit[naMe] <- NA
chrom[chrom == "X"] <- xchrom
chrom[chrom == "Y"] <- ychrom
#if necessary: chrom[chrom == "M"] <- ?

RG$genes$Chr = as.numeric(chrom)
RG$genes$Position = as.numeric(posit)
RG$genes$Block = as.numeric(block)
RG$genes$Column = as.numeric(RG$genes$Col)
RG$genes$Status = stat

RG$printer <- getLayout(RG$genes)
types <- readSpotTypes("resources/SpotTypes.txt")
RG$genes$Status <- controlStatus(types, RG)

# A value of 1 indicates that the Cy3 channel is the reference,
# whilst a value of -1 equates to Cy5 being the reference
# NB: WTSI traditionally uses Cy3 for reference
# (As many 1/-1's as number of samples)
RG$design <- c(rep(1, times=numSamples))

cat('\nBackground correction ...\n')
# all background correction methods are tested but only one is required
# choose the one that works the best for your own data
RG.nn <- backgroundCorrect(RG, method = "none")
RG.sb <- backgroundCorrect(RG, method = "subtract")
RG.mb <- backgroundCorrect(RG, method = "minimum")
RG.nb <- backgroundCorrect(RG, method = "normexp") # <-- most cases use this


cat('generate background correction comparison plot...\n')
# create a pdf to compare methods. JPEG files are faster to display
# pdf(file=paste(workingDir, outputDir, "rg_compare_bkgnd.pdf", sep="/"))
jpeg(file=paste(workingDir, outputDir, "rg_compare_bkgnd.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
par(mfrow = c(2, 2))
plotMA(RG.nn, array = 1, xlim = c(0, 16), ylim = c(-5, 5), main = "bkgcorr: none")
plotMA(RG.mb, array = 1, xlim = c(0, 16), ylim = c(-5, 5), main = "bkgcorr: minimum")
plotMA(RG.sb, array = 1, xlim = c(0, 16), ylim = c(-5, 5), main = "bkgcorr: subtract")
plotMA(RG.nb, array = 1, xlim = c(0, 16), ylim = c(-5, 5), main = "bkgcorr: normexp")
dev.off()

cat('\nData normalisation ... \n')
# Q = Which background method to choose? And then which normalization method?
# A = use norm=loess if bg=minimum|subtract|normexp ;; use norm=median if bg=normexp
#     Method "printtiploess fails due to discrepancies bewtween the M and the print layout.
#     However, it is not relevant anymore for Agilent arrays because of their design.
#     Same applies to method "robustspline", see: https://stat.ethz.ch/pipermail/bioconductor/2012-November/049520.html
# MA.p = normalizeWithinArrays(RG.mb, method="printtiploess")
# MA.c = normalizeWithinArrays(RG.nb, method = "composite")
# no norm (any RG.xx will do):
MA <- MA.RG(RG.mb)
# alternatives for loess:
MA.nb.l <- normalizeWithinArrays(RG.nb, method = "loess")
MA.mb.l <- normalizeWithinArrays(RG.mb, method = "loess")
MA.sb.l <- normalizeWithinArrays(RG.sb, method = "loess")
# alternatives for median:
MA.nb.m <- normalizeWithinArrays(RG.nb, method = "median")
MA.mb.m <- normalizeWithinArrays(RG.mb, method = "median")

cat('\ngenerate normalisation methods comparison plot ...\n')
# create a pdf to compare methods. JPEG files are faster to display
# pdf(file=paste(workingDir, outputDir, "ma_compare_norm.pdf", sep="/"))
jpeg(file=paste(workingDir, outputDir, "ma_compare_norm.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
par(mfrow = c(3, 2))
plotMA(MA, array = 1, main = "no norm")
plotMA(MA.nb.l, array = 1, main = "loess norm with normexp bg corr")
plotMA(MA.mb.l, array = 1, main = "loess norm with min bg corr")
plotMA(MA.sb.l, array = 1, main = "loess norm with subst bg corr")
plotMA(MA.nb.m, array = 1, main = "median norm with normexp bg corr")
plotMA(MA.mb.m, array = 1, main = "median norm with min bg corr")
dev.off()

# TODO: change code to ask for user input to select the MA to be used
#       otherwise the code has to be run once, change the code, then run again
#       Something like this:
#       MA <- readLine("Look at the MA plots and pick the best")
#       -- some code to choose MA from user input best goes here --
#       MA2 <- processCGH(MA, method.of.averaging = mean, ID = "ProbeName")

cat('\nProcess CGH data ...\n')
# Use the processCGH to 'tidy up' the MAList object [1]
MA2 <- processCGH(MA.sb.l, method.of.averaging = mean, ID = "ProbeName")


#--------------#
# MA Plots     #
#--------------#

cat('\nGenerate MA plots ...\n\n')
cat(' > whole genome ...\n')
#pdf(file = paste(workingDir, outputDir, "ma_genes_wholegenome.jpg", sep = "/"))
jpeg(file=paste(workingDir, outputDir,"ma_wholegenome.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
genomePlot(MA2)
dev.off()

# plot chrom by chrom, alternative view to previous but we end up with 22-24 plots instead of 1
cat(' > per chromosome ...\n')
pdf(paste(workingDir, outputDir, "ma_perchromosome.pdf", sep = "/"))
for (c in 1:ychrom) {
    cat('   chromosome',c)
    genomePlot(MA2, chrom.to.plot = c)
}
dev.off()


#--------------#
# Segmentation #
#--------------#

library("GLAD")

cat('\n\nData segmentation ...\n\n')
# Fit the segmentation method. There are 4 methods for calling segmentation algorithms:
# HomMMH, DNAcopy, GLAD and tilingArray:
SegInfo.GLAD = runGLAD(MA2)
SegInfo.DNACopy = runDNAcopy(MA2, undo.splits = "sdundo", undo.SD = 1.5)
SegInfo.Hom <- runHomHMM(MA2, criteria = "AIC")
# SegInfo.Bio <- runBioHMM(MA2)              # runs ok BUT takes ages!
# SegInfo.TilingArray <- runTilingArray(MA2) # Fails, and requires a lot of memory (+8Gb)

## Segmentation methods sometimes have a tendency to fit states
## whose means are very close together. We overcome this problem
## by merging states whose means are within a given threshold.
# There are two different methods for carrying out the merging process, use default (1)
SegInfo.GLAD.merged <- mergeStates(SegInfo.GLAD, MergeType = 1)
SegInfo.DNACopy.merged <- mergeStates(SegInfo.DNACopy, MergeType = 1)
SegInfo.Hom.merged <- mergeStates(SegInfo.Hom, MergeType = 1)
# SegInfo.Bio.merged <- mergeStates(SegInfo.Bio, MergeType = 1)  # runs ok but read comment above


#--------------#
# More Plots   #
#--------------#

cat('\nGenerate Segmentation plots ...\n')
cat('\n > whole genome ...\n')
#pdf(file = paste(workingDir, outputDir, "segments_wholegenome_homhmm.pdf", sep = "/"))
jpeg(file=paste(workingDir, outputDir,"segments_wholegenome_homhmm.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
plotSegmentedGenome(SegInfo.Hom.merged, array = 1)
dev.off()

# pdf(file = paste(workingDir, outputDir, "segments_wholegenome_dnacopy.pdf", sep = "/"))
jpeg(file=paste(workingDir, outputDir,"segments_wholegenome_dnacopy.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
plotSegmentedGenome(SegInfo.DNACopy.merged, array = 1)
dev.off()

# pdf(file = paste(workingDir, outputDir, "segments_wholegenome_glad.pdf", sep = "/"))
jpeg(file=paste(workingDir, outputDir,"segments_wholegenome_glad.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
plotSegmentedGenome(SegInfo.GLAD.merged, array = 1)
dev.off()

# pdf(file = paste(workingDir, outputDir, "segments_wholegenome_biohmm.pdf", sep = "/"))
# jpeg(file=paste(workingDir, outputDir,"segments_wholegenome_biohmm.jpg", sep="/"), width=6500, height=6500, pointsize=102, quality=8, bg="white", res=NA)
# plotSegmentedGenome(SegInfo.Bio.merged, array = 1)
# dev.off()

cat('\n > per chromosome ...\n')
# equally, we run chrom by chrom ...
pdf(paste(workingDir, outputDir, "segments_perchromosome_homhmm.pdf", sep = "/"))
for (c in 1:ychrom) {
  plotSegmentedGenome(SegInfo.Hom.merged, chrom.to.plot = c)
}
dev.off()

pdf(paste(workingDir, outputDir, "segments_perchromosome_dnacopy.pdf", sep = "/"))
for (c in 1:ychrom) {
  plotSegmentedGenome(SegInfo.DNACopy.merged, chrom.to.plot = c)
}
dev.off()

pdf(paste(workingDir, outputDir, "segments_perchromosome_glad.pdf", sep = "/"))
for (c in 1:ychrom) {
  plotSegmentedGenome(SegInfo.GLAD.merged, chrom.to.plot = c)
}
dev.off()

cat('\n > combined segmentation methods ...')
pdf(paste(workingDir, outputDir, "combined_segments_perchrom_homhmm_dnacopy.pdf", sep = "/"))
for (c in 1:ychrom) {
    plotSegmentedGenome(SegInfo.DNACopy.merged, SegInfo.Hom.merged, array = 1, chrom.to.plot = c, colors = c("blue", "red"))
}
dev.off()

pdf(paste(workingDir, outputDir, "combined_segments_perchrom_dnacopy_glad.pdf", sep = "/"))
for (c in 1:ychrom) {
    plotSegmentedGenome(SegInfo.DNACopy.merged, SegInfo.GLAD.merged, chrom.to.plot = c)
}
dev.off()

pdf(paste(workingDir, outputDir, "combined_segments_perchrom_homhmm_glad.pdf", sep = "/"))
for (c in 1:ychrom) {
    plotSegmentedGenome(SegInfo.Hom.merged, SegInfo.GLAD.merged, array = 1, chrom.to.plot = c, colors = c("blue", "red"))
}
dev.off()

cat('\n\nAll Done.\n\n')
