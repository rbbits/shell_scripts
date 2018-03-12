##-----------------------------------
## Some info
##-----------------------------------
# A stepped process controls the memory-hungry process of producing
# SummarizedExperiment objects by doing it a single bam file at a time.
#
# Step 1 creates SE objects and serialise them to save them as Rds files.
# Arguments: --ibam=<./deseq2/input/sample.bam> --bamid=<sampleId | runId> --step=1 [--merged] [--sedata=$SERdsfile]
#
# Step 2 reads these files and merges them into a singe SE object to be used by the DESeq2 function
# Arguments: --irds=<./deseq2/input> --contrasts=<./resources/deseq2_contrasts.txt> --targets=$targets --pooling=<percontrast | allsamples> --step=2 [--sedata=<./deseq2/input/se.Rds>]
# Example: args <- c("--irds=./deseq2/input", "--contrasts=./resources/deseq2_contrasts.txt", "--targets=./resources/targets_deseq2.txt", "--pooling=percontrast", "--step=2", "--varofinterest=genotype", "--merged")
# Example: args <- c("--irds=./deseq2/input", "--contrasts=./resources/deseq2_contrasts.txt", "--targets=./resources/targets_deseq2.txt", "--pooling=percontrast", "--step=2", "--varofinterest=genotype")
#
# No Steps runs the whole process in one go: args values are extracted from config.yml file
# Arguments: --idir=<./deseq2/input> --contrasts=<./resources/deseq2_contrasts.txt> --targets=<./resources/targets_deseq2.txt> --pooling=<percontrast | allsamples> --step=0 [--merged] [--sedata=<./deseq2/input/se.Rds>]
# Example: args<-c("--idir=./deseq2/input", "--contrast=./resources/deseq2_contrasts.txt", "--targets=./resources/targets_deseq2.txt", "--pooling=allsamples", "--merged=0", "--sedata=./deseq2/input/summarizedExperiment_2014-12-03_1537.Rds")
#
# R can be run in interactive mode like this to include the arguments (change accordingly):
#
# R-3.2.2 --args --irds=deseq2/input --contrasts=resources/deseq2_contrasts.txt --targets=resources/targets_deseq2.txt --pooling=allsamples --step=2 --varofinterest=genotype --species=mouse --odir=deseq2/output


##---------------------------------------------------------------------------------
## Read arguments and run the wrapper function for DESeq2
##---------------------------------------------------------------------------------
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

args <- commandArgs(TRUE)

argsDF <- as.data.frame(do.call(rbind, parseArgs(args)))

args <- as.list(as.character(argsDF$V2))

names(args) <- argsDF$V1

library("BiocParallel", quietly=TRUE, verbose=FALSE)

library("tools", quietly=TRUE, verbose=FALSE)

register(MulticoreParam(4))

step <- as.numeric(args$step)

merged <- exists('merged', where=args)

wd <- getwd()

formatDate <- format(Sys.time(), format="%Y-%m-%d")

if (exists('experimentid', where=args)) {
    expId <- args$experimentid
} else {
    expId <- "NoExpId"
}
    

##-----------------------------------
## DESeq2 Functions
##-----------------------------------
if (file.exists("bin/run_deseq2_functions.R")){
    source("bin/run_deseq2_functions.R")
} else {
    cat(paste("source: file bin/run_deseq2_funtions.R could not be found in current directory: ", cwd, "\n", sep=""))
    quit(save = "no", status = 1, runLast = FALSE)
}


library("DESeq2", quietly=TRUE, verbose=FALSE)

library("Rsamtools", quietly=TRUE, verbose=FALSE)

library("GenomicAlignments", quietly=TRUE, verbose=FALSE)

# If we are running the script several times, we can save time by saving the main
# SerializedExperiment object and pass the /path/to/rdsfile in list of arguments

if(!exists('sedata', where=args)) {

    if ( !step || step == 1) {

        if (exists('reposdir', where=args)) {
            reposdir <- args$reposdir
        } else {
            reposdir <- "/lustre/scratch117/core/sciops_repository/transcriptomes"
        }

        cat(paste("Using repository directory:", reposdir, "\n", sep=""))

        if (exists('datadir', where=args)) {
            datadir <- args$datadir
        } else {
            datadir <- "/nfs/gapi/users/rb11/data"
        }

        cat(paste("Using data directory: ", datadir, "\n", sep=""))
        
        if(args$species == "mouse") {
            
            transcriptsSpecies="Mus musculus"

            if (exists('genome', where=args)) {
                genomeVersion <- args$genome
            } else {
                genomeVersion <- "GRCm38"
            }

            if (exists('transcriptome', where=args)) {
                transcriptomeVersion <- args$transcriptome
            } else {
                transcriptomeVersion <- "ensembl_75_transcriptome"
            }

        } else if (args$species == "human") {

            transcriptsSpecies <- "Homo sapiens"

            if (exists('genome', where=args)) {
                genomeVersion <- args$genome
            } else {
                genomeVersion <- "1000Genomes_hs37d5"
            }

            if (exists('transcriptome', where=args)) {
                transcriptomeVersion <- args$transcriptome
            } else {
                transcriptomeVersion <- "ensembl_75_transcriptome"
            }
            
        }
        
        library("GenomicFeatures", quietly=TRUE, verbose=FALSE)
        
        transcriptsBasename <- paste(transcriptomeVersion, genomeVersion, sep="-")
        
        transcriptsDB <- paste(transcriptsBasename, "sqlite", sep=".")
        
        transcriptsGFF <- paste(transcriptsBasename, "gtf", sep=".")
        
        chromInfo <- paste(genomeVersion, "chrominfo", sep=".")
        
        gtfFile <- paste(reposdir, gsub(" ", "_", transcriptsSpecies), transcriptomeVersion, genomeVersion, "gtf", transcriptsGFF, sep="/")
        
        chrominfo <- read.table(paste(datadir, chromInfo, sep="/"), col.names=c("chrom", "length", "is_circular"))
            
        # Read Gene Model from GTF. 
        if(file.exists(paste("resources/", transcriptsDB, sep=""))) {
            txdb <- loadDb(paste("resources/", transcriptsDB, sep=""))
        } else if(file.exists(paste(datadir, transcriptsDB, sep="/"))) {
            txdb <- loadDb(paste(datadir, transcriptsDB, sep="/"))
        } else {
           
            # Read Gene Model from GTF
            txdb <- makeTxDbFromGFF(gtfFile,
                                    format="gtf",
                                    chrominfo=chrominfo,
                                    organism=transcriptsSpecies)
                                     
            if (file.access(datadir, mode = 2)) {
                saveDb(txdb, file=paste(datadir, transcriptsDB, sep="/"))
            } else {
                saveDb(txdb, file=paste(".", "resources", transcriptsDB, sep="/"))
            }

        }
    
        # List of exons grouped by gene
        exonsByGene <- exonsBy(txdb, by="gene")
    }

    if (!step) {  # NO STEPS

        if (args$pooling == "allsamples") {

            # Read targets and select samples to process
            samplesData <- read.table(args$targets, comment.char="", sep="\t", header=TRUE)

            # get correct run ids for colnames(se)
            bamfilenames <- paste(samplesData$RunId, "accepted_hits.bam", sep="_")
            
            # list of files using ids listed in targets
            fls <- paste(args$idir, bamfilenames, sep="/")

            # format runid names based on bamids
            runids <- gsub("#", "_", samplesData$RunId)

            # Define fls list as BAM files
            bamLst <- BamFileList(fls, yieldSize=100000, asMates=TRUE)
        
            # create counts object
            se <- summarizeOverlaps(exonsByGene, bamLst, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

            # match the colnames with the vector file names
            colnameIdx <- match(colnames(se), bamfilenames)

            # assign the runids that correspond to the filenames in colnames
            colnames(se) <- runids[colnameIdx]
            
            # create input directory
            dir.create(args$idir)
        
            # serialize object and save it
            saveRDS(se, file = paste(args$idir, "/", expId, "_", formatDate, "_summarizedExperiment_nosteps_allsamples.Rds", sep=""))
            
            runDESeq2(se, args$targets, args$contrasts, args$varofinterest, args$pooling, merged, args$species, args$odir)

        } else if (args$pooling == "percontrast") {

            # read file of contrasts
            contrastsData <- read.table(args$contrasts, sep="\t", comment.char="", header=FALSE)

            # create a list of contrasts
            contrastsList <- as.list(as.character(contrastsData$V1))

            # Read targets and select samples to process
            samplesData <- read.table(args$targets, comment.char="", sep="\t", header=TRUE)

            # create a list of SE objects one for each contrast
            seList <- lapply(contrastsList, FUN = function(x) createSumExp (x, samplesData, args$varofinterest, merged, "bam", args$idir, step, exonsByGene))

            # each SE-contrast object will be processed by DESeq2
            invisible(lapply(seList, FUN = function(x) runDESeq2 (x, args$targets, args$contrasts, args$varofinterest, args$pooling, merged, args$species, args$odir)))

        } else {

            cat(paste("Wrong value for pooling: '", args$pooling, "' for option 'nosteps'.", sep=""))
            
            quit(save = "no", status = 1, runLast = FALSE)
    
        }
        
    } else if (step == 1) {

        .Last <- function() {
           
           cat("Finished running step1.\n")
           
           cat("Run again with --step=2 to continue the analysis.\n")
           
        }

        # a single bam file passed as argumet
        fl <- args$ibam

        bamid <- sub("#", "_", args$bamid)
        
        bamFl <- BamFile(fl, yieldSize=100000, asMates=TRUE)
        
        # se for a single file
        se <- summarizeOverlaps(exonsByGene, bamFl, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

        colnames(se) <- bamid
        
        # save rds for se object
        bamFileName <- basename(file_path_sans_ext(fl))
        
        # create input directory
        dir.create(file.path(wd, args$irds))

        saveRDS(se, file = paste(args$irds, "/", bamFileName, ".Rds", sep="")) #serialize object

        # stop the script, run again with --step=2 to continue
        quit(save = "no", status = 0, runLast = TRUE)

    } else if (step == 2) {

        # Find the RDS files created in step 1
        fls <- list.files(args$irds, pattern="accepted_hits.Rds$", full=TRUE )

        # verify that all of the Rds exist
        if (!all(file.exists(fls))){
            cat(paste("Not all required RDS files were found in ", args$irds, sep=""))
            quit(save = "no", status = 1, runLast = FALSE)
        }

        # In per-contrast DE analysis each contrast is done considering
        # the component samples independently of other contrasts', 
        # as opposed to pooling 'allsamples'
        if (args$pooling == "percontrast") {

            # read file of contrasts
            contrastsData <- read.table(args$contrasts, sep="\t", comment.char="", header=FALSE)

            # create a list of contrasts
            contrastsList <- as.list(as.character(contrastsData$V1))

            # Read targets and select samples to process
            samplesData <- read.table(args$targets, sep="\t", comment.char="", header=TRUE)

            # create a list of SE objects one for each contrast
            seList <- lapply(contrastsList, FUN = function(x) createSumExp (x, samplesData, args$varofinterest, merged, "rds", args$irds, step))

            invisible(lapply(seList, FUN = function(x) runDESeq2 (x, args$targets, args$contrasts, args$varofinterest, args$pooling, merged, args$species, args$odir)))
            
        } else if (args$pooling == "allsamples") {
            
            # Read the RDS files into a list of SE's
            seList <- lapply(fls, readRDS)

            # Combine them into a single SE object containing
            # an assay of reads from all of the bams
            se <- do.call(cbind, seList)

            # create input directory
            dir.create(file.path(wd, args$irds))

            # serialize object
            saveRDS(se, file = paste(args$irds, "/", expId, "_", formatDate, "_summarizedExperiment_step2_allsamples.Rds", sep=""))

            # run DESeq2
            runDESeq2(se, args$targets, args$contrasts, args$varofinterest, args$pooling, merged, args$species, args$odir)
            
        } else {

            cat(paste("Wrong value for pooling: '", args$pooling, "' for option 'step2'.", sep=""))
            quit(save = "no", status = 1, runLast = FALSE)

        }
        
    }
        
} else if (file.exists(args$sedata)) {
    
    se <- readRDS(file=args$sedata)
    
    # run DESeq2
    runDESeq2(se, args$targets, args$contrasts, args$varofinterest, args$pooling, merged,  args$species, args$odir)
    
}






