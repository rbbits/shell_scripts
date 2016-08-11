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

register(MulticoreParam(4))

step <- as.numeric(args$step)

merged <- exists('merged', where=args)

wd <- getwd()


##-----------------------------------
## DESeq2 Functions
##-----------------------------------
if (file.exists("bin/run_deseq2_functions.R")){
    source("bin/run_deseq2_functions.R")
} else {
    cat(paste("source: file bin/run_deseq2_funtions.R could not be found in current directory: ", cwd, sep=""))
    quit(save = "no", status = 1, runLast = FALSE)
}


##-----------------------------------
## load libraries
##-----------------------------------

library("DESeq2")
library("GenomicAlignments")
library("Rsamtools")
library("BiocParallel")
library("tools")
library("biomaRt")


# If we are running the script several times, we can save time by saving the main
# SerializedExperiment object and pass the /path/to/rdsfile in list of arguments

if(!exists('sedata', where=args)) {

    if ( !step || step == 1) {

        if(args$species == "mouse") {
            
            transcriptsDB="transcriptsEns75GRCm38gtf.sqlite"

            transcriptsGFF="/lustre/scratch110/srpipe/transcriptomes/Mus_musculus/ensembl_75_transcriptome/GRCm38/gtf/ensembl_75_transcriptome-GRCm38.gtf"

            transcriptSpecies="Mus Musculus"

            chrominfo <- read.table("/nfs/gapi/users/rb11/data/chrominfo_Ens75GRCm38.txt", col.names=c("chrom", "length", "is_circular"))
            
        } else if (args$species == "human") {
            
            transcriptsDB="transcriptsEns751000Genomes_hs37d5.sqlite"

            transcriptsGFF="/lustre/scratch110/srpipe/transcriptomes/Homo_sapiens/ensembl_75_transcriptome/1000Genomes_hs37d5/gtf/ensembl_75_transcriptome-1000Genomes_hs37d5.gtf"

            transcriptSpecies="Homo Sapiens"

            chrominfo <- read.table("/nfs/gapi/users/rb11/data/chrominfo_Ens751000Genomes_hs37d5.txt", col.names=c("chrom", "length", "is_circular"))
            
        }
        
        library("GenomicFeatures")
        
        # Read Gene Model from GTF. 
        if(file.exists(paste("resources/", transcriptsDB, sep=""))) {
            
            hse <- loadDb(paste("resources/", transcriptsDB, sep=""))
            
        } else {
           
            # Read Gene Model from GTF
            hse <- makeTxDbFromGFF(transcriptsGFF,
                                   format="gtf",
                                   exonRankAttributeName="exon_number",
                                   chrominfo=chrominfo,
                                   species=transcriptSpecies)
                                     
            saveDb(hse, file=paste("/nfs/gapi/users/rb11/data/", transcriptsDB, sep=""))

        }
    
        # List of exons grouped by gene
        exonsByGene <- exonsBy(hse, by="gene")
    }

    if (!step) {  # NO STEPS
        
        if (args$pooling == "allsamples") {
            
            # list of files. Specify BAMs location
            fls <- list.files(args$idir, pattern="bam$", full=TRUE )
        
            # get correct run ids for colnames(se)
            bamfilenames <- basename(list_files_with_exts(args$idir, "bam"))
            
            runids <- gsub("#", "_", gsub("_accepted_hits.bam", "", bamfilenames))
        
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
            saveRDS(se, file = paste(args$idir, "/summarizedExperiment_nosteps_allsamples.Rds", sep=""))
            
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

        saveRDS(se, file = paste(wd, "/", args$irds, "/", bamFileName, ".Rds", sep="")) #serialize object

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
            dir.create(file.path(wd,args$irds))

            #serialize object. Save for debugging reasons.
            saveRDS(se, file = paste(wd, "/", args$irds, "/summarizedExperiment_step2_allsamples_", format(Sys.time(), format="%Y-%m-%d"), ".Rds", sep=""))

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






