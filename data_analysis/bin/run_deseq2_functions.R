##-----------------------------------
## Functions
##-----------------------------------
runDESeq2 <- function(sumExp, targets, contrast, varofinterest, pooling, merged, ddsSpecies, outputdir) {
    
    # Add gene names
    addAnnotation <- function(presentResults, annoSpecies) {
        
        presentResults$ensembl <- sapply(strsplit(rownames(presentResults), split="\\+"), "[", 1)
            
        if (annoSpecies == "human") {
            ensDataset="hsapiens_gene_ensembl"
                
            ensGeneSymbol="hgnc_symbol"
            
        } else {
                    
            ensDataset="mmusculus_gene_ensembl"
            
            ensGeneSymbol="mgi_symbol"
            
        }            

        library("biomaRt")
    
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset=ensDataset, host="www.ensembl.org")
            
        genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene", ensGeneSymbol),
                         filters = "ensembl_gene_id",
                         values = presentResults$ensembl,
                         mart = ensembl)
                     
        idx <- match(presentResults$ensembl, genemap$ensembl_gene_id)
        
        presentResults$entrez <- genemap$entrezgene[idx]
        
        presentResults$gene_symbol <- genemap[[ensGeneSymbol]][idx]
        
        return(presentResults)
            
    }
    
    # Print the topTable
    createTopTable <- function(presentDDS, contrast, variable, fileSuffix, topTableSpecies, oDir) {

        presentContrast <- unlist(contrast)
        
        # put the "control" or "untreated" level as the first element, so that the
        # log2 fold changes and results will be most easily interpretable.
        firstLevel <- gsub("(.*)-(.*)", "\\2", presentContrast)
        
        presentDDS[[variable]] <- relevel(presentDDS[[variable]], firstLevel)
        
        # Run the DESeq pipeline and extract the results
        writeLines(paste("Running DESeq2 for contrast: ", presentContrast, sep=""))
        
        dds <- DESeq(presentDDS)
        
        # define numerator and denominator
        smplDenom <- firstLevel
        
        smplNumer <- gsub("(.*)-(.*)", "\\1", presentContrast)
        
        # extract results
        res <- results(dds, contrast=c(variable, smplNumer, smplDenom), pAdjustMethod = "BH", independentFiltering = FALSE, alpha=0.05)
        
        # add the annotation
        res.anno <- addAnnotation(res, topTableSpecies)
        
        # sort by P.Adj.Value
        res.anno <- res.anno[order(res.anno$padj),]

        # create output directory
        dir.create(oDir)

        ## output file name
        filename <- paste(oDir, "/results_", presentContrast, "_", fileSuffix, sep="")
            
        ## Write toptable
        write.csv(as.data.frame(res.anno), file=paste(filename, "csv", sep="."))

        ## Save results object
        saveRDS(res.anno, file = paste(filename, "Rds", sep="."))

    }

    # Read targets and select samples to process
    targetsAll <- read.table(targets, sep="\t", comment.char="", header=TRUE)

    # Form the run id string
    runIds <- paste(targetsAll$idrun, targetsAll$position, targetsAll$tag, sep="_")

    # Deal with non-multiplexed samples
    runIds <- gsub("(.*)(_$)", "\\1", runIds)
    
    targetsAll$run <- factor(runIds, levels=unique(runIds))

    columnName <- ifelse(merged, "sampleName", "run")
    
    # Use the right number of rows from targets file and define the formula for each case
    if(pooling == "percontrast") {
        
        # Only the samples in the contrast are considered
        samplesIdx <- match(colnames(sumExp), targetsAll[ ,columnName])

        sampleInfoCols <- c("sampleName", "run", varofinterest)

        designFmla <- as.formula(paste("~ ", varofinterest, sep=""))

        # extract contrast string smuggled by createSE function 
        # (as a list so it can be processed with the same code)
        contrastList <- metadata(sumExp)
       
    } else if(pooling == "allsamples") {
        
        # all samples pooled
        samplesIdx <- c(1:nrow(targetsAll))

        # 2016-01-13: bioRep may not be needed at all, but this 
        #             if is to provide backward compatibility
        #             (cf. comment below about new version of DESeq2)
        if (exists("bioRep", where=targetsAll)) {
            
            targetsAll$bioRep <- factor(targetsAll$bioRep)

            sampleInfoCols <- c("sampleName", "bioRep", "run", varofinterest)
            
        } else {
            
            sampleInfoCols <- c("sampleName", "run", varofinterest)
            
        }
        
        # 2015-11-11: In light of DESeq2's latest documentation 
        #             Use same formula as percontrast until we figure this out for sure
        #designFmla <- as.formula(paste("~ bioRep + ", varofinterest, sep=""))
        designFmla <- as.formula(paste("~ ", varofinterest, sep=""))
        
        # for the next line to work 'contrast' must contain the path to a file listing contrasts
        contrastsData <- read.table(contrast, sep="\t", comment.char="", header=FALSE)
        
        contrastList <- as.list(as.character(contrastsData$V1))
        
    }

    # Complete the information about the samples using targets
    # keeping only the samples needed for current analysis (if percontrast)
    sampleInfo <- DataFrame(droplevels(targetsAll[samplesIdx, sampleInfoCols]))
    
    seIdx <- match(colnames(sumExp), sampleInfo[ ,columnName])

    colData(sumExp) <- cbind(colData(sumExp), sampleInfo[seIdx, ])
    
    # Create a DESeq2-ready data structure. Collapse technical reps
    ddsFull <- DESeqDataSet(sumExp, design = designFmla)

    # Group technical replicates by sampleName
    ddsCollapsed <- collapseReplicates(ddsFull, groupby=ddsFull$sampleName, run=ddsFull$run)

    # run DESEq2 on each contrast and print topTable
    invisible(lapply(contrastList, FUN = function(x) createTopTable(ddsCollapsed, x, varofinterest, pooling, ddsSpecies, outputdir)))
    
}


# Create a list of SummarizedExperimet objects based on a list of contrasts
createSumExp <- function(contrast, samples, varofinterest, merged, filetype, inputdir, stepped, exonsbygene=NULL) {

    # unlist the current contrast
    currentContrast <- as.character(unlist(contrast))

    # select the samples that belong only to that contrast
    samplesIdx <- which(samples[[varofinterest]] %in% unique(unlist(strsplit(as.character(gsub("-", " ", currentContrast)), split=" +"))))
    
    # choose the right name for the file
    if (merged) {
        
        bamids <- as.vector(samples$sampleName[samplesIdx])
        
    } else {
        
        bamids <- paste(samples[samplesIdx,"idrun"], "_", samples[samplesIdx,"position"], "#", samples[samplesIdx,"tag"], sep="")
        
        # Deal with non-multiplexed samples
        bamids <- gsub("(.*)(#$)", "\\1", bamids)
        
    }

    # rds for stepped ; bam for no-steps
    if (filetype == "rds") {
        
        # vector of file names to be read
        contrastFiles <- paste(inputdir, "/", bamids, "_accepted_hits.Rds", sep="")
    
        # read the Rds files into a list
        contrastSElist <- lapply(contrastFiles, readRDS)
    
        # bind them into a single SE object
        contrastSE <- do.call(cbind, contrastSElist)

    } else if (filetype == "bam") {

        # vector of file names to be read
        contrastFiles <- paste(inputdir, "/", bamids, "_accepted_hits.bam", sep="")
       
        # Define fls list as BAM files
        bamLst <- BamFileList(contrastFiles, yieldSize=100000, asMates=TRUE)
        
        # create counts object
        contrastSE <- summarizeOverlaps(exonsbygene, bamLst, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
        
        # match the colnames with the vector file names
        colnameIdx <- match(colnames(contrastSE), basename(contrastFiles))
        
        # format runid names based on bamids
        runids <- gsub("#", "_", bamids)

        # assign the runids that correspond to the filenames in colnames
        colnames(contrastSE) <- runids[colnameIdx]

    }
    
    # create a list that holds the current contrast info
    contrastMetadata <- list(contrast=currentContrast)

    # cheekily smuggle the contrast string inside the SE to use later
    metadata(contrastSE) <- contrastMetadata
    
    # right name for RDS file
    if (stepped) {
        
        fileStepName <- paste("step", stepped, "_", sep="")
        
    } else {
        
        fileStepName <- "nosteps_"

    }

    # create input directory (throws a warning if already exists: OK)
    dir.create(inputdir)

    # just in case, save it for later use
    saveRDS(contrastSE, file = paste(wd, "/", inputdir, "/summarizedExperiment_percontrast_", fileStepName, currentContrast, ".Rds", sep=""))
    
    # return the se object
    return(contrastSE)
}

# Sample PCA plot for transformed data
# This plot helps to check for batch effects and the like.
# Code extracted from https://github.com/Bioconductor-mirror/DESeq2/blob/master/R/plots.R
# for customisation
# Author: Wolfgang Huber
myPlotPCA <- function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- genefilter::rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) +
      geom_text(aes(label=group), hjust=0.5, vjust=2, size=2) +
          xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
              ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
                  coord_fixed()
}
