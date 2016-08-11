#!/bin/bash

if [ -e /software/npg/etc/profile.npg ]; then
    . /software/npg/etc/profile.npg
fi

if [ -e /software/sanger-samtools-refpath/etc/profile.sanger-samtools-refpath ]; then
    . /software/sanger-samtools-refpath/etc/profile.sanger-samtools-refpath
fi

# Function to parse a simple YAML file
# Author: Stefan Farestam. WWW: http://is.gd/aQeT2O
function parse_yaml {
    local prefix=$2
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
        awk -F$fs '{
           indent = length($1)/2;
           vname[indent] = $2;
           for (i in vname) {if (i > indent) {delete vname[i]}}
           if (length($3) > 0) {
             vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
             printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
           }
        }'
}


# useful config file in YAML format contains all the info
# needed to process the files in a single step.
confYaml=$1

eval $(parse_yaml $confYaml)

# Variables set inside the YAML file (unless otherwise set inside this script):
# targets, contrasts, bamsource, merged, inputDir, outputDir, pooling, tophatDir, irodsDir, varOfInterest, species

# A stepped way aims to contain the memory-hungry process of producing the
# SummarizedExperiment (SE) object by getting R to do it on a file by file basis.
# Stage 1 creates SE objects and serialise them to be saved as Rds files.
# Stage 2 reads these files and merges them into a single SE object to be used by DESeq2

ret_code=0
wd=$PWD

if [ $2 = "nosteps" ]; then # ... or process all bam at once

    
    
    # If we are running the script several times, we save time by saving the main
    # SerializedExperiment object (in the R script) in which case we provide the
    # /path/to/rdsfile in the command line args via the sedata variable
    if [ -n "$3" ]; then
        
        if [ "$pooling" = "allsamples" ]; then

            if [ ! -e "${inputDir}/summarizedExperiment_nosteps_allsamples.Rds" ]; then
                echo "You specified an RDS file must be used, but no such file exist: ${inputDir}/summarizedExperiment_nosteps_allsamples.Rds"
                exit 1
            fi
            
        elif [ "$pooling" = "percontrast" ]; then
            
            if [ ! -e "$3" ]; then
                echo "You specified an RDS file must be used, but no such file exist: ${3}"
                exit 1
            fi
            
        fi

        sedata_arg=" --sedata=$3"
        
    else

        # if the samples have been merged they have been pre-processed and made ready for the R script
        # otherwise we need to extract the accepted hits and make them available for the R script
        if [ "$merged" -eq 0 ]; then
            
            # read ech line in the targets file provided by the yaml file
            while read line; do
                
                if [[ $line =~ ^[0-9]+.*$ ]]; then
                    
                    # read run info from targets
                    run=`echo "$line" | awk -F'\t' '{print $2}'`
                    position=`echo "$line" | awk -F'\t' '{print $3}'`
                    tag=`echo "$line" | awk -F'\t' '{print $4}'`
                    
                    # deal with non-multiplexed lanes
                    if [ -z "$tag" ]; then
                        bamid=$run\_$position
                    else
                        bamid=$run\_$position\#$tag
                    fi
                    
                    # create DESeq2 input dir
                    mkdir -p $inputDir
                    
                    # If the source of the bams is irods and they were produced using Tophat,
                    # as most RNA-Seq samples are, then we need to extract the accepted hits out
                    # of them. The index is needed for viewing using IGV.
                    # If the source of the bams is an accepted hits bam file produced by Tophat
                    # run by us (e.g. from fastq files obtained from files downloaded from irods)
                    # then we create soft links to these inside the DESeq2 input directory
                    if [ $bamsource = "irodstophat" ]; then
                        
                        inputBam=$irodsDir/$run/${bamid}.cram
                        outputBam=$inputDir/${bamid}_accepted_hits.bam
                        if [ ! -e "${outputBam}.bai" ]; then
                            samtools1 view -bh -F 0x4 -o $outputBam $inputBam && samtools1 index $outputBam
                            ret_code=$?
                            if [ $ret_code -ne 0 ]; then
                                echo "Samtools failed!"
                                exit $ret_code
                            fi
                        fi
                        
                    elif [ $bamsource = "tophat" ]; then
                        
                        bamfile=${bamid}_accepted_hits.bam
                        bamdir=$inputDir
                        # go to deseq/input dir
                        . /nfs/users/nfs_r/rb11/local/bin/chdir $inputDir
                        ln -s -f ../../${tophatDir}/${run}/${bamid}/accepted_hits.bam $bamfile
                        . /nfs/users/nfs_r/rb11/local/bin/chdir $wd
                        
                    else
                        
                        echo "Process STAR or others here ... bye."
                        exit 1
                        
                    fi
                    
                fi
                
            done <$targets
        
        fi # merge

    fi # $3 (sedata)
     
    arg="--idir=$inputDir --odir=$outputDir --contrasts=$contrasts --targets=$targets --pooling=$pooling --step=0 --varofinterest=$varOfInterest --species=$species"

    [ $merged -eq 1 ] && arg+=" --merged"

    [ -n "$sedata_arg" ] && arg+="$sedata_arg"
        
    stage="no steps"


    
elif [ $2 = "step1" ] ; then

    

    let "line=$LSB_JOBINDEX+1"
    
    if [ "$merged" -eq 1 ]; then

        sample=`head -n ${line} ${targets} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $6}'`
        bamid=$sample
        bamdir=$tophatDir/merged
        bamfile=${sample}_accepted_hits.bam
        arg="--merged "
        
    else
        
        run=`head -n ${line} ${targets} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'` 
        lane=`head -n ${line} ${targets} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
        tag=`head -n ${line} ${targets} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'` 

        # deal with non-multiplexed lanes
        if [ -z $tag ]; then
            bamid=$run\_$lane
        else
            bamid=$run\_$lane\#$tag
        fi

        bamdir=$inputDir
        bamfile=${bamid}_accepted_hits.bam

        mkdir -p $inputDir
        
        if [ $bamsource = "irodstophat" ]; then

            # QUESTION: what happens if I make irodsDir=./irods/$run in the config file?
            inputBam=$irodsDir/$run/${bamid}.cram
            outputBam=$inputDir/${bamid}_accepted_hits.bam         
            if [ ! -e "${outputBam}.bai" ]; then
                samtools1 view -bh -F 0x4 -o $outputBam $inputBam && samtools1 index $outputBam
                ret_code=$?
                if [ $ret_code -ne 0 ]; then
                    echo "Samtools failed!"
                    exit $ret_code
                fi
            fi
           
        elif [ $bamsource = "tophat" ]; then
            
            # go to deseq/input dir
            . /nfs/users/nfs_r/rb11/local/bin/chdir $inputDir
            ln -s -f ../../${tophatDir}/${run}/${bamid}/accepted_hits.bam $bamfile
            . /nfs/users/nfs_r/rb11/local/bin/chdir $wd

        else
            
            echo "Process STAR or others here ... bye."
            exit 1
            
        fi

    fi
    
    inputBam=$bamdir/$bamfile
    arg+="--ibam=$inputBam --bamid=$bamid --irds=$inputDir --step=1 --species=$species"
    stage="_1step"

    
        
elif [ $2 = "step2" ] ; then

    arg="--irds=$inputDir --odir=$outputDir --contrasts=$contrasts --targets=$targets --pooling=$pooling --step=2 --varofinterest=$varOfInterest --species=$species"
    
    # in we are running the script several times, we save time by saving the main
    # SerializedExperiment object (in the R script) so as to save processing time
    # in which case we provide the /path/to/rdsfile in the command line args
    if [ -n "$3" ]; then
        arg+=" --sedata=$3"        
    fi

    if [ "$merged" -eq 1 ]; then
        arg+=" --merged"        
    fi
    
    stage="step 2"
    

else
    

    echo Unknown argument \"${2}\"
    exit  1

    
fi


echo Ready to run DESeq2 R script Stage: \"${stage}\"

echo "Command: /software/R-3.2.2/bin/Rscript ./bin/run_deseq2.R $arg"

echo Working directory: $PWD

/software/R-3.2.2/bin/Rscript ./bin/run_deseq2.R $arg

ret_code=$?

echo Done.

exit $ret_code
