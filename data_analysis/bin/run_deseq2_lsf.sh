#!/bin/bash

usage(){
    cat <<-EOF
	Usage: $0 options
	
	This script runs a DESeq2 analysis using LSF
	
	OPTIONS
	
	Required: 
	
	   -s <step>        A number specifying the step to run: 0,1,2. 
	                    Default: 0 (no steps).
	   -y <yaml file>   Path to configuration file in YAML format.
	
	Optional:
	
	   -J <indexes>   String indicating the list of indexes to be processed 
	                  in LSF command. E.g. -J 1-12 -J 10-11 -J 1,2,4,5 -J 4
	                  Default: process all records in targets file.
	   -r <rds file>  Path to RDS file.
	   -w <job id>    Only for step 2: Use this job ID as dependency for 
	                  DESeq2
	   -0             Dry-run.
	
	email to <rb11@sanger.ac.uk> to report bugs.
	EOF
}

while getopts ":0hJ:r:s:w:y:" OPTION; do
    case $OPTION in
        0)
            DRYRUN=1;;
        h)
            usage; exit 1;;
        J)
            regexN='^[[:digit:]]*$'
            regexNdN='^([[:digit:]]*)-([[:digit:]]*)$'
            regexNcN='^[[:digit:]]*(,?[[:digit:]]*)*$'
            if [[ "$OPTARG" =~ $regexN ]]; then
                ARRAY_IDXS=$OPTARG
            elif [[ "$OPTARG" =~ $regexNdN ]]; then
                d1="${BASH_REMATCH[1]}"
                d2="${BASH_REMATCH[2]}"
                [ "$d1" -gt "$d2" ] && echo "Bad -J: ${d1} bigger than ${d2}." >&2; exit 1
                ARRAY_IDXS=$OPTARG
            elif [[ "$OPTARG" =~ $regexNcN ]]; then
                ARRAY_IDXS=$OPTARG
            else
                echo "Bad -J: \"${OPTARG}\"." >&2; exit 1
            fi;;
        r)
            RDS_FILE=$OPTARG;;
        s)
            STEP=$OPTARG;;
        w)
            regexN='^[[:digit:]]*$'
            [[ "$OPTARG" =~ $regexN ]] && JOBID_PARAM=$OPTARG || echo "Bad -w: ${OPTARG} must be a number." >&2; exit 1;;
        y)
            # useful config file in YAML format contains all the info
            # needed to process the files in a single step.
            YAML_FILE=$OPTARG;;
        \?)
            echo "Invalid option: -$OPTARG" >&2; exit 1;;
        :)
            echo "Option -$OPTARG requires an argument." >&2; exit 1;;
    esac
done


# Function to parse a simple YAML file
# Author: Stefan Farestam. WWW: http://is.gd/aQeT2O
parse_yaml (){
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

# function to print messages
# to STDOUT or STDERR
exitmessage(){
    declare EXITMESG=$1
    declare EXITCODE=$2
    if [ $EXITCODE -eq 0 ]; then
        >&1 printf '%s\n' "$EXITMESG"
    else
        >&2 printf '%s\n' "$EXITMESG"
    fi
    exit $EXITCODE
}


[[ -z $YAML_FILE || -z $STEP ]] && exitmessage "$0: A YAML file and a Step are required use: -y /path/to/yaml/file -s 0|1|2" 1


# Variables set inside the YAML file (unless otherwise set inside this script):
# targets, contrasts, bamsource, merged, inputDir, outputDir, pooling, tophatDir, irodsDir, varOfInterest, species
eval $(parse_yaml $YAML_FILE "YML_")

# A stepped way aims to contain the memory-hungry process of producing the
# SummarizedExperiment (SE) object by getting R to do it on a file by file basis.
# Stage 1 creates SE objects and serialise them to be saved as Rds files.
# Stage 2 reads these files and merges them into a single SE object to be used by DESeq2

DRYRUN=${DRYRUN-0}
RET_CODE=0
WD=$PWD
ANALYSIS_DATE=`date +%Y-%m-%d`
ANALYSIS_USER=`whoami`
EXPERIMENT_ID=$YML_experimentId

printf "[INFO] Working directory: %s\n" "${WD}"


[[ "$DRYRUN" -eq "0" ]] && VERBOSELBL="INFO" || VERBOSELBL="DRYRUN"

# If we are running the script several times, we save time by saving the main
# SerializedExperiment object (in the R script) in which case we provide the
# /path/to/rdsfile in the command line args via the sedata variable    

if [[ "$STEP" -eq "0" &&  -n "$RDS_FILE" ]] || [[ "$STEP" -eq "2" &&  -n "$RDS_FILE" ]]; then
        
    # This is only useful if:
    # a) the script was run with -s 0 (no steps) and the R script died
    # b) step 1 was run successfully, then 2 was run the R script died
    if [ ! -e "$RDS_FILE" ]; then
        exitmessage "You specified an RDS file must be used, but no such file exist: ${RDS_FILE}" 1
    else
        printf "[${VERBOSELBL}] Using RDS file: %s\n" "${RDS_FILE}"
        RDS_FILE_ARG="--sedata=${RDS_FILE}"
    fi

elif [ "$STEP" -eq "0" ] || [ "$STEP" -eq "1" ]; then
            
    # if the samples have been merged they have been pre-processed and made ready for the R script
    # otherwise we need to extract the accepted hits and make them available for the R script
    if [ "$YML_bamsource" = "irodstophat" ]; then
                
        # If the source of the bams is irods and they were produced using Tophat,
        # as most RNA-Seq samples are, then we need to extract the accepted hits out
        # of them. The index is needed for viewing using IGV
        OUT_SH_FILE="./resources/extract_accepted_hits.sh"
        
        (cat <<-'EOF'
		#!/bin/bash
		. /software/npg/etc/profile.npg
		IDIR=$1; ODIR=$2; TARGETS=$3;
		let "LINE=$LSB_JOBINDEX+1" ##header offset
		RUN=`head -n ${LINE} ${TARGETS} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
		POS=`head -n ${LINE} ${TARGETS} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
		TAG=`head -n ${LINE} ${TARGETS} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
		[ -n $TAG ] && BAMID="${RUN}_${POS}#${TAG}" || BAMID="${RUN}_${POS}"
		INPUT_BAM="${IDIR}/${RUN}/${BAMID}.cram"
		OUTPUT_BAM="${ODIR}/${BAMID}_accepted_hits.bam"
		RET_CODE=0
		mkdir -pv $ODIR
		if [ ! -e "${OUTPUT_BAM}.bai" ]; then
		   CMD="samtools1 view -bh -F 0x4 -o $OUTPUT_BAM $INPUT_BAM"
		   printf "[INFO] Running view command:\n[INFO] %s\n" "${CMD}"
		   SAMTOOLSCMD=$($CMD)
		   RET_CODE=$?
		   if [ "$RET_CODE" -eq "0" ]; then
		       CMD="samtools1 index $OUTPUT_BAM"
		       printf "[INFO] Running index command:\n[INFO] %s\n" "${CMD}"
		       SAMTOOLSCMD=$($CMD)
		       RET_CODE=$?
		       [ "$RET_CODE" -ne "0" ] && printf -- "[ERROR] Samtools Index for RunId [${BAMID}] failed!\n" >&2
		   else       
		       printf -- "[ERROR] Samtools View for RunId [${BAMID}] failed!\n" >&2
		   fi
		else
		   printf -- "[INFO] Accepted hits for [$BAMID] have been extracted already\n"
		fi
		[ "$RET_CODE" -eq "0" ] && printf -- "[INFO] Done\n" || printf -- "[ERROR] Exited\n" >&2
		exit $RET_CODE
		EOF
        ) > $OUT_SH_FILE

        if [ -f "$OUT_SH_FILE" ]; then
                
            chmod 755 $OUT_SH_FILE
            printf "[${VERBOSELBL}] Auxiliary script file %s was created successfully\n" "${OUT_SH_FILE}"
            NUM_TARGETS=`grep -c ^[0-9] $YML_targets`
            JOB_ARRAY_LENGTH=${ARRAY_IDXS:-"1-${NUM_TARGETS}"}
            BSUBCMD="bsub -J ${EXPERIMENT_ID}_${ANALYSIS_USER}_${ANALYSIS_DATE}_XAHITS[${JOB_ARRAY_LENGTH}] "
            BSUBCMD+="-oo log/extract_accepted_hits.%I.o -eo log/extract_accepted_hits.%I.e "
            BSUBCMD+="-R select[mem>8000] -R rusage[mem=8000] -M 8000 $OUT_SH_FILE ${YML_irodsDir} ${YML_inputDir} ${YML_targets}"
            
            printf "[${VERBOSELBL}] BSUB COMMAND: %s\n" "${BSUBCMD}"
            
            if [ "$DRYRUN" -eq 1 ]; then
                BAM_JOBID="dummy"
            else
                LSFCMD=$($BSUBCMD) && LSFCMD_RET_CODE=$?
                BAM_JOBID=`grep -oP "(?<=\<)[[:digit:]]*(?=\>)" <<< $LSFCMD`
                [ "$LSFCMD_RET_CODE" -ne "0" ] && exitmessage "Bsub command failed to execute: \"$LSFCMD\"" $LSFCMD_RET_CODE || printf "[INFO] LSF command output: %s\n" "${LSFCMD}"
            fi
            
            # Give LSF time to submit the job
            printf "[${VERBOSELBL}] Waiting for job <%s> to be submitted... " "${BAM_JOBID}"
            [ "$DRYRUN" -eq "0" ] && sleep 2
            printf -- "[${VERBOSELBL}] Done\n"
            printf "[${VERBOSELBL}] Job id <%s> will be used as dependency for DESeq2 job\n" "${BAM_JOBID}"
            
        else
            
            exitmessage "Problem in creating bash file: \"$OUTSHFILE\"" 1
            
        fi
            
    elif [ "$YML_bamsource" = "tophat" ]; then
            
        # If the source of the bams is an accepted hits bam file produced by Tophat
        # run by us (e.g. from fastq files obtained from files downloaded from irods)
        # then we create soft links to these inside the DESeq2 input directory
        
        [ "$DRYRUN" -eq "0" ] && mkdir -p "${YML_inputDir}"
        printf "[${VERBOSELBL}] Created dir: %s\n" "${YML_inputDir}"
        [ "$DRYRUN" -eq "0" ] && . /nfs/users/nfs_r/rb11/local/bin/chdir "${YML_inputDir}"
        printf "[${VERBOSELBL}] Changed dir to: %s\n" "${YML_inputDir}"
            
        # Read ech line in the targets file provided by the yaml file
        while read line; do
            # read run info from targets
            RUN=`echo "$line" | awk -F'\t' '{print $2}'`
            POS=`echo "$line" | awk -F'\t' '{print $3}'`
            TAG=`echo "$line" | awk -F'\t' '{print $4}'`
            SMP=`echo "$line" | awk -F'\t' '{print $6}'`
            # deal with non-multiplexed lanes: the combination of run_pos#tag
            # is used as bam identifier for both merged and not merged bams            
            [ -n "$TAG" ] && BAMID="${RUN}_${POS}#${TAG}" || BAMID="${RUN}_${POS}"
            
            if [ "$YML_merged" -eq "1" ]; then
                # as of today (this may change in the future) it's assumed that the merged
                # files were originated from an irods->fastq->tophat->merge process and
                # therefore we already have accepted hits that need to be soft linked to;
                # we use the sample name as input bam name, this value should be at the 6th
                # column of the targets file                
                BAM_DIR="${YML_tophatDir}/merged"
                TARGET_BAM="${SMP}_accepted_hits.bam"
                INPUT_BAM=$TARGET_BAM
                MERGED_ARG="--merged"
            else
                BAM_DIR="${YML_tophatDir}/${RUN}/${BAMID}"
                TARGET_BAM="accepted_hits.bam"
                INPUT_BAM="${BAMID}_accepted_hits.bam"
            fi
            [ "$DRYRUN" -eq "0" ] && ln -s -f "../../${BAMDIR}/${TARGET_BAM} ${INPUT_BAM}"
            printf "[${VERBOSELBL}] Created soft link: %s\n" "${INPUT_BAM}"
        done <$YML_targets
            
        [ "$DRYRUN" -eq "0" ] && . /nfs/users/nfs_r/rb11/local/bin/chdir "${WD}"
        printf "[DRYRUN] Changed dir back to %s\n" "${WD}"
        BAM_JOBID=""
        
    else
        
        echo "Process STAR or others here ... bye."
        exit 1
        
    fi # (bam source)

    # pass different arguments to the R script depending on the step
    [ "$STEP" -eq "0" ] && RSCRIPT_ARGS="--idir=${YML_inputDir} --odir=${YML_outputDir} --contrasts=${YML_contrasts} --targets=${YML_targets} --pooling=${YML_pooling} --varofinterest=${YML_varOfInterest} "
    [ "$STEP" -eq "1" ] && RSCRIPT_ARGS="--ibam=${YML_inputDir}/${INPUT_BAM} --bamid=${BAMID} --irds=${YML_inputDir} "

    # pass common arguments for steps 0 and 1
    RSCRIPT_ARGS+="--step=${STEP} --species=${YML_species} ${MERGED_ARG} --experimentid=${EXPERIMENT_ID} "

    # pass common arguments only if they exist
    [[ -n $YML_genome ]] && RSCRIPT_ARGS+="--genome=${YML_genome} "
    [[ -n $YML_transcriptome ]] && RSCRIPT_ARGS+="--transcriptome=${YML_transcriptome} "

elif [ "$STEP" -eq "2" ] ; then

    [ "$YML_merged" -eq "1" ] && MERGED_ARG="--merged"
    
    RSCRIPT_ARGS="--irds=${YML_inputDir} --odir=${YML_outputDir} --contrasts=${YML_contrasts} "
    RSCRIPT_ARGS+="--targets=${YML_targets} --pooling=${YML_pooling} --step=${STEP} "
    RSCRIPT_ARGS+="--varofinterest=${YML_varOfInterest} --species=${YML_species} --experimentid=${EXPERIMENT_ID} "
    RSCRIPT_ARGS+="${MERGED_ARG}"

    if [[ -n $JOBID_PARAM ]]; then

        printf "[${VERBOSELBL}] Using job ID <%s> as dependency for the DESeq2 job\n" "${JOBID_PARAM}"
        BAM_JOBID=$JOBID_PARAM

    else
        
        printf -- "[INFO] Looking for Accepted Hits extraction jobs ...\n"
        XAHITS_IDS=(`bjobs -J "${EXPERIMENT_ID}_${ANALYSIS_USER}_${ANALYSIS_DATE}_XAHITS" | awk '{print $1}' | grep -v "JOBID"`)

        if [ "${#XAHITS_IDS[@]}" -gt "0" ]; then
            REPLY=""
            while [[ ! $REPLY =~ 1|2 ]]; do
                printf -- "At least one LSF job for extracting accepted hits was found. Do you want to:\n"
                printf -- "[1] Use it as dependency for DESeq2 job\n"
                printf -- "[2] Investigate what's going on, abort DESeq2 for now\n"
                printf -- "(NOTE: avoid this prompt by using option -w JOBID when running this script)\n"
                printf -- "Your choice: "
                read -n 1 REPLY
                if [ "$REPLY" -eq "1" ]; then
                    BAM_JOBID="${XAHITS_IDS[0]}"
                    printf "[INFO] Job id <%s> will be used as dependency for DESeq2 job\n" "${BAM_JOBID}"
                else
                    exitmessage "[INFO] Aborted execution of DESeq2 R script. Bye!" 1
                fi
            done
        else
            printf -- "[INFO] No other jobs found\n"
            BAM_JOBID=""
        fi

    fi
    
else
    
    exitmessage "[ERROR] Unknown step \"${STEP}\"" 1
    
fi

RSCRIPT_ARGS+=" ${SE_DATA_ARG}"

BSUBCMD="bsub -J ${EXPERIMENT_ID}_${ANALYSIS_USER}_${ANALYSIS_DATE}_DESEQ2 "
BSUBCMD+="-oo log/deseq2_${ANALYSIS_DATE}.o -eo log/deseq2_${ANALYSIS_DATE}.e "
BSUBCMD+="-R select[mem>22000] -R rusage[mem=22000] -M 22000 "
BSUBCMD+="-n 4 -R span[hosts=1] "

[[ -n $BAM_JOBID ]] && BSUBCMD+="-w done(${BAM_JOBID}) "
    
printf "[${VERBOSELBL}] Ready to run DESeq2 R script step: %s with arguments:\n" "${STEP}"

BSUBCMD+="/software/R-3.3.0/bin/Rscript ./bin/run_deseq2.R ${RSCRIPT_ARGS}"

printf "[${VERBOSELBL}] Bsub commmand : %s\n" "${BSUBCMD}"

if [ "$DRYRUN" -eq 1 ]; then

    LSFCMD_RET_CODE=0
    
else

    LSFCMD=$($BSUBCMD) && let "LSFCMD_RET_CODE+=$?"
    [ "$LSFCMD_RET_CODE" -ne "0" ] && exitmessage "Bsub command failed to execute: \"$LSFCMD\"" $LSFCMD_RET_CODE || printf "[INFO] LSF command output: %s\n" "${LSFCMD}"

fi

printf -- "[${VERBOSELBL}] Done\n"

exit $LSFCMD_RET_CODE
