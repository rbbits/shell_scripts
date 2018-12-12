#!/bin/bash

usage(){
    cat <<-EOF
	Wrapper script to run Salmon.
	Output, input, json and staging directories are defined relative to ./ or a working directory provided by -w, and
	if -m runfolder is used, their values will be generated by the script and -[ijos] will be ignored.
	
	Usage: 
	
	$0 [options] targets_file.txt
	
	Options:
	   -0               Dry-run: only print commands to be executed - verbose mode by default.
	   -c <number>      Do not use WTSI composite id (run_position[#tag]), use insted the contents of this column in the targets file (must be unique).
	   -f <format>      Input file format: cram | bam. Default: auto-detect.
	   -h               Show usage message.
	   -i <directory>   Input directory. Default: ./input/.
	   -M <method>      Aligning Method, one of: <hisat2 | star | tophat2>.
	   -m <method hint> Shortcut for specific directory structure: <runfolder | reanalysis>.
	   -n <number>      Line number in targets to be processed
	   -o <directory>   Output directory. Default: ./output/<method>/<run>/<run_pos#tag>/.
	   -r <directory>   Absolute path to repository for reference genome/transcriptome/.
	   -s <directory>   Staging directory. Default: ./staging/<method>/<run>/<run_pos#tag>/.
	   -t <number>      If -m runfolder is used, numeric part of tmp_XXXXX folder
	   -w <directory>   Absolute path to working directory. Default: $PWD.
	   -z               Zip these output files and save them in the output directory: quant.sf, quant.genes.sf, lib_format_counts.json, libParams/, cmd_info.json.
	EOF

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

#set -x

if [ $# -lt 1 ]; then
    usage; 
    exit 1
fi

while getopts ":0c:f:hi:M:m:n:o:r:s:t:vzw:" OPTION; do
    case $OPTION in
        0)
            DRYRUN=1;;
        c)
            TARGETSCOLUMN=$OPTARG
            [[ ! $TARGETSCOLUMN =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -c: not a digit: ${TARGETSCOLUMN}" 1;;
        f)
            FORMAT=$OPTARG;;
        h)
            usage; exit 1;;
        i)
            INPUTDIR=$OPTARG;;
        M)
            METHOD=$OPTARG
            METHODREGEX="^bam2cram|star|tophat2|bwa\_aln|bwa\_mem|hs\_split|y\_split|salmon|hisat2$"
            [[ ! $METHOD =~ $METHODREGEX ]] && exitmessage "[ERROR] -M: invalid method $METHOD: try '$0 -h' for more information" 1;;
        m)
            METHODHINT=$OPTARG
            METHODHINTREGEX="^runfolder|reanalysis$"
            [[ ! $METHODHINT =~ $METHODHINTREGEX ]] && exitmessage "[ERROR] -m: invalid method hint $METHODHINT: try '$0 -h' for more information" 1;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${RECNO}" 1;;
        o)
            OUTPUTDIR=$OPTARG;;
        r)
            REPOSITORY=$OPTARG
            [[ ! -d $REPOSITORY ]] && exitmessage "[ERROR] -w: Cannot access ${REPOSITORY}: no such directory" 2
            [[ ! $REPOSITORY = /* ]] && exitmessage "[ERROR] -w: Not an absolute path" 1;;
        s)
            STAGINGDIR=$OPTARG;;
        t)
            TMPDIRNUM=$OPTARG
            [[ ! $TMPDIRNUM =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${TMPDIRNUM}" 1;;
        v)
            VERBOSE=1;;
        w)
            CWD=$OPTARG
            [[ ! -d $CWD ]] && exitmessage "[ERROR] -w: Cannot access ${CWD}: no such directory" 2
            [[ ! $CWD = /* ]] && exitmessage "[ERROR] -w: Not an absolute path" 1;;
        z)
            ZIPOUTPUT=1;; 
        \?)
            exitmessage "[ERROR] Invalid option: -$OPTARG" 1;;
        :)
            exitmessage "[ERROR] Option -$OPTARG requires an argument." 1;;
    esac
done

shift $((OPTIND-1))
TARGETSFILE=$1

DRYRUN=${DRYRUN-0}
VERBOSE=${VERBOSE-0}
ZIPOUTPUT=${ZIPOUTPUT-0}
INFO_LBL="INFO"

if [ "$DRYRUN" -eq "1" ]; then
    printf "=========\n[DRY-RUN]\n=========\n"
    VERBOSE=1
    INFO_LBL="DRYRUN"
fi
    

BINARY=$(readlink -f `which salmon`)
SALMON=$BINARY

WORKINGDIR="${CWD:-"$PWD"}"
[ "$VERBOSE" -eq 1 ] && printf -- "[INFO] Working directory: ${WORKINGDIR}\n"


# Read info from targets file
if [ -n "$TARGETSFILE" ]; then

    [ ! -f "$TARGETSFILE" ] && exitmessage "[ERROR] Cannot access targets file ${TARGETSFILE}: No such file or directory" 2

    [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1

    LINE=${LSB_JOBINDEX:-$RECNO}

    if [ -n "$TARGETSCOLUMN" ]; then
        BAMID=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk -v column=$TARGETSCOLUMN 'BEGIN { FS = "\t" } ; {print $column}'`
    else
        RUN=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
        POS=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
        TAG=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
        # deal with non-multiplexed lanes
        if [ -z $TAG ]; then
            BAMID="${RUN}_${POS}"
        else
            BAMID="${RUN}_${POS}#${TAG}"
        fi
    fi

    ALIGNMENTSINBAM=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $5}'`
    ALIGNREFGENOME=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $6}'`
    REFDICTNAME=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $7}'`
    REFNAMEFASTA=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $8}'`
    TRANSCRIPTOME=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $9}'`
    TRANSCRIPTANNO=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $10}'`
    LIBRARYTYPE=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $11}'`
    LIBRARYLAYOUT=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $12}'`
    REFTRANSCRIPTFASTA=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $13}'`
else
    exitmessage "[ERROR] A targets file is needed" 1
fi

if [ "$METHODHINT" = "reanalysis" ]; then
    IBAMDIR="input/${RUN}"
    OBAMDIR="output/${METHOD}/${RUN}/${BAMID}"
    SBAMDIR="staging/${METHOD}/${RUN}/${BAMID}"
elif [ "$METHODHINT" = "runfolder" ]; then
    IBAMDIR="no_cal/lane${POSITION}"
    OBAMDIR="no_cal/archive/lane${POSITION}"
    if [ -z "$TMPDIRNUM" ]; then
        exitmessage "[ERROR] -n: a numeric value is required for -n when -m runfolder is being used (numbers in tmp_XXXXXX directory)"
    else
        SBAMDIR="no_cal/archive/tmp_${TMPDIRNUM}/${BAMID}"
    fi
    OUTDATADIR="${WORKINGDIR}/${OBAMDIR}"
    [ -d "$OUTDATADIR" ] || exitmessage "[ERROR] Cannot access ${OUTDATADIR}: No such directory" 2
    OUTSTAGINGDIR="${WORKINGDIR}/${SBAMDIR}"
    [ -d $OUTSTAGINGDIR ] || exitmessage "[ERROR] Cannot access ${OUTSTAGINGDIR}: No such directory" 2
else
    IBAMDIR="${INPUTDIR:-"input/$RUN"}"
    OBAMDIR="${OUTPUTDIR:-"output/$METHOD/$RUN/$BAMID"}"
    SBAMDIR="${STAGINGDIR:-"staging/$METHOD/$RUN/$BAMID"}"
fi

OUTDATADIR="${WORKINGDIR}/${OBAMDIR}"
[ ! -d "$OUTDATADIR" ] && printf -- "[WARNING] Cannot access ${OUTDATADIR}: No such directory\n" 2
[ "$VERBOSE" -eq 1 ] && printf "[INFO] Output directory: %s\n" "${OUTDATADIR}"

INDATADIR="${WORKINGDIR}/${IBAMDIR}"
[ ! -d "$INDATADIR" ] && exitmessage "[ERROR] Cannot access ${INDATADIR}: No such directory" 2
[ "$VERBOSE" -eq 1 ] && printf "[INFO] Input directory: %s\n" "${INDATADIR}"

OUTSTAGINGDIR="${WORKINGDIR}/${SBAMDIR}"
[ ! -d "$OUTSTAGINGDIR" ] && printf -- "[WARNING] Cannot access ${OUTSTAGINGDIR}: No such directory\n" 2
[ "$VERBOSE" -eq 1 ] &&  printf "[INFO] Staging directory: %s\n" "${OUTSTAGINGDIR}"

SRCINPUT="${INDATADIR}/${BAMID}"
REPOSDIR="${REPOSITORY-"/lustre/scratch117/core/sciops_repository"}"
ALIGNMENTMETHOD=$METHOD

SALMON_ARGS="--no-version-check quant --libType A -p 8 "
    
[[ ! $TRANSCRIPTANNO = /* ]] && SALMON_ARGS+="--geneMap ${REPOSDIR}/transcriptomes/${TRANSCRIPTANNO} " || SALMON_ARGS+="--geneMap ${TRANSCRIPTANNO} "

if [[ ! $TRANSCRIPTOME = */salmon ]]; then
    SALMON_TRANSCRIPTOME="$(dirname $TRANSCRIPTOME)"
    SALMON_TRANSCRIPTOME="$(dirname $SALMON_TRANSCRIPTOME)"
    SALMON_TRANSCRIPTOME+="/salmon"
else
    SALMON_TRANSCRIPTOME="$TRANSCRIPTOME"
fi

[[ ! $TRANSCRIPTOME = /* ]] && SALMON_ARGS+="--index ${REPOSDIR}/transcriptomes/${SALMON_TRANSCRIPTOME} " || SALMON_ARGS+="--index ${SALMON_TRANSCRIPTOME} "

SALMON_ARGS+="--mates1 ${OUTSTAGINGDIR}/intfile_1_${BAMID}.fq "
SALMON_ARGS+="--mates2 ${OUTSTAGINGDIR}/intfile_2_${BAMID}.fq "

SALMON_OUTPUT_DIR="${OUTSTAGINGDIR}/salmon_quant_${BAMID}_${LSB_JOBID}"

SALMON_ARGS+="--output ${SALMON_OUTPUT_DIR} "

SALMON_CMD="$SALMON $SALMON_ARGS"

[ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${SALMON_CMD}"

if [ "$DRYRUN" -eq 0 ]; then
    CMD="$(${SALMON_CMD})"
    RET_CODE=$?
else
    RET_CODE=0
fi

if [ "$RET_CODE" -eq "0" ]; then
    
    printf "[${INFO_LBL}] Done Salmon for [%s] = [OK]\n" "${BAMID}"

    if [ "$ZIPOUTPUT" -eq 1 ]; then

        ZIP=`which zip`
        ZIP_CMD="${ZIP} -r ${OUTDATADIR}/${BAMID}_salmon.quant.zip "
        ZIP_CMD+="${SALMON_OUTPUT_DIR}/quant.sf "
        ZIP_CMD+="${SALMON_OUTPUT_DIR}/quant.genes.sf "
        ZIP_CMD+="${SALMON_OUTPUT_DIR}/lib_format_counts.json "
        ZIP_CMD+="${SALMON_OUTPUT_DIR}/libParams "
        ZIP_CMD+="${SALMON_OUTPUT_DIR}/cmd_info.json"

        [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${ZIP_CMD}"

        if [ "$DRYRUN" -eq 0 ]; then
            CMD="$(${ZIP_CMD})"
            RET_CODE=$?
        else
            RET_CODE=0
        fi

        if [ "$RET_CODE" -eq "0" ]; then
            exitmessage "[${INFO_LBL}] Done Zip for [${BAMID}] = [OK]" 0
        else
            exitmessage "[ERROR] Failed Zip for [${BAMID}] exited with exit code ${RET_CODE}: ${CMD}" $RET_CODE
        fi

    else

        exit 0

    fi
    
else
    
    exitmessage "[ERROR] Failed Salmon for [${BAMID}] exited with exit code ${RET_CODE}: ${CMD}" $RET_CODE

fi
