#!/bin/bash

#set -x

usage(){
    cat <<-EOF
	NAME

	   run_aspera.sh 

	SYNOPSIS

	   Wrapper script for ascp. Accepts a specific accession number or
	   a targets file with multiple run accession numbers.

	   $0 [options] <accession>

	OPTIONS

	   -h               Help
	   -o <directory>   Output directory. Default is ./

	   If you want to use a targets file then the following options are required:

	   -c <number>      Column number in targets file that contains SSR accession numbers
	   -n <number>      Line number in targets to be processed
	   -t <file>        File containing targets or an SSR accession number

	EOF
}

exitmessage(){
    declare EXITMESG=$1
    declare EXITCODE=$2
    if [ $EXITCODE -eq 0 ]; then
        >&1 printf '\n%s\n' "$EXITMESG"
    else
        >&2 printf '\n%s\n' "$EXITMESG"
    fi
    exit $EXITCODE
}

while getopts ":0c:hn:o:t:" OPTION; do
    case $OPTION in
        0)
            DRYRUN=1;;
        c)
            TARGETSCOLUMN=$OPTARG
            [[ ! $TARGETSCOLUMN =~ ^[[:digit:]]+$ ]] && exitmessage "-c: not a digit: ${TARGETSCOLUMN}" 1;;
        h)
            usage; exit 1;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${RECNO}" 1;;
        o)
            OUTPUT_DIR=$OPTARG
            [ -d $OUTPUT_DIR ] || exitmessage "[ERROR] Cannot access ${OUTPUT_DIR}: No such file or directory" 2;;
        t)
            TARGETSFILE=$OPTARG;;
        \?)
            exitmessage "[ERROR] Invalid option: -${OPTARG}" 1;;
        :)
            exitmessage "[ERROR] Option -${OPTARG} requires an argument." 1;;
    esac    
done

shift $((OPTIND-1))
ACCESSION=$1
ACCESSION_CPY=$ACCESSION
BINARY="/software/aspera/bin/ascp"
DSA_KEY="/software/aspera/etc/asperaweb_id_dsa.openssh"
OUTPUT_DIR=${OUTPUT_DIR:-"./"}

[ ! -e "$BINARY" ] || [ ! -e "$DSA_KEY" ] && exitmessage "[ERROR] Cannot access ${BINARY} or ${DSA_KEY}: No such file or directory" 2

if [ -n "$TARGETSFILE" ]; then

    [ ! -f "$TARGETSFILE" ] && exitmessage "[ERROR] Cannot access ${TARGETSFILE}: No such file or directory" 2
    [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1
    [ -z "$TARGETSCOLUMN" ] && exitmessage "[ERROR] -c: a column number is required when using a targets file" 1
    LINE=${LSB_JOBINDEX:-$RECNO}
    ACCESSION=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk -v COLUMN=${TARGETSCOLUMN} 'BEGIN { FS = "\t" } ; {print $COLUMN}'`
    OUTPUT_DIR+="/${ACCESSION}/"
    mkdir -p "${OUTPUT_DIR}" && printf -- "\n[INFO] Output directory ${OUTPUT_DIR} was created\n"
    [ -n "$ACCESSION_CPY" ] && [ "$ACCESSION" != "$ACCESSION_CPY" ] && printf -- "\n[WARNING] $0: Requested accession number ${ACCESSION_CPY} was changed to ${ACCESSION}\n"
    
fi

[ -z "$ACCESSION" ] && exitmessage "[ERROR] $0: an accession number is required" 1

SSR="${ACCESSION:0:6}"
SRA_FILE="/sra/sra-instant/reads/ByRun/sra/SRR/${SSR}/${ACCESSION}/${ACCESSION}.sra"

#------------#
# command
#------------#
ASCP_CMD="$BINARY "
ASCP_CMD+="-i ${DSA_KEY} "
ASCP_CMD+="-k 1 -T -l200m "
ASCP_CMD+="anonftp@ftp.ncbi.nlm.nih.gov:${SRA_FILE} "
ASCP_CMD+="$OUTPUT_DIR"

ASCP="$($ASCP_CMD 2>&1)"

RET_CODE=$?

[ "$RET_CODE" -eq 0 ] && printf -- "\n[INFO] ${ASCP}\n\n[INFO] Done\n" || exitmessage "[INFO] Exited with exit code ${RET_CODE}: ${ASCP}" $RET_CODE

