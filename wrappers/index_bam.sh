#!/bin/bash

usage(){
    cat <<-EOF
	Usage: 

	$0 [options] targets_file.txt
	
	Options:
	   -c <number>      Do not use WTSI composite id (run_position[#tag]), use insted the
	                    contents of this column in the targets file (must be unique).
	   -f <format>      Input file format: cram | bam. Default: auto-detect in this order: cram, bam, sam.
	   -h               Show usage message.
	   -i <directory>   Input directory where the [b|cr]am file is located. Default: <current|working directory>.
	   -n <number>      Line number in targets to be processed
	   -w <directory>   Absolute path to working directory. Default: $PWD
	   -x <extra args>  Extra arguments passed to samtools in a quoted string.
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

while getopts ":c:f:hi:n:w:x:" OPTION; do
    case $OPTION in
        c)
            TARGETSCOLUMN=$OPTARG
            [[ ! $TARGETSCOLUMN =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -c: not a digit: ${TARGETSCOLUMN}" 1;;
        f)
            FORMAT=$OPTARG;;
        h)
            usage; exit 1;;
        i)
            INPUTDIR=$OPTARG
            [ -d $INPUTDIR ] || exitmessage "[ERROR] Cannot access ${INDIR}: No such directory" 2;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${RECNO}" 1;;
        w)
            CWD=$OPTARG
            [[ ! -d $CWD ]] && exitmessage "[ERROR] -w: Cannot access ${CWD}: no such directory" 2
            [[ ! $CWD = /* ]] && exitmessage "[ERROR] -w: Not an absolute path" 1;;
        x)
            EXTRAKEYVALS=$OPTARG;;
        \?)
            exitmessage "[ERROR] Invalid option: -$OPTARG" 1;;
        :)
            exitmessage "[ERROR] Option -$OPTARG requires an argument." 1;;
    esac
done

shift $((OPTIND-1))
TARGETSFILE=$1

BINARY=$(readlink -f `which samtools`)
SAMTOOLS=$BINARY

WORKINGDIR=${CWD:-"$PWD"}
printf -- "[INFO] Working directory: ${WORKINGDIR}\n"

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
else
    exitmessage "[ERROR] A targets file is needed" 1
fi

IBAMDIR=${INPUTDIR:-"."}
OUTPUTDIR=$IBAMDIR
IBAM="${IBAMDIR}/${BAMID}"

printf "[INFO] Output directory: %s\n" "${IBAMDIR}"
printf "[INFO] Input directory: %s\n" "${OUTPUTDIR}"

if [ -z $FORMAT ]; then
    if [ -e "${IBAM}.cram" ]; then
        IBAM+=".cram"
        FORMAT="cram"
        OBAM="${IBAMDIR}/${BAMID}.crai"
    elif [ -e "${IBAM}.bam" ]; then
        IBAM+=".bam"
        FORMAT="bam"
        OBAM="${IBAM}.bai"
    elif [ -e "${IBAM}.sam" ]; then
        IBAM+=".sam"
        FORMAT="sam"
        OBAM="${IBAM}.sai"
    else
        exitmessage "[ERROR] ${BAMID}: No cram or bam or sam file was found in ${IBAMDIR}" 1
    fi
else
    IBAM+=".${FORMAT}"
fi

printf "[INFO] Input format: %s\n" "${FORMAT}"

[ ! -z "$EXTRAOPTS" ] && printf -- "[INFO] Using extra arguments [ ${EXTRAOPTS} ]\n"

CMD="${SAMTOOLS} index ${IBAM} ${OBAM}"

printf "[COMMAND] %s\n" "${CMD}"

INDEX="$($CMD 2>&1)"
        
RET_CODE=$?

if [ "$RET_CODE" -eq "0" ]; then
    exitmessage "[INFO] Done index for [${BAMID}] = [OK]" 0
else
    exitmessage "[ERROR] Failed index for [${BAMID}] exited with exit code ${RET_CODE}: ${INDEX}" $RET_CODE
fi
