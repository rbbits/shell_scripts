#!/bin/bash

usage(){
    cat <<-EOF
	This script runs a qc check using the 2nd, 3rd and 4th columns of a targets file for RUN, POSITION and TAG values correspondingly.
	npg_qc uses the st::api::lims library to retrieve LIMS information (such as the reference genome); if different LIMS information
	is required, you may want to use a custome sample sheet and direct the library to it by setting the NPG_CACHED_SAMPLESHEET_FILE env var.

	Usage:

	$0 [options] targets_file.txt

	Options:
	   -0                  Dry-run: only print commands to be executed - verbose mode by default.
	   -c <check>          QC check to run.
	   -h                  Show usage message.
	   -i <dir>            Input directory; default is ./input/<RUN>.
	   -n <number>         Line number in targets to be processed
	   -o <dir>            Output directory; default is ./output/<METHOD>/<RUN>/<ID_RUN>/qc/. The check may create its own.
	   -v                  (More) verbose.
	   -w <path>           Absolute path to working directory. Default: $PWD.
	   -x <extra qc args>  Extra arguments passed to qc in a quoted string.
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

while getopts ":0c:hi:n:o:vw:x:" OPTION; do
    case $OPTION in
        0)
            DRYRUN=1;;
        c)
            QC_CHECK=$OPTARG;;
        h)
            usage; exit 1;;
        i)
            IDIR=$OPTARG
            if [[ $IDIR != . ]]; then
                [ -d $IDIR ] || exitmessage "[ERROR] -o: Cannot access input directory ${SDIR}: No such directory" 2
            fi;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${RECNO}" 1;;
        o)
            ODIR=$OPTARG
            if [[ $ODIR != . ]]; then
                if [ -h $ODIR ]; then
                    LINKDIR="$(readlink $ODIR)"
                    [ -d $LINKDIR ] || exitmessage "[ERROR] -o: Cannot access output directory ${LINKDIR}: No such directory" 2
                else
                    [ -d $ODIR ] || exitmessage "[ERROR] -o: Cannot access output directory ${ODIR}: No such directory" 2
                fi
            fi;;
        v)
            VERBOSE=1;;
        w)
            CWD=$OPTARG
            [[ ! -d $CWD ]] && exitmessage "[ERROR] -w: Cannot access ${CWD}: no such directory" 2
            [[ ! $CWD = /* ]] && exitmessage "[ERROR] -w: Not an absolute path" 1;;
        x)
            EXTRAOPTS=$OPTARG;;
        \?)
            exitmessage "[ERROR] Invalid option: -$OPTARG" 1;;
        :)
            exitmessage "[ERROR] Option -$OPTARG requires an argument." 1;;
    esac
done

shift $((OPTIND-1))

TARGETSFILE=$1
DRYRUN=${DRYRUN-0}
[ "$DRYRUN" -eq 1 ] && VERBOSE=1 || VERBOSE=${VERBOSE-0}
[ "$DRYRUN" -eq "1" ] && printf "=========\n[DRY-RUN]\n=========\n" && VERBOSE=1

if env | grep -q ^QC_PATH=; then
    # use whatever version of npg_qc defined by:
    # export QC_PATH=/nfs/users/nfs_r/rb11/dev/perl/rube14-npg_qc
    echo "[INFO] Using QC_PATH=$QC_PATH"
    BINARY="${QC_PATH}/bin/qc"
else
    BINARY=$(readlink -f `which qc`)
    export QC_PATH=${BINARY%/bin*}
fi

QC_EXEC=$BINARY

WORKINGDIR=${CWD:-"$PWD"}
[ "$VERBOSE" -eq 1 ] && printf -- "[INFO] Working directory: ${WORKINGDIR}\n"

# Read info from targets file
if [ -n "$TARGETSFILE" ]; then
    [ ! -f "$TARGETSFILE" ] && exitmessage "[ERROR] Cannot access targets file ${TARGETSFILE}: No such file or directory" 2
    [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1
    LINE=${LSB_JOBINDEX:-$RECNO}
    # Read info from targets file
    RUN=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
    POS=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
    TAG=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
    # deal with non-multiplexed lanes
    if [ -z $TAG ]; then
        BAMID=$RUN\_$POS
    else
        BAMID=$RUN\_$POS\#$TAG
        TAG_INDEX_ARG="--tag_index $TAG"
    fi
fi

QC_OUT=${ODIR:-"${WORKINGDIR}/output/qc"}
QC_IN=${IDIR:-"${WORKINGDIR}/input"}

[ "$VERBOSE" -eq 1 ] && printf "[INFO] Output directory: %s\n" "${QC_OUT}"
[ "$VERBOSE" -eq 1 ] && printf "[INFO] Input directory: %s\n" "${QC_IN}"
[ "$VERBOSE" -eq 1 ] && [ ! -z "$EXTRAOPTS" ] && printf -- "[INFO] Using extra arguments [ ${EXTRAOPTS} ]\n"

QC_CMD="${QC_EXEC} --check ${QC_CHECK} --id_run ${RUN} --position ${POS} ${TAG_INDEX_ARG} --qc_in ${QC_IN} --qc_out ${QC_OUT} ${EXTRAOPTS}"

[ "$VERBOSE" -eq 1 ] && printf -- "[INFO] Running qc command:\n"
[ "$VERBOSE" -eq 1 ] && printf -- "[COMMAND] ${QC_CMD}\n"

[ $DRYRUN -eq 0 ] && QC="$($QC_CMD 2>&1)"

RET_CODE=$?

if [ "$RET_CODE" -eq "0" ]; then
    exitmessage "[INFO] Done QC Check '${QC_CHECK}' for [${BAMID}] = [OK]" 0
else
    exitmessage "[ERROR] QC check for [${BAMID}] exited with exit code ${RET_CODE}: $QC" $RET_CODE
fi
