#!/bin/bash

# load srpipe environment
#. /nfs/gapi/users/rb11/viv_reprocessing/bin/load_srpipe_env.sh
. /software/npg/etc/profile.npg
. /software/sanger-samtools-refpath/etc/profile.sanger-samtools-refpath

###≈ßset -x

usage(){
    cat <<-EOF
	This script runs VIV using the information provided by a targets file.
	[M]ethods available: bam2salmon | bwa_mem | bwa_aln | hisat2| hs_split | salmon | star | tophat2 | y_split
	If -m runfolder is used, paths for output and staging directories are taken from the json file and it's assumed they exist.
	
	Usage: 

	$0 -M <METHOD> [options] targets_file.txt
	
	Options:
	   -c <number>      Do not use WTSI composite id (run_position[#tag]), use insted the
	                    contents of this column in the targets file (must be unique).
	   -h               Show usage message.
	   -m <method hint> Shortcut for specific directory structure: <runfolder | reanalysis>. Default: reanalysis.
	   -n <number>      Line number in targets to be processed
	   -t <number>      If -m runfolder is used, numeric part of tmp_XXXXX folder
	   -w <directory>   Absolute path to working directory. Default: $PWD
	   -x <extra args>  Extra arguments passed to viv in a quoted string
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

while getopts ":c:hM:m:n:t:w:x:" OPTION; do
    case $OPTION in
        c)
            TARGETSCOLUMN=$OPTARG
            [[ ! $TARGETSCOLUMN =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -c: not a digit: ${TARGETSCOLUMN}" 1;;
        h)
            usage; exit 1;;
        M)
            METHOD=$OPTARG
            METHODREGEX="^tophat2|star|bwa\_aln|bwa\_mem|hs\_split|y\_split|bam2salmon|salmon|hisat2$"
            [ -z "$METHOD" ] && exitmessage "[ERROR] -M: a method is required: try '$0 -h' for more information" 1
            [ -n "$METHOD" ] && [[ ! $METHOD =~ $METHODREGEX ]] && exitmessage "[ERROR] -M: invalid method $METHOD: try '$0 -h' for more information" 1;;
        m)
            METHODHINT=$OPTARG
            METHODHINTREGEX="^runfolder|reanalysis$"
            [[ ! $METHODHINT =~ $METHODHINTREGEX ]] && exitmessage "[ERROR] -m: invalid method hint $METHODHINT: try '$0 -h' for more information" 1;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${RECNO}" 1;;
        t)
            TMPDIRNUM=$OPTARG
            [[ ! $TMPDIRNUM =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: not a digit: ${TMPDIRNUM}" 1;;
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

if env | grep -q ^P4_PATH=
then
    # use whatever version of p4 you defined like this:
    # export P4_PATH=/nfs/users/nfs_r/rb11/dev/perl/wtsi-npg_p4
    printf -- "[INFO] Using P4_PATH=$P4_PATH\n"
    BINARY="${P4_PATH}/bin/viv.pl"
else
    # for production stuff e.g. remapping, use production versions
    # for development/experimental/non-deployed stuff e.g. bam2cram use own repo (see path above)
    if [ $METHOD = "bam2cram" ]; then
        exitmessage "[ERROR] Not P4_PATH env variable: try 'export P4_PATH=/path/to/p4'" 1
    else
        BINARY=$(readlink -f `which viv.pl`)
        export P4_PATH=${BINARY%/bin*}
    fi
fi

VIVEXECUTABLE="$BINARY"
WORKINGDIR="${CWD:-"$PWD"}"
METHODHINT="${METHODHINT:-"reanalysis"}"

printf -- "[INFO] Working directory: ${WORKINGDIR}\n"

if [ -n "$TARGETSFILE" ]; then
    [ ! -f "$TARGETSFILE" ] && exitmessage "[ERROR] Cannot access ${TARGETSFILE}: No such file or directory" 2
    [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1
    LINE=${LSB_JOBINDEX:-$RECNO}
    # Read info from targets file
    if [ ! -z $TARGETSCOLUMN ]; then
        BAMID=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk -v column=$TARGETSCOLUMN -F'\t' '{print $column}'`
    else
        RUN=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
        POSITION=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
        TAG=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
        # deal with non-multiplexed lanes
        if [ -z $TAG ]; then
            BAMID="${RUN}_${POSITION}"
        else
            BAMID="${RUN}_${POSITION}#${TAG}"
        fi
    fi
    if [ $METHOD = "bam2cram" ]; then
        SUBMETHOD=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $14}'`
    fi
fi
        
LOGFILE="${BAMID}_viv_${METHOD}.log"
JSONFILE="${BAMID}_${METHOD}.json"

if [ "$METHODHINT" = "reanalysis" ]; then
    OUTDATADIR="${WORKINGDIR}/output/${METHOD}/${RUN}/${BAMID}"
    STAGINGDIR="${WORKINGDIR}/staging/${METHOD}/${RUN}/${BAMID}"
    # create staging and output (including the qc dir) directories
    mkdir -pv "${STAGINGDIR}/tmpdir"
    mkdir -pv "${OUTDATADIR}/qc"
    IJSONDIR="${WORKINGDIR}/json"
elif [ "$METHODHINT" = "runfolder" ]; then
    OUTDATADIR="${WORKINGDIR}/no_cal/archive/lane${POSITION}"
    if [ -z "$TMPDIRNUM" ]; then
        exitmessage "[ERROR] -t: a numeric value is required for -t when -m runfolder is being used (numbers in tmp_XXXXXX directory)"
    else
        STAGINGDIR="${WORKINGDIR}/no_cal/archive/tmp_${TMPDIRNUM}/${BAMID}"
        IJSONDIR="${WORKINGDIR}/no_cal/archive/tmp_${TMPDIRNUM}/${BAMID}"
    fi
fi

printf "[INFO] Output directory: %s\n" "${OUTDATADIR}"

printf -- "[INFO] Changing to staging directory\n"

. /nfs/users/nfs_r/rb11/local/bin/chdir $STAGINGDIR

printf "[INFO] Staging directory: %s\n" "${PWD}"

[ ! -z "$EXTRAOPTS" ] && printf -- "[INFO] Using extra arguments [ ${EXTRAOPTS} ]\n"

VIVCMD="${VIVEXECUTABLE} -s -v 3 -x -o ${LOGFILE} ${EXTRAOPTS} ${IJSONDIR}/${JSONFILE}"

printf -- "[INFO] Running VIV command:\n"
printf -- "[COMMAND] ${VIVCMD}\n"

VIV="$($VIVCMD 2>&1)"
RET_CODE=$?

if [ "$RET_CODE" -eq 0 ]; then
    exitmessage "[INFO] Done VIV Method '${METHOD}' for [${BAMID}] = [OK]" 0
else
    exitmessage "[ERROR] Exited with exit code ${RET_CODE}: ${VIV}" $RET_CODE
fi
