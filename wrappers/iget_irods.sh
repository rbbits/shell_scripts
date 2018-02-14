#!/bin/bash

usage(){
    cat <<-EOF
	This script downloads bam|cram files from iRODS
	Use $0 -H for more details about options.
	
	SYNOPSYS
	
	$0  -i  RUN_POSITION[#TAG[_SPLIT]]      [-0]            [-d VALUE] [-k VALUE]            [-o VALUE] [-R VALUE] [-u VALUE] [-v] [-x VALUE]
	$0  -m  MERGED_TARGETS_FILE             [-0]            [-d VALUE] [-k VALUE] [-n VALUE]            [-R VALUE] [-u VALUE] [-v] [-x VALUE]
	$0  -r  RUN[[_POSITION][#TAG][_SPLIT]]  [-0]            [-d VALUE] [-k VALUE]                                  [-u VALUE] [-v]
	$0  -s  STUDYID                         [-0]            [-d VALUE] [-k VALUE]                                  [-u VALUE] [-v]
	$0  -t  TARGETS_FILE                    [-0] [-c VALUE] [-d VALUE] [-k VALUE] [-n VALUE] [-o VALUE] [-R VALUE] [-u VALUE] [-v] [-x VALUE]
	
	EXAMPLES

	$0 -t targets.txt        # a targets file containing specific id's and other info (e.g. save dir)
	$0 -m merged_targets.txt # a targets file containing merged id's: run and tag, or library id columns and other info
	$0 -i 12345_1#1_yhuman   # specific target
	$0 -r 12345              # everything for run 12345
	$0 -r 12345_3            # everything in lane 3 in run 12345
	$0 -r 12345_#2           # only tag index 2 of every lane in run 12345
	$0 -r 12345_#_yhuman     # only yhuman alignment filter for every tag index and every lane in run 12345

	EOF
}


morehelp(){
    cat <<-EOF
	usage: $0 options
	
	OPTIONS
	
	Required:
	
	   -i  RUN_POSITION[#TAG[_SPLIT]]  |  -m  MERGED_TARGETS_FILE  |  -r  RUN[[_POSITION][#TAG][_SPLIT]]  |  -s  STUDYID  |  -t  TARGETS_FILE
	
	Options:
	
	   -0               Dry-run: only print commands to be executed - verbose mode by default.
	   -c <number  >    Use to indicate the column in the targets file to be used as source of ids.
	   -d <savedir>     Download directory; default is ./<RUN>; If -t is used you can use -c <column n>.
	   -h               Show usage message.
	   -H               Show this message.
	   -i <run_id>      A full valid run id (e.g. 17550_1#3, 16660_3, 17550_8#0_phix), if it doesn't exist in iRODS an error is thrown.
	   -k <full path>   Full path to keytabs directory. Use env KEYTABPATH for default path otherwise one is required.
	   -m <file>        Path to targets file with records grouped by library_id (default) or run_id and tag_index. The following columns are assumed to contain the id values: 2nd = run or library id, 3rd = tag index. Use -x to specify which case should be assumed.
	   -n <number>      Record (line number) inside targets file to be processed.
	   -o <out format>  Download output in this format.
	                      -o b        BAM file, fail if doest't exist or can't be accessed.
	                      -o c        CRAM file. This is the default.
	                      -o f        Files in FASTQ format (.fq.gz).	
	                      -o i        Download index file only (.cram.crai or .bai) if both files are wanted then use -o c|b|B and -x i.
	                      -o t        Tophat files only, extracting accepted_hits.bam from CRAM or BAM plus other files as specified by option -o.
	                      -o B        BAM file, if only CRAM exist convert CRAM->BAM.
	   -p <path>        Absolute path to repository for reference genome/transcriptome
	   -r <regex>       Regular expression must contain at least one of the following components in their expected order to work: <RUN><POSITION><TAG><SPLIT>.
	                      <RUN>       Numeric. A run id (e.g. 17550, 80001).
	                      <POSITION>  Alphaumeric. A '_' followed by one number from 1 to 8 (e.g. _1, _2, _3).
	                      <TAG>       Alphanumeric. A '#' followed by a number (e.g. #1, #888, #0).
	                      <SPLIT>     Alpha. A '_' followed by a string denoting a split (e.g. human, phix, y_human).
	   -R <number|file> Use with -o = f if BAM is not present in iRODS and REF_PATH is not present in the env (you may want to use sanger-samtools-refpath for that).
	                    Use with -m to pass a reference genome to bammerge as an argument.
	                      -R <number> Column in targets file that contains the path to the reference genome file.
	                      -R <file>   Path to reference genome file.
	   -t <file>        Path to targets file. The following columns are assumed to contain the id values: 2nd = run_id, 3rd = position, 4th = tag_index. Use -c to specify a different column.
	   -u <auth mode>   IRODS user authentication method.
	                      -u k        Kerberos credentials (provide path to keytab with -k if necessary).
	                      -u p        Password authentication (iinit), this is the default.
	   -v               (More) verbose.
	   -x <options>     Specify extra files to be downloaded or extra options to be used. 
	                    Use with -o = t to download (in all cases files will be saved into ./tophat/<runid>/<rpt> and -d is ignored).
	                      -x a        All available Tophat files: accepted_hits.bam, unmapped.bam, *.bed.
	                      -x b        Tophat's BAMs: accepted_hits.bam and unmapped.bam.
	                      -x h        Tophat's accepted_hits.bam only (default).
	                    Use with -o = c|b|B to download along with BAM/CRAM.
	                      -x a        All target tags if available (No PhiX/split files).
	                      -x i        Index file: .cram.crai or .bai correspondingly.
	                    Use with -m to indicate how many columns to use as merging values (starting from column 2).
	                      -x 1        CRAM files merged by library id (column 2 only) - this is the default.
	                      -x 2        CRAM files merged by run_id and tag_index (columns 2 and 3 respectively).
	
	AUTHOR
	
	Ruben Bautista <rb11@sanger.ac.uk>.
	EOF
}

# function to print messages
# to STDOUT or STDERR
exitmessage(){
    declare EXITMESG=$1
    declare EXITCODE=$2
    if [ $EXITCODE -eq 0 ]; then
        >&1 printf '%b\n' "$EXITMESG"
    else
        >&2 printf '%b\n' "$EXITMESG"
    fi
    exit $EXITCODE
}

#set -x

#TODO: Not implemented yet: <s>tats, seqch<k>sum
while getopts ":0c:d:hHi:k:m:n:o:p:r:R:s:t:u:vx:" OPTION; do
    case $OPTION in
        0)
            DRYRUN=1;;
        c)
            TARGETSCOLUMN=$OPTARG
            [[ ! $TARGETSCOLUMN =~ ^[[:digit:]]?$ ]] && exitmessage "[ERROR] -c: targets column not a number: ${TARGETSCOLUMN}" 1;;
        d)
            SDIR=$OPTARG
            if [[ $SDIR != . ]]; then
                [ -d $SDIR ] || exitmessage "[ERROR] -o: Cannot access ${SDIR}: No such directory" 2
            fi;;
        h)
            usage; exit 1;;
        H)
            morehelp; exit 1;;
        i)
            ITARGET=$OPTARG
            IOPT=1
            INPUTMODE="i";;
        k)
            KEYTABPATH=$OPTARG;;
        m)
            MTARGET=$OPTARG
            [ ! -f "$MTARGET" ] && exitmessage "[ERROR] -m: Cannot access ${MTARGET}: No such file or directory" 2
            MOPT=1
            INPUTMODE="m"
            ;;
        n)
            RECNO=$OPTARG
            [[ ! $RECNO =~ ^[0-9]+$ ]] && exitmessage "[ERROR] -n: targets row not a number: ${RECNO}" 1;;
        o)
            OUTPUTFORMAT=$OPTARG
            [[ ! $OUTPUTFORMAT =~ ^[bBcfit]$ ]] && exitmessage "[ERROR] -o: output format is invalid: ${OUTPUTFORMAT}. Try ${0} -H for help" 1;;
        p)
            REPOSITORY=$OPTARG
            [[ ! -d $REPOSITORY ]] && exitmessage "[ERROR] -p: Cannot access repository ${REPOSITORY}: no such directory" 2
            [[ ! $REPOSITORY = /* ]] && exitmessage "[ERROR] -p: Not an absolute path" 1
            ;;
        r)
            RTARGET=$OPTARG
            ROPT=1
            INPUTMODE="r";;
        R)
            if [[ ! $OPTARG =~ ^[0-9]+$ ]]; then
                [ ! -f "$OPTARG" ] && exitmessage "[ERROR] -R: Cannot access ${REFERENCE}: No such file or directory" 2
                REFERENCE=$OPTARG
            else
                REF_TARGETSCOLUMN=$OPTARG
            fi;;
        s)
            STARGET=$OPTARG
            SOPT=1
            INPUTMODE="s";;
        t)
            TTARGET=$OPTARG
            [ ! -f "$TTARGET" ] && exitmessage "[ERROR] -t: Cannot access ${TTARGET}: No such file or directory" 2
            TOPT=1
            INPUTMODE="t"
            ;;
        u)
            IRODSUSERAUTH=$OPTARG
            [[ ! $IRODSUSERAUTH =~ ^[kp]$ ]] && exitmessage "[ERROR] -u: iRODS user authentication method is invalid: ${IRODSUSERAUTH}. Try ${0} -H for help" 1;;
        v)
            VERBOSE=1;;
        x)
            FMTXTRAOPT=$OPTARG
            [[ ! $FMTXTRAOPT =~ ^[12abhi]$ ]] && exitmessage "[ERROR] -x: extra option is invalid: ${FMTXTRAOPT}. Try ${0} -H for help" 1;;
        \?)
            echo "Invalid option: -$OPTARG" >&2; exit 1;;
        :)
            echo "Option -$OPTARG requires an argument." >&2; exit 1;;
    esac
done


if [[ ! $INPUTMODE =~ [imrst] ]]; then
    exitmessage "$0: at least one target option must be used: -i|m|r|s|t" 1
else
    let MULTIOPTS=$IOPT+$MOPT+$ROPT+$SOPT+$TOPT+0
    [ $MULTIOPTS -gt 1 ] && exitmessage "$0: only one target option can be used: -i|m|r|s|t" 1
fi


DRYRUN=${DRYRUN-0}
VERBOSE=${VERBOSE-0}
DEFAULTDIR="."
OUTPUTFORMAT=${OUTPUTFORMAT-c}
REPOSDIR=${REPOSITORY-"/lustre/scratch117/core/sciops_repository"}


if [ "$OUTPUTFORMAT" = "t" ]; then
    FMTXTRAOPT=${FMTXTRAOPT-h}
    [[ ! $FMTXTRAOPT =~ ^[abh]$ ]] && exitmessage "[ERROR] -x: Wrong extra option for output format: -o ${OUTPUTFORMAT} -x ${FMTXTRAOPT}. Try ${0} -H for help" 1
elif [[ $OUTPUTFORMAT =~ ^[cb]$ ]]; then
    if [ ! -z "$FMTXTRAOPT" ]; then
        [[ ! $FMTXTRAOPT =~ ^[ai12]$ ]] && exitmessage "[ERROR] -x: Wrong extra option for output format: -o ${OUTPUTFORMAT} -x ${FMTXTRAOPT}. Try ${0} -H for help" 1
    fi
fi


IMETA_EXE=`which imeta`
[ "$?" -ne "0" ] && exitmessage "[ERROR] imeta: command not found" 1
IGET_BIN=`which iget`
[ "$?" -ne "0" ] && exitmessage "[ERROR] iget: command not found" 1
BAMSORT_BIN=`which bamsort`
[ "$?" -ne "0" ] && printf -- "[WARNING] bamsort: command not found\n"
BAM2FQ_BIN=`which bamtofastq`
[ "$?" -ne "0" ] && printf -- "[WARNING] bamtofastq: command not found\n"
SAMTOOLS_BIN=`which samtools`
[ "$?" -ne "0" ] && printf -- "[WARNING] samtools: command not found\n"
SCRAMBLE_BIN=`which scramble`
[ "$?" -ne "0" ] && printf -- "[WARNING] scramble: command not found\n"
BATON_BIN=`which baton`
[ "$?" -ne "0" ] && printf -- "[WARNING] baton: command not found\n"
BATON_METAQUERY_BIN=`which baton-metaquery`
[ "$?" -ne "0" ] && printf -- "[WARNING] baton-metaquery: command not found\n"
BAMMERGE_BIN=`which bammerge`
[ "$?" -ne "0" ] && printf -- "[WARNING] bammerge: command not found\n"
JQ_BIN=`which jq`
[ "$?" -ne "0" ] && printf -- "[WARNING] jq: command not found\n"


USER=`whoami`
[ "$USER" = "srpipe" ] && umask 0002
IRODSUSERAUTH=${IRODSUSERAUTH-u}
if [ "$IRODSUSERAUTH" = "k" ]; then
    if env | grep -q ^KEYTABPATH=; then
        [ "$VERBOSE" -eq 1 ] && echo "[INFO] Using keytabs found in $KEYTABPATH"
        KEYTABFILE=${USER}.keytab
        if [ -d "$KEYTABPATH" ] && [ -e "${KEYTABPATH}/${KEYTABFILE}" ]; then
            [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] /usr/bin/kinit  ${USER}\@INTERNAL.SANGER.AC.UK -k -t ${KEYTABPATH}/${KEYTABFILE}"
            /usr/bin/kinit  ${USER}\@INTERNAL.SANGER.AC.UK -k -t  "${KEYTABPATH}/${KEYTABFILE}"
        else
            exitmessage "-k: No keytab found: cannot access ${KEYTABPATH}/${KEYTABFILE}: no such file or directory. Try export KEYTABPATH=/path/to/keytabs?" 1
        fi
    else
        if [ -z "$KEYTABPATH" ]; then
            exitmessage "-k: No keytab found: cannot access ${KEYTABPATH}/${KEYTABFILE}: no such file or directory. Try KEYTABPATH=/path/to/keytabs?" 1
        fi
    fi
else
    [ "$VERBOSE" -eq 1 ] && echo "[INFO] Using password authentication"
    #[ "$VERBOSE" -eq 1 ] && echo "[INFO] If download fails due to authentication run \`iinit' then try again"
    #[ "$VERBOSE" -eq 1 ] && echo "[INFO] If problems persist contact helpdesk for help logging in to iRODS"
fi


case "$INPUTMODE" in

    i|t)
        if [ "$INPUTMODE" = "i" ]; then
            XAMID="${ITARGET}"
        elif [ "$INPUTMODE" = "t" ]; then
            [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1
            LINE=${LSB_JOBINDEX:-$RECNO}
            if [ ! -z "$TARGETSCOLUMN" ]; then
                XAMID=`head -n ${LINE} ${TTARGET} | tail -1 | awk -v col=$TARGETSCOLUMN 'BEGIN { FS = "\t" } ; {print $col}'`
                [ "$VERBOSE" -eq 1 ] && echo "[INFO] Using column ${TARGETSCOLUMN} to find id: ${XAMID}"
            else
                RUN=`head -n ${LINE} ${TTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
                POS=`head -n ${LINE} ${TTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
                TAG=`head -n ${LINE} ${TTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
                XAMID="${RUN}_${POS}"
                [ -n "$TAG" ] && XAMID="${RUN}_${POS}#${TAG}"
                [ -n "$SPL" ] && XAMID="${RUN}_${POS}#${TAG}_${SPL}"
            fi
        fi
        RPT_REGEX='^([[:digit:]]+)_([[:digit:]]+)(_([[:alpha:]]+))?(#([[:digit:]]+))?(_([[:alpha:]]+))?$'
        if [[ "$XAMID" =~ $RPT_REGEX ]]; then
            RUN="${BASH_REMATCH[1]}";
            POS="${BASH_REMATCH[2]}";
            TAG="${BASH_REMATCH[6]}";
            SPL="${BASH_REMATCH[8]}";
        else
            exitmessage "[ERROR] -i|-t: target '$XAMID' doesn't look right (e.g. 14940_3#34_human) Other formats not supported" 1
        fi
        ;;
    m)
        [ -z "$LSB_JOBINDEX" -a -z "$RECNO" ] &&  exitmessage "[ERROR] -n: Env variable LSB_JOBINDEX not set and no record number was specified in its place" 1
        LINE=${LSB_JOBINDEX:-$RECNO}
        MERGEBYCOLUMNS=${FMTXTRAOPT-1}
        if [ "$MERGEBYCOLUMNS" -eq 1 ]; then
            [ "$VERBOSE" -eq 1 ] && printf -- "[INFO] Using column 2 - library_id as search and merging value\n"
            LIBRARY=`head -n ${LINE} ${MTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
            XAMID="${LIBRARY}"
            BATON_AVU_OPTS="-a library_id -v ${LIBRARY}"
        elif [ "$MERGEBYCOLUMNS" -eq 2 ]; then
            [ "$VERBOSE" -eq 1 ] && printf -- "[INFO] Using columns 2 - run_id and 3 - tag_index as search and merging values\n"
            RUN=`head -n ${LINE} ${MTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
            TAG=`head -n ${LINE} ${MTARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
            XAMID="${RUN}_${TAG}"
            BATON_AVU_OPTS="-a id_run -v ${RUN} -a tag_index -v ${TAG}"
        else
            exitmessage "[ERROR] -x: Wrong extra option for input mode 'm': supported values are 1|2. Try ${0} -H for help" 1
        fi
        if [ -n "$REF_TARGETSCOLUMN" ]; then
            REFERENCE=`head -n ${LINE} ${MTARGET} | tail -1 | awk -v col=$REF_TARGETSCOLUMN 'BEGIN { FS = "\t" } ; {print $col}'`
            [[ ! $REFERENCE = /* ]] && REFERENCE="${REPOSDIR}/references/${REFERENCE}"
            [ "$VERBOSE" -eq 1 ] && printf "[INFO] Using reference from targets: %s\n" "$REFERENCE"
            [ ! -f "$REFERENCE" ] && exitmessage "[ERROR] -R: Cannot access reference ${REFERENCE}: No such file or directory" 2
        fi
        ;;
    r)
        regex='^([[:digit:]]+)(_([[:digit:]]*))?(#([[:digit:]]*))?(_([[:alpha:]]*))?'
        if [[ "$RTARGET" =~ $regex ]]; then
            RUN="${BASH_REMATCH[1]}";
            POS="${BASH_REMATCH[3]}";
            TAG="${BASH_REMATCH[5]}";
            SPL="${BASH_REMATCH[7]}";
        else
            exitmessage "[ERROR] -r: target '$RTARGET' doesn't look right" 1
        fi
        [ ! -n "$RUN" ] && exitmessage "[ERROR] -i: At least a run id must be provided" 1
        [ -n "$POS" ] && QU_CMD+=" and lane "'>='" $POS and lane "'<='" $POS "
        [ -n "$TAG" ] && QU_CMD+=" and tag_index = $TAG "
        [ -n "$SPL" ] && QU_CMD+=" and alignment_filter = $SPL "
        META_QU_CMD="$IMETA_EXE qu -z seq -d id_run = $RUN $QU_CMD and target "'>='" 1 "
        ;;

    s)
        [ ! -n "$STARGET" ] && exitmessage "[ERROR] -s: A study id must be provided" 1
        [[ ! "$STARGET" =~ ^[[:digit:]]+$ ]] && exitmessage "-s: target study id '$STARGET' doesn't look right" 1
        META_QU_CMD="$IMETA_EXE qu -z seq -d study_id = $STARGET and target = 1 "
        ;;
esac


if [[ $OUTPUTFORMAT =~ [bB] ]]; then
    FTYPE_CMD=' and type = bam '
else #c,f,i,t
    FTYPE_CMD=' and type = cram '
fi
META_QU_CMD+="$FTYPE_CMD"

[ "$DRYRUN" -eq "1" ] && printf "=========\n[DRY-RUN]\n=========\n" && VERBOSE=1

SDIR=${SDIR-$DEFAULTDIR}
[ "$OUTPUTFORMAT" = "t" ] && SDIR+="/tophat"
CMD="mkdir -pv ${SDIR}"
[ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
[ $DRYRUN -eq 0 ] && mkdir -pv ${SDIR}

RET_CODE=1
INFO_MESSAGE=""

case "$INPUTMODE" in
    
    i|t)
        
        LSBAM="$(ils /seq/${RUN}/${XAMID}.bam 2>&1)"
        NOBAM=$?
        case $NOBAM in
            0)
                [[ $OUTPUTFORMAT =~ [bB] ]] && [ "$VERBOSE" -eq 1 ] && printf -- "[INFO] BAM file FOUND: /seq/${RUN}/${XAMID}.bam\n"
                [ "$OUTPUTFORMAT" = "B" ] && OUTPUTFORMAT="b"
                XAM_FILE_EXT="bam"
                XAM_IXFILE_EXT="bai";;
            4)
                [[ $OUTPUTFORMAT =~ [bB] ]] && [ "$VERBOSE" -eq 1 ] && printf -- "[INFO] BAM file NOT FOUND: /seq/${RUN}/${XAMID}.bam\n" && printf -- "[INFO] Exit message: ${LSBAM}\n"
                XAM_FILE_EXT="cram"
                XAM_IXFILE_EXT="cram.crai";;
            3)
                [[ $OUTPUTFORMAT =~ [bB] ]] && [ "$VERBOSE" -eq 1 ] && printf -- "[ERROR] BAM file NOT ACCESSIBLE: /seq/${RUN}/${XAMID}.bam\n"
                exitmessage "[ERROR] ${LSBAM}" $NOBAM;;
            *)
                [[ $OUTPUTFORMAT =~ [bB] ]] && [ "$VERBOSE" -eq 1 ] && printf -- "[ERROR] Error trying to get bam file for ${XAMID}\n"
                exitmessage "[ERROR] Error trying to get bam file for [${XAMID}] with exit code ${NOBAM}: ${LSBAM}" $NOBAM;;
        esac

        RET_CODE=0
        case "$OUTPUTFORMAT" in
            b|c)
                CMD="${IGET_BIN} -vf /seq/${RUN}/${XAMID}.${XAM_FILE_EXT} ${SDIR}/${XAMID}.${XAM_FILE_EXT}"
                [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" " $CMD"
                [ $DRYRUN -eq 0 ] && IGETCMD="$($CMD 2>&1)" && RET_CODE=$? && INFO_MESSAGE="[INFO] [IGET] ${IGETCMD}\n"
                if [ "$FMTXTRAOPT" = "i" ]; then
                    RET_CODE_IDX=0
                    CMD="${IGET_BIN} -KPvf /seq/${RUN}/${XAMID}.${XAM_IXFILE_EXT} ${SDIR}/${XAMID}.${XAM_IXFILE_EXT}"
                    [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                    [ $DRYRUN -eq 0 ] && IGETIXCMD="$($CMD 2>&1)" && RET_CODE_IDX=$? && INFO_MESSAGE+="[INFO] [IGET] ${IGETIXCMD}\n"
                    if [ "$RET_CODE_IDX" -ne "0" ]; then
                        printf -- "[INFO] Failed to download index file for [${SDIR}/${XAMID}.${XAM_FILE_EXT}]:\n"
                        printf -- "[ERROR] $IGETIXCMD\n"
                    fi
                fi
                ;;
            B)
                echo "[INFO] No bam format file present in iRODS for [${XAMID}]: $LSBAM"
                echo "[INFO] Converting cram->bam ..."
                CMD="${SAMTOOLS_BIN} view -bh -o ${SDIR}/${XAMID}.bam irods:/seq/${RUN}/${XAMID}.cram"
                [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] $CMD"
                [ $DRYRUN -eq 0 ] && CRAM2BAMCMD="$($CMD 2>&1)" && RET_CODE=$?
                if [ "$FMTXTRAOPT" = "i" ]; then
                    RET_CODE_IDX=0
                    CMD="${SAMTOOLS_BIN} index ${SDIR}/${XAMID}.bam"
                    [ "$VERBOSE" -eq 1 ] && echo "[INFO] Indexing cram or bam ..."
                    [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] $CMD"
                    [ $DRYRUN -eq 0 ] && IXCMD="$($CMD 2>&1)" && RET_CODE_IDX=$?
                    [ "$RET_CODE_IDX" -ne "0" ] && echo "[WARNING] Failed to create index file for ${SDIR}/${XAMID}.bam:" && echo "[ERROR] $IGETIXCMD"
                fi
                ;;
            f)
                REF_OPT=""
                if [ "$NOBAM" -eq "4" ]; then
                    if ! printenv | grep -q ^REF_PATH=; then
                        if [ ! -z "$REFERENCE" ]; then
                            REF_OPT="reference=$REFERENCE"
                        else
                            exitmessage "[ERROR] -f: A reference genome is required: use -R" 1
                        fi
                    fi
                fi
                FQ1="${SDIR}/${XAMID}_1.fq.gz"
                FQ2="${SDIR}/${XAMID}_2.fq.gz"
                CMD="${IGET_BIN} /seq/${RUN}/${XAMID}.${XAM_FILE_EXT} - | ${BAMSORT_BIN} SO=queryname level=0 inputformat=${XAM_FILE_EXT} ${REF_OPT} outputformat=bam index=0 tmpfile=./bamfq/tmp_${XAMID} | ${BAM2FQ_BIN} inputformat=bam exclude=QCFAIL T=./bamfq/tmp_b2fq_${XAMID} gz=1 level=9 F=${FQ1} F2=${FQ2}"
                [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                [ $DRYRUN -eq 0 ] && IGETCMD="$(iget /seq/${RUN}/${XAMID}.${XAM_FILE_EXT} - | bamsort SO=queryname level=0 inputformat=${XAM_FILE_EXT} ${REF_OPT} outputformat=bam index=0 tmpfile=./bamfq/tmp_${XAMID} | bamtofastq inputformat=bam exclude=QCFAIL T=./bamfq/tmp_b2fq_${XAMID} gz=1 level=9 F=${FQ1} F2=${FQ2} 2>&1)" && RET_CODE=$?
                ;;
            i)
                CMD="${IGET_BIN} -KPvf /seq/${RUN}/${XAMID}.${XAM_IXFILE_EXT} ${SDIR}/${XAMID}.${XAM_IXFILE_EXT}"
                [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                [ $DRYRUN -eq 0 ] && IGETCMD="$($CMD 2>&1)" && RET_CODE=$?
                ;;
            t)
                ret_code_junc=0
                ret_code_dels=0
                ret_code_hits=0
                ret_code_unmp=0
                CMD="${SAMTOOLS_BIN} view -bh -F 0x4 -o ${SDIR}/${XAMID}_accepted_hits.bam irods:/seq/${RUN}/${XAMID}.${XAM_FILE_EXT}"
                [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                [ $DRYRUN -eq 0 ] && ACCEPTEDHITSFILE="$($CMD 2>&1)" && ret_code_hits=$?
                if [[ $FMTXTRAOPT =~ [ab] ]] ; then
                    CMD="${SAMTOOLS_BIN} view -bh -f 0x4 -o ${SDIR}/${XAMID}_unmapped.bam irods:/seq/${RUN}/${XAMID}.${XAM_FILE_EXT}"
                    [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                    [ $DRYRUN -eq 0 ] && UNMAPPEDFILE="$($CMD 2>&1)" && ret_code_unmp=$?
                fi
                if [ "$FMTXTRAOPT" = "a" ]; then
                    CMD="${IGET_BIN} -KPvf /seq/${RUN}/${XAMID}.junctions.bed  ${SDIR}/${XAMID}.junctions.bed"
                    [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] ${CMD}"
                    [ $DRYRUN -eq 0 ] && JUNCTIONSFILE="$($CMD 2>&1)" && ret_code_junc=$?
                    CMD="${IGET_BIN} -KPvf /seq/${RUN}/${XAMID}.deletions.bed  ${SDIR}/${XAMID}.deletions.bed"
                    [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] ${CMD}"
                    [ $DRYRUN -eq 0 ] && DELETIONSFILE="$($CMD 2>&1)" && ret_code_dels=$?
                fi
                let "RET_CODE = ret_code_junc + ret_code_dels + ret_code_hits + ret_code_unmp"
                if [ "$RET_CODE" -ne "0" ]; then
                    [ "$VERBOSE" -eq 1 ] && echo "[WARNING] -f t: One or more files failed to download"
                    [ "$ret_code_junc" -ne "0" ] && echo "[WARNING] Failed to download junctions:" && echo "  $JUNCTIONSFILE"
                    [ "$ret_code_dels" -ne "0" ] && echo "[WARNING] Failed to download deletions:" && echo "  $DELETIONSFILE"
                    [ "$ret_code_hits" -ne "0" ] && echo "[WARNING] Failed to download accepted hits:" && echo "  $ACCEPTEDHITSFILE"
                    [ "$ret_code_unmp" -ne "0" ] && echo "[WARNING] Failed to download unmapped:" && echo "  $UNMAPPEDFILE"
                fi
                ;;
        esac
        ;;
    m)
        BATON_COMMON_OPTS="-a type -v cram -a target -v 1"
        BATON_CMD="${BATON_BIN} ${BATON_COMMON_OPTS} ${BATON_AVU_OPTS} | ${BATON_METAQUERY_BIN} --zone seq --avu | ${JQ_BIN} -r '.[] | .collection  + \"/\" + .data_object'"
        [ "$VERBOSE" -eq 1 ] && printf -- "[COMMAND] $BATON_CMD\n"
        if [ "$DRYRUN" -eq 0 ]; then
            META_INFO=$(${BATON_BIN} ${BATON_COMMON_OPTS} ${BATON_AVU_OPTS} | ${BATON_METAQUERY_BIN} --zone seq --avu | ${JQ_BIN} -r '.[] | .collection  + "/" + .data_object')
        else
            META_INFO=`printf "/seq/${RUN}/${RUN}_1#${TAG}.cram\n/seq/${RUN}/${RUN}_2#${TAG}.cram\n"`
        fi
        if [ -n "$REFERENCE" ]; then
            REF_OPT="reference=$REFERENCE"
        fi
        if [ "$?" -eq "0" ]; then
            HIT_LIST=($META_INFO)
            if [ "${#HIT_LIST[@]}" -ge "1" ]; then
                [ "$VERBOSE" -eq 1 ] && printf "[INFO] %s files found\n" "${#HIT_LIST[@]}"
                BAMMERGE_INPUTS=`printf 'I=irods:%s ' ${META_INFO[@]}`
                CMD="${BAMMERGE_BIN} level=1 inputformat=cram outputformat=cram tmpfile=staging/${RUN}_${TAG}_tmpfile ${BAMMERGE_INPUTS} ${REF_OPT} > ${SDIR}/${XAMID}.cram"
                [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] ${CMD}"
                [ "$DRYRUN" -eq 0 ] && IGETCMD=$(${BAMMERGE_BIN} level=1 inputformat=cram outputformat=cram tmpfile=staging/${RUN}_${TAG}_tmpfile ${BAMMERGE_INPUTS} ${REF_OPT} > ${SDIR}/${XAMID}.cram) && RET_CODE=$?
                
            else
                [ "$VERBOSE" -eq 1 ] && printf -- "[INFO] 1 file found: nothing to merge\n"
                CMD="${IGET_BIN} -KPvf ${META_INFO[@]} ${SDIR}/${XAMID}.cram"
                [ "$VERBOSE" -eq 1 ] && echo "[COMMAND] ${CMD}"
                [ "$DRYRUN" -eq 0 ] && IGETCMD="$($CMD)" && RET_CODE=$?
            fi
        else
            exitmessage "[ERROR] Error trying to query metadata: $META_INFO" 2
        fi
        ;;     
    r)    
        [ "$VERBOSE" -eq 1 ] && printf -- "[COMMAND] $META_QU_CMD\n"
        META_INFO="$($META_QU_CMD)"
        if [ "$?" -eq "0" ]; then
            declare -a HIT_LIST=(`grep -oP "(?<=dataObj: ).*$" <<< "$META_INFO"`)
            declare -i NUM_HITS="${#HIT_LIST[@]}"
            if [ "$NUM_HITS" -ge "1" ]; then
                printf "[INTERACTIVE] Found %d file(s). Do you want to download all/some/none/list? " "$NUM_HITS"
                read ANSWER
                case "$ANSWER" in
                    a*)
                        printf -- "[INFO] Downloading all of the files ...\n"
                        RET_CODE=0
                        for IRODS_FILE in "${HIT_LIST[@]}"; do
                            CMD="${IGET_BIN} -KPvf /seq/${RUN}/${IRODS_FILE} ${SDIR}/${IRODS_FILE}"
                            [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                            ret_code=0
                            [ $DRYRUN -eq 0 ] && IGETCMD="$($CMD 2>&1)" && ret_code=$?
                            if [ "$ret_code" -eq 0 ]; then
                                [ "$VERBOSE" -eq 1 ] && printf "[INFO] Succesfully downloaded %s\n" "$IRODS_FILE"
                            else
                                printf "[INFO] Failed to download file %s with error message:\n" "$IRODS_FILE"
                                printf "[ERROR] %s\n" "${IGETCMD}"
                            fi
                            let "RET_CODE += $ret_code"
                        done;;
                    l*)
                        printf '%s |' "${HIT_LIST[@]}";
                        exitmessage "Bye!" 0;;
                    n*)
                        exitmessage "Bye!" 0;;
                    s*)
                        printf -- "[INFO] Downloading some of the files only. Press [Q] or [q] to quit ...\n"
                        RET_CODE=0
                        for IRODS_FILE in "${HIT_LIST[@]}"; do
                            REPLY=""
                            while [[ ! $REPLY =~ [nNyYqQ] ]]; do
                                printf "[INTERACTIVE] Download %s? [Yy/Nn/Qq] " "$IRODS_FILE"
                                read -n 1 REPLY
                                if [[ $REPLY =~ y|Y ]]; then
                                    CMD="${IGET_BIN} -KPvf /seq/${RUN}/${IRODS_FILE} ${SDIR}/${IRODS_FILE}"
                                    [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "${CMD}"
                                    ret_code=0
                                    [ $DRYRUN -eq 0 ] && IGETCMD="$($CMD 2>&1)" && ret_code=$?
                                    if [ "$ret_code" -eq 0 ]; then
                                        [ "$VERBOSE" -eq 1 ] && printf "[INFO] Succesfully downloaded %s\n" "$IRODS_FILE"
                                    else
                                        printf "[INFO] Failed to download file %s with error message:\n" "$IRODS_FILE"
                                        printf "[ERROR] %s\n" "${IGETCMD}"
                                    fi
                                    let "RET_CODE += $ret_code"
                                fi
                                [[ $REPLY =~ n|N ]] && printf "[INFO] Skipping %s\n" "$IRODS_FILE"
                                [[ $REPLY =~ q|Q ]] && printf -- "[INFO] Stopping downloads. Bye!\n" && break
                            done
                            [[ $REPLY =~ Q|q ]] && break
                        done
                esac
            else
                exitmessage "$META_INFO" 0
            fi
        else
            exitmessage "[ERROR] $META_INFO" 2 #I couldn't find any files with the following information:  run: $RUN  pos: $POS  tag: $TAG  spl: $SPL" 2
        fi
        ;;

    s)
        [ "$VERBOSE" -eq 1 ] && printf -- "[COMMAND] $META_QU_CMD\n"
        META_INFO="$($META_QU_CMD 2>&1)"
        if [ "$?" -eq "0" ]; then
            declare -a HITS_LIST=( `grep -oP "(?<=dataObj: ).*$" <<< "$META_INFO"` )
            declare -a COLL_LIST=( `grep -oP "(?<=collection: ).*$" <<< "$META_INFO"` )
            NUM_HITS=${#HITS_LIST[@]}
            if [ ! "$NUM_HITS" -ge "1" ]; then
                RET_CODE=1
                exitmessage "[WARNING] $META_INFO" $RET_CODE
            else
                [ "$VERBOSE" -eq 1 ] && printf "Found %d file(s)" "$NUM_HITS"
            fi
            [ "$VERBOSE" -eq 1 ] && printf -- "INFO: Downloading all of the files ...\n"
            RET_CODE=0
            let "NUM_HITS-=1"
            for IDX in $(seq  0  $NUM_HITS); do
                IRODS_FILE="${COLL_LIST[$IDX]}/${HITS_LIST[$IDX]}"
                if [[ ${COLL_LIST[$IDX]} =~ ([[:digit:]]+) ]]; then
                    RUN="${BASH_REMATCH[1]}";
                fi
                CMD="${IGET_BIN} -KPvf ${IRODS_FILE} ${SDIR}/${RUN}/${HITS_LIST[$IDX]}"
                [ "$VERBOSE" -eq 1 ] && printf "[COMMAND] %s\n" "$CMD"
                ret_code=0
                [ $DRYRUN -eq 0 ] && IGETCMD="$($CMD 2>&1)" && ret_code=$?
                if [ "$ret_code" -eq 0 ]; then
                    [ "$VERBOSE" -eq 1 ] && printf "[INFO] Succesfully downloaded %s\n\n" "$IRODS_FILE"
                else
                    printf "[INFO] Failed to download file %s with error message:\n" "$IRODS_FILE"
                    printf "[ERROR] %s\n" "$IGETCMD"
                fi
                let "RET_CODE += $ret_code"
            done
        fi
        ;;
esac



if [ "$RET_CODE" -ne 0 ]; then
    exitmessage "[ERROR] Error while trying to download [${XAMID}] Exited with exit code ${RET_CODE}:\n ${IGETCMD}" $RET_CODE
else
    INFO_MESSAGE+="[INFO] Done [${XAMID}]"
    exitmessage "${INFO_MESSAGE}" 0
fi
