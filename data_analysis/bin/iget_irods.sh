#!/bin/bash

usage(){
    cat <<-EOF
	usage: $0 options
	
	This script downloads bam|cram files from iRODS
	
	OPTIONS
	
	Required:

	   -i <rpt|file>    Input id ( RUN_POSITION[#TAG[_SPLIT]] ) or a targets file.
	   -k <full path>   Full path to keytabs directory
	
	Optional:
	
	   -d <savedir>     Download directory; defaults: If -i rpt Then ./irods/<RUN>; If -i targetsfile Then 7th column.
	   -f <format|file> Download format: <b>am, <c>ram, <f>astq (.fq.gz), <t>ophat; default: c.
	                       If -f = b & bam not present in iRODS then cram->bam will be enforced unless -n is used.
	                       If -f = t accepted_hits.bam will be extracted from source, see -t for options.
	   -n               Don't carry out cram->bam operation if there's no bam (exit with error).
	   -r <reference>   Only needed if -f = f & bam not present in iRODS.
	   -t <options>     Only needed if -f = t: <a>ll, <b>ams, <h>its; default: h. Downloads:
	                       <a>ll available files: accepted_hits.bam, unmapped.bam, *.bed
	                       <b>am files: accepted_hits.bam and unmapped.bam
	                       <h>its: accepted_hits.bam only.
	                       In all cases files will go to a ./tophat/<runid>/<rpt> directory (-d is ignored).
	   -h               Show this message.
	
	EOF
}
#	   -f <format/file> Download format: <b>am, <c>ram, <f>astq (.fq.gz), <t>ophat, <o>ther; default: c.
# 	   -o <options>     Only needed if -f = o: <m>d5. Not implemented yet:  <i>ndex (bai/crai), <p>hiX, <s>tats, seqch<k>sum

#die() { echo >&2 -e "\nERROR: $@\n"; exit 1; }
#printcmd_and_run() { echo "COMMAND: $@"; "$@"; code=$?; [ $code -ne 0 ] && die "command [$*] failed with error code $code"; }

while getopts ":hi:d:f:nr:t:o:k:" OPTION; do
    case $OPTION in
        h)
            usage; exit 1;;
        i)
            TARGET=$OPTARG;;
        d)
            SDIR=$OPTARG;;
        f)
            FORMAT=$OPTARG;;
        n)
            CONVERT=0;;
        r)
            REFERENCE=$OPTARG;;
        t)
            TFOPTION=$OPTARG;;
        o)
            OTHERFILE=$OPTARG;;
        k)
            KEYTABPATH=$OPTARG;;
        \?)
            echo "Invalid option: -$OPTARG" >&2; exit 1;;
        :)
            echo "Option -$OPTARG requires an argument." >&2; exit 1;;
    esac
done


formatregex='^[bcfto]$'
tfoptionregex='^[abh]$'
otherfilergx='^[mipsk]$'
defformat=c
deftfoption=h
wrongformat=0
wrongtfoption=0
wrongofoption=0
cram2bam=1
defaultgzlevel=9

[ -z "$FORMAT" ] && FORMAT=$defformat
[ -n "$FORMAT" ] && [[ $FORMAT =~ $formatregex ]] || wrongformat=1
[ -z "$SDIR" ] && SDIR=./irods/$RUN
[ -z "$CONVERT" ] && CONVERT=$cram2bam
[ -z "$TFOPTION" ] && TFOPTION=$deftfoption
[ -n "$TFOPTION" ] && [[ $TFOPTION =~ $tfoptionregex ]] || wrongtfoption=1
#[ -z "$OTHERFILE" ] && OTHERFILE=""
#[ -n "$OTHERFILE" ] && [[ $OTHERFILE =~ $otherfilergx ]] || wrongofoption=1

[ "$wrongformat" -eq "1" ] && echo "-f: Wrong output format: $FORMAT" && exit 2
[ "$wrongtfoption" -eq "1" ] && echo "-t: Wrong option for -f t: $TFOPTION" && exit 2
[ "$wrongtfoption" -eq "1" ] && echo "-o: Wrong file type option for -o: $OTHERFILE" && exit 2



if [ -z "$KEYTABPATH" ]; then
    usage
    exit 1
fi

USER=`whoami`
KEYTABFILE=${USER}.keytab

if [ -d "$KEYTABPATH" ] && [ -e "${KEYTABPATH}/${KEYTABFILE}" ]; then
    echo "COMMAND: /usr/bin/kinit  ${USER}\@INTERNAL.SANGER.AC.UK -k -t ${KEYTABPATH}/${KEYTABFILE}"
    /usr/bin/kinit  ${USER}\@INTERNAL.SANGER.AC.UK -k -t  "${KEYTABPATH}/${KEYTABFILE}"
else
    echo "-k: cannot access ${KEYTABPATH}/${KEYTABFILE}: no such file or directory"
    exit 1
    # example: /nfs/gapi/data/keytabs/${USER}.keytab
fi



if [ -z "$TARGET" ]; then

    usage
    exit 1

elif [ -e "$TARGET" ]; then

    RUN=`head -n ${LSB_JOBINDEX} ${TARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $2}'`
    POS=`head -n ${LSB_JOBINDEX} ${TARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
    TAG=`head -n ${LSB_JOBINDEX} ${TARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`
    SDIR=`head -n ${LSB_JOBINDEX} ${TARGET} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $7}'`
        
else

    regex='^([[:digit:]]*)_([[:digit:]]*)(_([[:alpha:]]*))?(#([[:digit:]]*))?(_([[:alpha:]]*))?'

    if [[ "$TARGET" =~ $regex ]]; then
        RUN="${BASH_REMATCH[1]}";
        POS="${BASH_REMATCH[2]}";
        TAG="${BASH_REMATCH[6]}";
        SPL="${BASH_REMATCH[8]}";
    else
        echo "The run id '$TARGET' doesn't look right. Example: 14940_3#34_human"
        exit 1
    fi
    
fi

#echo run: $RUN  pos: $POS  tag: $TAG  spl: $SPL

XAMID=$RUN\_$POS

[ -n "$TAG" ] && XAMID=$RUN\_$POS\#$TAG
[ -n "$SPL" ] && XAMID=$RUN\_$POS\#$TAG\_$SPL



LSBAM="$(ils /seq/${RUN}/${XAMID}.bam 2>&1)"
NOBAM=$?

if [ "$NOBAM" -eq "4" ]; then
    iformat=cram
elif [ "$NOBAM" -eq "0" ]; then
    iformat=bam
else
    echo "Error trying to get bam/cram file for ${XAMID}: $LSBAM"
    exit 1
fi



RET_CODE=0

mkdir -p $SDIR

case "$FORMAT" in
    b)
        if [ "$NOBAM" -eq "4" ]; then
            echo "No bam format file present in iRODS for ${XAMID}: $LSBAM"
            if [ "$CONVERT" -eq "0" ]; then
                echo "$0 run with -n flag, cram->bam operation will not be done."
                exit 1
            else
                echo "Converting cram->bam ..."
                echo "COMMAND: /software/irods/icommands/bin/iget -K /seq/${RUN}/${XAMID}.cram - | scramble -I cram -O bam - ${SDIR}/${XAMID}.bam"
                /software/irods/icommands/bin/iget -K /seq/${RUN}/${XAMID}.cram - | scramble -I cram -O bam - ${SDIR}/${XAMID}.bam
                RET_CODE=$?
            fi
        elif [ "$NOBAM" -eq "0" ]; then
            echo "COMMAND: /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.bam ${SDIR}/${XAMID}.bam"
            /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.bam ${SDIR}/${XAMID}.bam
            RET_CODE=$?
        fi
        ;;
    c)
        echo "COMMAND: /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.cram ${SDIR}/${XAMID}.cram"
        /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.cram ${SDIR}/${XAMID}.cram
        RET_CODE=$?
        ;;
    f)
        if [ "$NOBAM" -eq "4" ]; then
            REFERENCEOPTION=""
            if [ -z "$REFERENCE" ]; then
                echo "-f: A reference is needed: -r $REFERENCE ."
                exit 1
            else
                REFERENCEOPTION="reference=$REFERENCE"
            fi
        fi
        SDIR=./bamfq/$RUN
        FQ1=$SDIR/${XAMID}_1.fq.gz
        FQ2=$SDIR/${XAMID}_2.fq.gz
        mkdir -p $SDIR
        echo "COMMAND: /software/irods/icommands/bin/iget /seq/${RUN}/${XAMID}.$iformat - | /software/npg/bin/bamsort SO=queryname level=0 inputformat=$iformat $REFERENCEOPTION outputformat=bam index=0 tmpfile=./bamfq/ | /software/npg/bin/bamtofastq exclude=QCFAIL gz=1 level=$defaultgzlevel F=$FQ1 F2=$FQ2"
        /software/irods/icommands/bin/iget /seq/${RUN}/${XAMID}.$iformat - | \
            /software/npg/bin/bamsort SO=queryname level=0 inputformat=$iformat $REFERENCEOPTION outputformat=bam index=0 tmpfile=./bamfq/ | \
            /software/npg/bin/bamtofastq exclude=QCFAIL gz=1 level=$defaultgzlevel F=$FQ1 F2=$FQ2
        RET_CODE=$?
        ;;
    t)
        SDIR=./tophat/$RUN/$XAMID
        mkdir -p $SDIR
        ret_code_junc=0
        ret_code_dels=0
        ret_code_hits=0
        ret_code_unmp=0
        echo "COMMAND: /software/npg/bin/samtools1 view -bh -F 0x4 -o ${SDIR}/accepted_hits.bam irods:/seq/${RUN}/${XAMID}.$iformat"
        ACCEPTEDHITSFILE="$(/software/npg/bin/samtools1 view -bh -F 0x4 -o ${SDIR}/accepted_hits.bam irods:/seq/${RUN}/${XAMID}.$iformat 2>&1)"
        ret_code_hits=$?
        if [ "$TFOPTION" = "a" ] || [ "$TFOPTION" = "b" ] ; then
            echo "COMMAND: /software/npg/bin/samtools1 view -bh -f 0x4 -o ${SDIR}/unmapped.bam irods:/seq/${RUN}/${XAMID}.$iformat"
            UNMAPPEDFILE="$(/software/npg/bin/samtools1 view -bh -f 0x4 -o ${SDIR}/unmapped.bam irods:/seq/${RUN}/${XAMID}.$iformat 2>&1)"
            ret_code_unmp=$?
        fi
        if [ "$TFOPTION" = "a" ]; then
            echo "COMMAND: /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.junctions.bed  ${SDIR}/junctions.bed"
            JUNCTIONSFILE="$(/software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.junctions.bed  ${SDIR}/junctions.bed 2>&1)"
            ret_code_junc=$?
            echo "COMMAND: /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.deletions.bed  ${SDIR}/deletions.bed"
            DELETIONSFILE="$(/software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.deletions.bed  ${SDIR}/deletions.bed 2>&1)"
            ret_code_dels=$?
        fi
        let "RET_CODE = ret_code_junc + ret_code_dels + ret_code_hits + ret_code_unmp"
        if [ "$RET_CODE" -ne "0" ]; then
            echo "-f t: One or more files failed to download"
            [ "$ret_code_junc" -ne "0" ] && echo "iget: $JUNCTIONSFILE"
            [ "$ret_code_dels" -ne "0" ] && echo "iget: $DELETIONSFILE"
            [ "$ret_code_hits" -ne "0" ] && echo "iget: $ACCEPTEDHITSFILE"
            [ "$ret_code_unmp" -ne "0" ] && echo "iget: $UNMAPPEDFILE"
        fi
        ;;
#    o)
#        case "$OTHERFILE" in
#            m)
#                echo "COMMAND: /software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.deletions.bed  ${SDIR}/deletions.bed"
#                DELETIONSFILE="$(/software/irods/icommands/bin/iget -KPvf /seq/${RUN}/${XAMID}.deletions.bed  ${SDIR}/deletions.bed 2>&1)"
#                ret_code_dels=$?
#                ;;
#            *)
#                echo "-o: Wrong file type option for -o: $OTHERFILE" && exit 2 ;;
#        esac
        
esac

# NOTE: Exit code 3 seems to be the one for failed download.

exit $RET_CODE
