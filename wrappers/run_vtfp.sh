#!/bin/bash

. /software/npg/etc/profile.npg
. /software/sanger-samtools-refpath/etc/profile.sanger-samtools-refpath


usage(){
    cat <<-EOF
	This script generates VTFP templates using the information provided by a targets file.
	[M]ethods available:  bwa_mem | bwa_aln | tophat2 | bam2cram | y_split | hs_split | salmon
		
	Usage: 

	$0 -M <METHOD> [options] <targets_file.txt
	
	Options:
	   -c <number>      Do not use WTSI composite id (run_position[#tag]), use insted the
	                    contents of this column in the targets file (must be unique).
	   -f <format>      Input file format: cram | bam. Default: auto-detect.
	   -h               Show usage message.
	   -i <directory>   Input directory. Default: <current|working directory>/irods/<run>.
	   -r <directory>   Absolute path to repository for reference genome/transcriptome
	   -w <directory>   Absolute path to working directory. Default: $PWD
	   -x <extra args>  Extra arguments passed to vtfp in a quoted string (-k/-v pairs)
	EOF

}

# [2017-11-25] Not supported yet:
#	   -o <directory>   Output directory. Folder structure will be created if it does not exist,
#	                    relative to ./ or a working directory provided by -w unless an absolute path is provided.
#	                    Default: <. or working directory>/output/<method>/<run>/<bamid>

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

while getopts ":c:i:f:hm:M:r:w:x:" OPTION; do
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
        o)
            # [2017-11-25] not supported yet
            OUTPUTDIR=$OPTARG
            [ -d $OUTPUTDIR ] || printf -- "[WARNING] Cannot access ${INDIR}: No such directory and it will be created\n";;
        M)
            METHOD=$OPTARG
            METHODREGEX="^bam2cram|star|tophat2|bwa\_aln|bwa\_mem|hs\_split|y\_split|salmon$"
            [ -z "$METHOD" ] && exitmessage "[ERROR] -M: a method is required: try '$0 -h' for more information" 1
            [ -n "$METHOD" ] && [[ ! $METHOD =~ $METHODREGEX ]] && exitmessage "[ERROR] -M: invalid method $METHOD: try '$0 -h' for more information" 1;;
        r)
            REPOSITORY=$OPTARG
            [[ ! -d $REPOSITORY ]] && exitmessage "[ERROR] -w: Cannot access ${REPOSITORY}: no such directory" 2
            [[ ! $REPOSITORY = /* ]] && exitmessage "[ERROR] -w: Not an absolute path" 1;;
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



if env | grep -q ^P4_PATH=
then
    # use whatever version of p4 you defined like this:
    # export P4_PATH=/nfs/users/nfs_r/rb11/dev/perl/wtsi-npg_p4
    printf -- "[INFO] Using P4_PATH=$P4_PATH\n"
    BINARY="${P4_PATH}/bin/vtfp.pl"
else
    # for production stuff e.g. remapping, use production versions
    # for development/experimental/non-deployed stuff e.g. bam2cram use own repo (see path above)
    if [ $METHOD = "bam2cram" ]; then
        exitmessage "[ERROR] Not P4_PATH env variable: try 'export P4_PATH=/path/to/p4'" 1
    else
        BINARY=$(readlink -f `which vtfp.pl`)
        export P4_PATH=${BINARY%/bin*}
    fi
fi

CFGDATADIR="${P4_PATH}/data/vtlib"
VTFPEXECUTABLE=$BINARY

[ -e "json/vtfp_commands_${METHOD}.log" ] && rm -f json/vtfp_commands_${METHOD}.log

COUNTTOTAL=0
COUNTOK=0
COUNTFAIL=0

while read line; do

    # Read info from targets file
    if [ ! -z $TARGETSCOLUMN ]; then
        BAMID=`echo "$line" | awk -v column=$TARGETSCOLUMN -F'\t' '{print $column}'`
    else
        RUN=`echo "$line" | awk -F'\t' '{print $2}'`
        POSITION=`echo "$line" | awk -F'\t' '{print $3}'`
        TAG=`echo "$line" | awk -F'\t' '{print $4}'`
        # deal with non-multiplexed lanes
        if [ -z $TAG ]; then
            BAMID="${RUN}_${POSITION}"
        else
            BAMID="${RUN}_${POSITION}#${TAG}"
        fi
    fi
    ALIGNMENTSINBAM=`echo "$line" | awk -F'\t' '{print $5}'`
    ALIGNREFGENOME=`echo "$line" | awk -F'\t' '{print $6}'`
    REFDICTNAME=`echo "$line" | awk -F'\t' '{print $7}'`
    REFNAMEFASTA=`echo "$line" | awk -F'\t' '{print $8}'`
    TRANSCRIPTOME=`echo "$line" | awk -F'\t' '{print $9}'`
    TRANSCRIPTANNO=`echo "$line" | awk -F'\t' '{print $10}'`
    LIBRARYTYPE=`echo "$line" | awk -F'\t' '{print $11}'`
    LIBRARYLAYOUT=`echo "$line" | awk -F'\t' '{print $12}'`
    REFTRANSCRIPTFASTA=`echo "$line" | awk -F'\t' '{print $13}'`

    OBAMDIR=${OUTPUTDIR:-"output/$METHOD/$RUN/$BAMID"}
    IBAMDIR=${INPUTDIR:-"irods/$RUN"}
    WORKINGDIR=${CWD:-"$PWD"}
    INDATADIR="${WORKINGDIR}/${IBAMDIR}"
    OUTDATADIR="${WORKINGDIR}/${OBAMDIR}"
    OUTSTAGINGDIR="{WORKINGDIR}/staging/${METHOD}/${RUN}/${BAMID}"
    OUTJSONDIR="${WORKINGDIR}/json"
    SRCINPUT="${INDATADIR}/${BAMID}"
    REPOSDIR=${REPOSITORY-"/lustre/scratch117/core/sciops_repository"}
    PHIXDICTNAME="PhiX/default/all/picard/phix_unsnipped_short_no_N.fa.dict"
    PHIXREFNAME="PhiX/Sanger-SNPs/all/fasta/phix_unsnipped_short_no_N.fa"
    ALIGNMENTFILTERJAR="/software/solexa/pkg/illumina2bam/1.19/AlignmentFilter.jar"

    ALIGNMENTMETHOD=$METHOD
    SALMON_TRANSCRIPTOME="$(dirname $TRANSCRIPTOME)"
    SALMON_TRANSCRIPTOME="$(dirname $SALMON_TRANSCRIPTOME)"
    SALMON_TRANSCRIPTOME+="/salmon"

    case $METHOD in
        tophat2)
            if [[ ! $TRANSCRIPTOME = NoTranscriptome && ! $TRANSCRIPTANNO = NoTranscriptome ]]; then
                TOPHAT2_ARGS="-keys transcriptome_val "
                # full or relative path should be OK
                [[ ! $TRANSCRIPTOME = /* ]] && TOPHAT2_ARGS+="-vals ${REPOSDIR}/transcriptomes/${TRANSCRIPTOME} " || TOPHAT2_ARGS+="-vals ${TRANSCRIPTOME} "
                TOPHAT2_ARGS+="-keys salmon_transcriptome_val "
                [[ ! $TRANSCRIPTOME = /* ]] && TOPHAT2_ARGS+="-vals ${REPOSDIR}/transcriptomes/${SALMON_TRANSCRIPTOME} " || TOPHAT2_ARGS+="-vals ${SALMON_TRANSCRIPTOME} "
                #TOPHAT2_ARGS+="-keys annotation_val " #only when there's no transcriptome index for tophat
                #[[ ! $TRANSCRIPTOME = /* ]] && TOPHAT2_ARGS+="-vals ${REPOSDIR}/transcriptomes/${TRANSCRIPTANNO} " || TOPHAT2_ARGS+="-vals ${TRANSCRIPTANNO} "
            else
                exitmessage "[ERROR] NoTranscriptome for Tophat2 alignment" 1
            fi
            [[ $LIBRARYTYPE =~ dUTP ]] && TOPHAT2_ARGS+="-keys library_type -vals fr-firststrand" || TOPHAT2_ARGS+="-keys library_type -vals fr-unstranded"
            ;;
        star)
            STAR_ARGS="-keys annotation_val -vals ${REPOSDIR}/transcriptomes/${TRANSCRIPTANNO} "
            STAR_ARGS+="-keys sjdb_overhang_val -vals 74 "
            STAR_ARGS+="-keys star_executable -vals star "
            #STAR_ARGS+="-keys reference_transcriptome_fasta -vals ${REPOSDIR}/transcriptomes/${REFTRANSCRIPTFASTA} "
            #STAR_ARGS+="-keys scramble_reference_fasta -vals ${REPOSDIR}/references/${REFNAMEFASTA} "
            if [[ ! $ALIGNREFGENOME =~ .*/star$ ]]; then
                ALIGNREFGENOME="$(dirname $ALIGNREFGENOME)"
                ALIGNREFGENOME="$(dirname $ALIGNREFGENOME)"
                ALIGNREFGENOME+="/star"
            fi
            # this is only for Salmon!:
            if [[ ! $TRANSCRIPTOME = NoTranscriptome && ! $TRANSCRIPTANNO = NoTranscriptome ]]; then
                STAR_ARGS+=" -keys salmon_transcriptome_val "
                # full or relative path should be OK
                [[ ! $TRANSCRIPTOME = /* ]] && STAR_ARGS+="-vals ${REPOSDIR}/transcriptomes/${SALMON_TRANSCRIPTOME}" || STAR_ARGS+="-vals ${SALMON_TRANSCRIPTOME}"
            else
                exitmessage "[ERROR] NoTranscriptome for Salmon quantification" 1
            fi            
            ;;
        salmon)
            if [[ ! $TRANSCRIPTOME = NoTranscriptome && ! $TRANSCRIPTANNO = NoTranscriptome ]]; then
                SALMON_ARGS="-keys transcriptome_val "
                # full or relative path should be OK
                [[ ! $TRANSCRIPTOME = /* ]] && SALMON_ARGS+="-vals ${REPOSDIR}/transcriptomes/${SALMON_TRANSCRIPTOME}" || SALMON_ARGS+="-vals ${SALMON_TRANSCRIPTOME} "
                SALMON_ARGS+="-keys annotation_val "
                [[ ! $TRANSCRIPTOME = /* ]] && SALMON_ARGS+="-vals ${REPOSDIR}/transcriptomes/${TRANSCRIPTANNO} " || SALMON_ARGS+="-vals ${TRANSCRIPTANNO}"
            else
                exitmessage "[ERROR] NoTranscriptome for Salmon quantification" 1
            fi
            ;;
        hs_split)
            ALIGNMENTMETHOD=bwa_mem
            ALIGNMENTTEMPLATE="realignment_wtsi_stage2_humansplit_template.json"
            HS_SPLIT_ARGS="-keys reference_dict_hs -vals ${REPOSDIR}/references/Homo_sapiens/1000Genomes/all/picard/human_g1k_v37.fasta.dict "
            HS_SPLIT_ARGS+="-keys hs_reference_genome_fasta -vals ${REPOSDIR}/references/Homo_sapiens/1000Genomes/all/fasta/human_g1k_v37.fasta "
            HS_SPLIT_ARGS+="-keys hs_alignment_reference_genome -vals ${REPOSDIR}/references/Homo_sapiens/1000Genomes/all/bwa0_6/human_g1k_v37.fasta "
            HS_SPLIT_ARGS+="-keys alignment_filter_jar -vals /software/solexa/pkg/illumina2bam/1.17/AlignmentFilter.jar "
            HS_SPLIT_ARGS+="-keys alignment_hs_method -vals bwa_aln "
            ;;
        y_split)
            ALIGNMENTMETHOD=bwa_mem
            # 20170609: in p4 0.18.6 bambi is used instead of illumina2bam. Is this by default?
            Y_SPLIT_ARGS="-keys split_bam_by_chromosomes_jar -vals /software/solexa/pkg/illumina2bam/1.17/SplitBamByChromosomes.jar "
            Y_SPLIT_ARGS+="-keys final_output_prep_target_name -vals split_by_chromosome "
            Y_SPLIT_ARGS+="-keys split_indicator -vals _yhuman "
            Y_SPLIT_ARGS+="-keys split_bam_by_chromosome_flags -vals S=Y "
            Y_SPLIT_ARGS+="-keys split_bam_by_chromosome_flags -vals V=true "
            Y_SPLIT_ARGS+="-keys s2b_mt_val -vals 7"
            ;;
    esac

    if [ -e "${SRCINPUT}.cram" ]; then
        SRCINPUTEXT="cram"
    elif [ -e "${SRCINPUT}.bam" ]; then
        SRCINPUTEXT="bam"
    elif [ -e "${SRCINPUT}.sam" ]; then
        SRCINPUTEXT="sam"
    else
        printf -- "[ERROR] ${BAMID}: No bam or cram or sam file was found in ${INDATADIR}/\n"
        exitmessage "[INFO] Use option -i to specify an input directory relative to ${WORKINGDIR}" 1
    fi

    ######################
    # Generate JSON files
    ######################
    printf -- "[INFO] Generating [${METHOD}] json file for [${ALIGNMENTSINBAM}] [${BAMID}] in [./json/${BAMID}_${METHOD}.json] with source format [${SRCINPUTEXT}]\n"
    [ ! -z "$EXTRAKEYVALS" ] && printf -- "[INFO] Using extra arguments [ ${EXTRAKEYVALS} ]\n"

    # start building vtfp command
    VTFP_CMD="${VTFPEXECUTABLE} -l json/vtfp_${BAMID}_${METHOD}.log "
    VTFP_CMD+="-ve 3 "
    VTFP_CMD+="-o ${OUTJSONDIR}/${BAMID}_${METHOD}.json "
    VTFP_CMD+="-keys rpt -vals ${BAMID} "
    VTFP_CMD+="-keys src_input_ext -vals ${SRCINPUTEXT} "
    VTFP_CMD+="-keys src_input_format -vals ${SRCINPUTEXT} "
    VTFP_CMD+="-keys outdatadir -vals ${OUTDATADIR} "
    VTFP_CMD+="-keys indatadir -vals ${INDATADIR} "
    VTFP_CMD+="-keys cfgdatadir -vals ${CFGDATADIR} "
    
    case $METHOD in
        bwa_aln|bwa_mem|tophat2|star|hs_split|y_split)
            # bwa args apply to all of these methods
            BWA_ARGS="-keys bwa_executable -vals bwa0_6 "
            [ "$LIBRARYLAYOUT" = "SINGLE" ] && BWA_ARGS+="-nullkeys bwa_mem_p_flag"

            # templates used for these methods
            [ "$ALIGNMENTSINBAM" = "aligned" ] && ALIGNMENTTEMPLATE="realignment_wtsi_template.json" || ALIGNMENTTEMPLATE="alignment_wtsi_stage2_template.json"

            # for realignment at least this prunning has to be there, more can be included in the -x option
            [[ $EXTRAKEYVALS =~ prune ]] && PRUNE_NODES_ARGS="" || PRUNE_NODES_ARGS="-prune_nodes fop.*samtools_stats_F0.*00_bait.*-"

            [[ $METHOD =~ tophat2|star ]] && QUANT_METHOD_ARGS="-keys quant_method -vals salmon"
            
            VTFP_CMD+="-keys samtools_executable -vals samtools1 "
            VTFP_CMD+="-keys alignment_method -vals ${ALIGNMENTMETHOD} "
            VTFP_CMD+="-keys af_metrics -vals ${BAMID}.bam_alignment_filter_metrics.json "
            VTFP_CMD+="-keys reference_dict -vals ${REPOSDIR}/references/${REFDICTNAME} "
            VTFP_CMD+="-keys reference_genome_fasta -vals ${REPOSDIR}/references/${REFNAMEFASTA} "
            VTFP_CMD+="-keys alignment_reference_genome -vals ${REPOSDIR}/references/${ALIGNREFGENOME} "
            VTFP_CMD+="-keys phix_reference_genome_fasta -vals ${REPOSDIR}/references/${PHIXREFNAME} "
            VTFP_CMD+="-keys alignment_filter_jar -vals ${ALIGNMENTFILTERJAR} "
            VTFP_CMD+="-keys aligner_numthreads -vals 16 "
            VTFP_CMD+="-keys br_numthreads_val -vals 7 "
            VTFP_CMD+="-keys b2c_mt_val -vals 7 "
            VTFP_CMD+="-keys s2b_mt_val -vals 7 "
            VTFP_CMD+="${TOPHAT2_ARGS} ${STAR_ARGS} ${QUANT_METHOD_ARGS} ${HS_SPLIT_ARGS} ${Y_SPLIT_ARGS} ${BWA_ARGS} "

            EXPORT_PV_JSON="-export_param_vals ${OUTJSONDIR}/${BAMID}_p4_${ALIGNMENTMETHOD}_realignment_pv_out.json "
            ;;
        salmon)
            # exclusive template for this method
            ALIGNMENTTEMPLATE="bam_to_salmon.json"

            VTFP_CMD+="-keys b2c_mt_val -vals 4 "
            VTFP_CMD+="${SALMON_ARGS} "

            EXPORT_PV_JSON="-export_param_vals ${OUTJSONDIR}/${BAMID}_p4_${METHOD}_pv_out.json "
            ;;
        bam2cram)
            exitmessage "[ERROR] -M: option not supported: $METHOD" 1
            # 20170308 rb11: bam2cram command needs to be updated to run on the latest version of P4
            # TODO: new command
            # prioritise use value from targets file
            ALIGNMENTSINBAM=${ALIGNMENTSINBAM:-$ALIGNMENTSINBAM_OPT}
            ALIGNMENTSINBAMREGEX="^unaligned|aligned$"
            [ -n "$ALIGNMENTSINBAM" ] && [[ ! $ALIGNMENTSINBAM =~ $ALIGNMENTSINBAMREGEX ]] && exitmessage "[ERROR] -m: invalid sub-method ${ALIGNMENTSINBAM}: possible values are: aligned | unaligned" 1
            JSONFILE=${METHOD}_${ALIGNMENTSINBAM}.json
            ;;
        *)
           exitmessage "[ERROR] -M: option not supported: $METHOD" 1 
    esac

    # final touches to the vtfp command
    VTFP_CMD+="${EXTRAKEYVALS} ${PRUNE_NODES_ARGS} ${EXPORT_PV_JSON}"
    VTFP_CMD+="${CFGDATADIR}/${ALIGNMENTTEMPLATE} "
    
    echo $VTFP_CMD >> json/vtfp_commands_${METHOD}.log
        
    RET_CODE=0
    VTFP="$($VTFP_CMD 2>&1)"
    RET_CODE=$?
    
    let COUNTTOTAL+=1
    COUNTBAMID+=("$BAMID")
 
    if [ "$RET_CODE" -eq 0 ]; then
        let COUNTOK+=1
    else
        printf -- "[INFO] VTFP command for [ ${BAMID} ] exited with exit code ${RET_CODE}\n"
        printf -- "[ERROR] $VTFP\n"
        let COUNTFAIL+=1
    fi
     
done < /dev/stdin

COUNTUNIQ=`tr ' ' '\n' <<< "${COUNTBAMID[@]}" | sort -u | wc -l`

if [ "$COUNTTOTAL" -gt "$COUNTUNIQ" ]; then
    MESSAGEDUPS=" or are duplicated"
    COUNTOK=$COUNTUNIQ
    let COUNTFAIL="$COUNTFAIL + ($COUNTTOTAL - $COUNTUNIQ)"
fi

if [ "$COUNTTOTAL" -eq "$COUNTOK" ]; then
    MESSAGE="Done. [ ${COUNTOK} ] command(s) executed successfully"
    RET_CODE=0
else
    MESSAGE="[ ${COUNTFAIL} ] command(s) exited with errors ${MESSAGEDUPS}"
    RET_CODE=1
fi

exitmessage "[INFO] $MESSAGE" $RET_CODE
