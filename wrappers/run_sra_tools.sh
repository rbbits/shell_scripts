#!/bin/bash

#set -x

if [ -f /software/gapi/Modules/default/init/bash ]; then
    . /software/hgi/modules/default/init/bash
    module load hgi/ncbi-sra-tools/2.8.0
fi

TARGETSFILE=$1
TARGETSCOLUMN=13
LINE=$LSB_JOBINDEX

ACCESSION=`head -n ${LINE} ${TARGETSFILE} | tail -1 | awk -v COLUMN=${TARGETSCOLUMN} 'BEGIN { FS = "\t" } ; {print $COLUMN}'`

INPUTDIR="./input/sra/${ACCESSION}"
INPUTFILE="${INPUTDIR}/${ACCESSION}.sra"
OUTPUTDIR="./input/sra_bam/${ACCESSION}"
OUTPUTFILE="${OUTPUTDIR}/${ACCESSION}.bam"
SRA_TOOL_CMD="sam-dump -r ${INPUTFILE} | samtools view -b -o ${OUTPUTFILE} -"

mkdir -p $OUTPUTDIR && echo "[INFO] Output directory ${OUTPUTDIR} was created"

echo "[INFO] Converting sra -> bam"
echo "[COMMAND] ${SRA_TOOL_CMD}"

SRA_TOOL="$(sam-dump -r -u ${INPUTFILE} | samtools view -b -o ${OUTPUTFILE} - 2>&1)"

RET_CODE=$?

if [ "$RET_CODE" -eq "0" ]; then
    echo "[INFO] Done successfully"
else
    echo "[INFO] Exited with exit code ${RET_CODE}"
    echo "[ERROR] ${SRA_TOOL}"
fi

exit $RET_CODE

