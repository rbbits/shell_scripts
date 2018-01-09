#!/bin/bash

. /software/gapi/Modules/default/init/bash

module load samtools

sample=`head -n ${LSB_JOBINDEX} ${1} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $3}'`
study=`head -n ${LSB_JOBINDEX} ${1} | tail -1 | awk 'BEGIN { FS = "\t" } ; {print $4}'`

if [ $2 = "star" ]; then

    stardir=star/$study/$sample
    bamfile=$study.$sample.bam

    samtools index $stardir/$bamfile

elif [ $2 = "tophat" ]; then

    tophatdir=tophat/$study/$sample
    bamfile=accepted_hits.bam

    samtools index $tophatdir/$bamfile

fi
