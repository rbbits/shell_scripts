#!/bin/bash
#old way: runfolder=`perl -le 'use srpipe::runfolder; print srpipe::runfolder->new(id_run=>$ARGV[0])->runfolder_path' $1`
runfolder=`perl -le 'use npg_tracking::illumina::runfolder; print npg_tracking::illumina::runfolder->new(id_run=>$ARGV[0])->runfolder_path' $1`
if [ -d $runfolder/Data ]; then
    if [ -d $runfolder/Data/Intensities ]; then
        cd $runfolder/Data/Intensities
        ls -l | grep -E 'Bustard|BAM_basecalls'
    else
        cd $runfolder/Data
        ls -l
    fi
    cd $runfolder
    ls -l | grep Latest_Summary
else
    cd $runfolder
fi
log st $*
