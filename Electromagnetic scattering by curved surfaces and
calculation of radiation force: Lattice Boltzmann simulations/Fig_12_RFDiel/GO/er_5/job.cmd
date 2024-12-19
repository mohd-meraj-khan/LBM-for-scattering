#!/bin/bash
#PBS -e errorfile.err
#PBS -o logfile.log
#PBS -q anubhab_q
#PBS -l select=1:ncpus=1
tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp -R $PBS_O_WORKDIR/* .
module load anaconda3_2019
module load ffmpeg4.3
./bash.sh
mv ../job$tpdir $PBS_O_WORKDIR/.

