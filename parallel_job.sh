#!/bin/bash

#PBS -N v9_test			 # job name
#PBS -M E.C.Nio@student.tudelft.nl  	 # specify e-mail address to send message to
#PBS -m ae                               	 # mail when job is aborted (a) or terminates (e) 
#PBS -l nodes=5:ppn=10                           # number of nodes & 


ID=${PBS_JOBID:0:6}		# get number of job-idb
WORKDIR=$PBS_O_WORKDIR		# get the directory you're working in
cd $WORKDIR			# go to this working directory

# Create output directory
name=${PBS_JOBNAME%-*}
OUTDIR=$WORKDIR/${ID}_${name}

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR;
fi

# Redirect stdout & stderr
exec 1>$OUTDIR/${ID}.out	# create output file for systemt
exec 2>$OUTDIR/${ID}.err	# create error file

# Copy files to output directory
scp $WORKDIR/*.m $OUTDIR
scp -r $WORKDIR/Data $OUTDIR

cd $OUTDIR

module load matlab

myfunction() {
	matlab -nosplash -nodisplay -noFigureWindows -nodesktop -r "Gillespie_func $1"
}
export -f myfunction

parallel myfunction ::: {10001..20000}

# execute code



