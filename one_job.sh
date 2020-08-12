#!/bin/bash

#PBS -M C.L.M.vandenHeuvel@student.tudelft.nl  	 # specify e-mail address to send message to
#PBS -m ae                               	 # mail when job is aborted (a) or terminates (e)
#PBS -l nodes=1:ppn=1 		 # number of nodes & 

ID=${PBS_JOBID:0:6}		# get number of job-idb
WORKDIR=$PBS_O_WORKDIR		# get the directory you're working in

exec 1>$WORKDIR/${ID}.out       # create output file for systemt
exec 2>$WORKDIR/${ID}.err       # create error file

# ececute code

cd $WORKDIR

module load matlab
matlab -nosplash -nodisplay -noFigureWindows -nodesktop -r "$SCRIPT_NAME $x"

