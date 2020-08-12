#!/bin/bash 

echo "Hello!"

#WORKDIR=$PBS_O_WORKDIR		# get the directory you're working in
#cd $WORKDIR			# go to this working directory
#echo $PBS_O_WORKDIR

# Create output directory
#name="test_directory"
OUTDIR="./0812_1_100_sims"
JOB_NAME="0812_1_100_sims"
MATLAB_SCRIPT="Gillespie_func"
MATLAB_FILE="${MATLAB_SCRIPT}.m"

if [ ! -d $OUTDIR ] 
then
  mkdir -p $OUTDIR;
else
  echo "Fatal error: Outdir already exists."
  exit 1
fi

# Redirect stdout & stderr
#exec 1>$OUTDIR/${ID}.out	# create output file for systemt
#exec 2>$OUTDIR/${ID}.err	# create error file

# Copy files to output directory
scp one_job.sh $MATLAB_FILE $OUTDIR
scp -r Functions/* $OUTDIR

cd $OUTDIR
mkdir -p "Output"

# ececute code
for i in {0..99}
do
	NAME="${JOB_NAME}${i}"
	qsub -N $NAME -v SCRIPT_NAME=$MATLAB_SCRIPT,x=$i one_job.sh
	echo "I submitted job ${NAME}"
done
