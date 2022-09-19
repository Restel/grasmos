#!/bin/bash
# 2022-06-10

# To run:
# ./color_models.sh
# OR
# bash color_models.sh

PARAMS="./parameters/"
OUTPUT="./output/"

for f in $(ls -1 $PARAMS*.txt); do
  LOGFILE=$OUTPUT$(basename $f .log)
  sbatch -o $LOGFILE mle_run.sh $f $LOGFILE  # RC call
  #/bin/bash mle_run_RC.sh $f $LOGFILE # laptop call
done
