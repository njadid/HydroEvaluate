#!/bin/sh
#$ -N Metrics_SM
#$ -j y
#$ -cwd
#$ -pe smp 56
####$ -pe 56cpn 56
####$ -l mf=16G
#$ -q IFC

/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

source /Users/njadidoleslam/virtenvs/smap_assim/bin/activate
python SM_metrics_MP.py