#!/bin/sh
#$ -N H2542019
#$ -j y
#$ -cwd
#$ -pe smp 56
####$ -pe 56cpn 56
####$ -l mf=16G
#$ -q IFC

source /Users/njadidoleslam/virtenvs/smap_assim/bin/activate
python multithread_test.py