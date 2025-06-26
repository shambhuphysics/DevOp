#!/bin/bash -l
#$ -N ZnO_001
#$ -pe mpi 40                # 40 physical cores (1 node) using 2 nodes
#$ -l h_rt=40:30:00
#$ -l mem=4.8G               # 192 GB/node (40 x 4.8G) = total 2 nodes
#$ -cwd                      # Use the current working directory
#$ -o $JOB_NAME.$JOB_ID.out  # Output file includes job name and ID
#$ -e $JOB_NAME.$JOB_ID.err  # Error file includes job name and ID
#$ -M ucfbsbh@ucl.ac.uk      # Email notifications to your UCL address
#$ -m bea                    # Email at (b)egin, (e)nd, (a)bort


module unload -f compilers mpi
module load compilers/intel/2019/update5
module load mpi/intel/2019/update5/intel
module load vasp/6.3.0-24Jan2022/intel-2019-update5

./elec1.1.sh
