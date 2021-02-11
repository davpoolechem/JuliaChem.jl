#!/bin/bash
#
#SBATCH --job-name=juliachem
#SBATCH --output=juliachem.log
#SBATCH --error=juliachem.err
#
#SBATCH --partition=compute
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#
export JULIA_NUM_THREADS=1
#

srun ./builddir/minimal-rhf-benchmark example_inputs/S22/01.inp 
