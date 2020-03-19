#!/bin/bash
#
#SBATCH --job-name=juliachem_runtest
#SBATCH --output=juliachem_runtest.log
#SBATCH --error=juliachem_runtest.err
#
#SBATCH --partition=haswell
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#
export JULIA_NUM_THREADS=1
#
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes runtest.jl 
