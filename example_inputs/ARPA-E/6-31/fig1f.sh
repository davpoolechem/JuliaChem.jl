#!/bin/bash
#
#SBATCH --job-name=example_inputs/ARPA-E/6-31/fig1f
#SBATCH --output=example_inputs/ARPA-E/6-31/fig1f.log
#SBATCH --error=example_inputs/ARPA-E/6-31/fig1f.err
#
#SBATCH --partition=compute
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
#
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes example_scripts/minimal-rhf.jl example_inputs/ARPA-E/6-31/fig1f.json
