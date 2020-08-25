#!/bin/bash
#
#SBATCH --job-name=example_inputs/ARPA-E/6-311/fig1i
#SBATCH --output=example_inputs/ARPA-E/6-311/fig1i.log
#SBATCH --error=example_inputs/ARPA-E/6-311/fig1i.err
#
#SBATCH --partition=compute
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#
export JULIA_NUM_THREADS=6
export OMP_NUM_THREADS=6
#
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes example_scripts/minimal-rhf.jl example_inputs/ARPA-E/6-311/fig1i.json
