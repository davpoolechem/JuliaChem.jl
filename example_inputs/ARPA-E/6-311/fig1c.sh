#!/bin/bash
#
#SBATCH --job-name=example_inputs/ARPA-E/6-311/fig1c
#SBATCH --output=example_inputs/ARPA-E/6-311/fig1c.log
#SBATCH --error=example_inputs/ARPA-E/6-311/fig1c.err
#
#SBATCH --partition=haswell
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#
export JULIA_NUM_THREADS=8
export OMP_NUM_THREADS=8
#
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes example_scripts/minimal-rhf-benchmark.jl example_inputs/ARPA-E/6-311/fig1c.json
