#!/bin/bash
#
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1
#
cd /home/davpoolechem/shared/Sandbox/Julia/JuliaChem.jl
mpirun -np 16 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes example_scripts/minimal-rhf-benchmark.jl example_inputs/Scaling-Memory/w50-16rank.json
