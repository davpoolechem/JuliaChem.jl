#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export JULIA_NUM_THREADS=56
#
cd /home/davpoolechem/shared/projects/Julia/JuliaChem.jl
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes example_scripts/minimal-rhf-benchmark.jl example_inputs/ARPA-E/6-311/fig1a.json
