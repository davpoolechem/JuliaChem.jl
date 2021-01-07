#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export JULIA_NUM_THREADS=56
export OPENBLAS_NUM_THREADS=1
#
mpirun -np 1 julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes /home/davpoolechem/shared/projects/Julia/JuliaChem.jl/example_scripts/full-rhf.jl example_inputs/ARPA-E/6-311/fig1a_56thread.json
