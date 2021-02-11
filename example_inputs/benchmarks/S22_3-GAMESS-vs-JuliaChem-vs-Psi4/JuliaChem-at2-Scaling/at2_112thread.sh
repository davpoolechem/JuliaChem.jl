#!/bin/bash
#
#COBALT -n 3 
#COBALT -t 30
#COBALT -q skylake_8180
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export JULIA_NUM_THREADS=112
export OPENBLAS_NUM_THREADS=1
#
${MPI_HOME}/bin/orterun --mca btl vader,self,tcp -np 1 --map-by ppr:1:node --bind-to none --report-bindings julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes /home/davpoolechem/shared/projects/Julia/JuliaChem.jl/example_scripts/minimal-rhf-benchmark.jl adenine_thymine_2.json
