#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export JULIA_NUM_THREADS=56
export OPENBLAS_NUM_THREADS=1
#
export RANKS_PER_NODE=1
#
${MPI_HOME}/bin/orterun --mca btl vader,self,tcp -np ${RANKS_PER_NODE} -map-by ppr:${RANKS_PER_NODE}:node --bind-to none --report-bindings julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes /home/davpoolechem/shared/projects/Julia/JuliaChem.jl/example_scripts/minimal-rhf-benchmark.jl example_inputs/ARPA-E/6-311/fig1e.json
