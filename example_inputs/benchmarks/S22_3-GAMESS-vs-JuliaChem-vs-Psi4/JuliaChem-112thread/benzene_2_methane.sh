#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export JULIA_NUM_THREADS=112
export OPENBLAS_NUM_THREADS=1
#
/home/davpoolechem/programs/install/openmpi/bin/orterun --mca btl vader,self,tcp -np 1 -map-by ppr:1:node --bind-to none --report-bindings julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes /home/davpoolechem/shared/projects/Julia/JuliaChem.jl/example_scripts/minimal-rhf-benchmark.jl example_inputs/S22_3/6-311++G_2d_2p/benzene_2_methane.json
