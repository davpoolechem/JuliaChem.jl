#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export OMP_NUM_THREADS=111
export MKL_NUM_THREADS=1
#
julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes 6-311++G_2d_2p/indole_benzene_2_t_shaped.jl
