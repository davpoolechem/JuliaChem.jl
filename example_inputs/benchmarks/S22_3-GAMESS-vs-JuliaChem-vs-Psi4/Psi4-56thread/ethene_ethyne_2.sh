#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export MKL_NUM_THREADS=1
#
julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes 6-311++G_2d_2p/ethene_ethyne_2.jl
