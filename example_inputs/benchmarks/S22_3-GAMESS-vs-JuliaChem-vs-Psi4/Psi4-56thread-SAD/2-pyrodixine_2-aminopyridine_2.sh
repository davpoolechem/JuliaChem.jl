#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=1
#
julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes 6-311++G_2d_2p/2-pyrodixine_2-aminopyridine_2.jl
