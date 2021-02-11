#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export OMP_NUM_THREADS=84
export OMP_STACKSIZE=6000M
#
ulimit -s unlimited
cd /home/davpoolechem/shared/gms-dp/gamess
./rungms-dev S22_3/6-311++G_2d_2p/at2/adenine_thymine_2.inp 00 1 1 
