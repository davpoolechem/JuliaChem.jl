#!/bin/bash
#
source /etc/profile.d/modules.sh
source ~/.bashrc
#
export OMP_NUM_THREADS=112
#
cd /home/davpoolechem/shared/gms-dp/gamess
./rungms-dev S22_3/6-311++G_2d_2p/benzene_hcn_2.inp 00 1 1
