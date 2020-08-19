#/bin/bash

#== download openmpi ==#
cd /home/travis
wget https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-2.1.6.tar.gz
tar -xzvf openmpi-2.1.6.tar.gz

#== build openmpi ==# 
cd openmpi-2.1.6
./configure --prefix=/home/travis/openmpi
make -j2
make -j2 install

