#/bin/bash

#== download openmpi ==#
cd /home/travis 
wget https://download.open-mpi.org/release/open-mpi/v3.0/openmpi-3.0.5.tar.gz
tar -xzvf openmpi-3.0.5.tar.gz

#== build openmpi ==# 
cd openmpi-3.0.5
./configure --prefix=/home/travis/openmpi
make
make install

