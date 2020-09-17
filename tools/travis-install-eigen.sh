#/bin/bash

#== download eigen ==#
cd /home/travis 
wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
tar -xzf eigen-3.3.7.tar.gz

#== install eigen ==#
cd eigen-3.3.7
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/travis/eigen-install ../
make -j2 -s install
