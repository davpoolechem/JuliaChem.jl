#/bin/bash

#== download boost ==#
cd /home/travis
wget https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.gz
tar -xzf boost_1_74_0.tar.gz

#== install boost ==#
cd boost_1_74_0
./bootstrap.sh
./b2 install --prefix=/home/travis/boost-install -d0
