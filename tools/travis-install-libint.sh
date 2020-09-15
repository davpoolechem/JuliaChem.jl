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

#== download boost ==#
cd /home/travis
wget https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.gz
tar -xzf boost_1_74_0.tar.gz

#== install boost ==#
cd boost_1_74_0
./bootstrap.sh
./b2 install --prefix=/home/travis/boost-install -d0

#== download libint code generator ==#
cd /home/travis 
git clone --quiet https://github.com/evaleev/libint.git

#== install libint ==#
cd libint
./configure CPPFLAGS="-I$EIGEN_ROOT/include" --with-boost=$BOOST_ROOT --enable-shared=yes --prefix=/home/travis/libint-install
make -j2 -s
make -j2 install
