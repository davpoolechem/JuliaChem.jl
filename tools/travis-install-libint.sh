#/bin/bash

#== download libint code generator ==#
cd /home/travis 
git clone --quiet https://github.com/evaleev/libint.git

#== install libint ==#
cd libint
./configure \
  --enable-shared=yes --prefix=/home/travis/libint-install \
  --enable-1body=0 --enable-eri=0 --with-max-am=3
make -j2 -s
sed -i "s/install_cmake //g" Makefile
make -j2 -s install
