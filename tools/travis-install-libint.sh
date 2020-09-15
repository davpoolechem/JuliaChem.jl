#/bin/bash

#== download libint code generator ==#
cd /home/travis 
git clone --quiet https://github.com/evaleev/libint.git

#== install libint ==#
cd libint
./configure \
  --enable-shared=yes --prefix=/home/travis/libint-install \
  --enable-1body=0 --enable-eri=0 --with-max-am=3 \
  --with-multipole-max-order=3 --enable-eri3=no \
  --enable-eri2=no --enable-g12=no

make -j2 -s
set -x
sed -i "s/install\_cmake //g" Makefile
set +x
make -j2 install
