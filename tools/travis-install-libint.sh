#/bin/bash

#== download libint code generator ==#
cd /home/travis 
git clone --quiet https://github.com/evaleev/libint.git

#== install libint ==#
cd libint
./configure CPPFLAGS="-I$EIGEN_ROOT/include" --with-boost=$BOOST_ROOT --enable-shared=yes --prefix=/home/travis/libint-install
make -j2 -s
make -j2 install
