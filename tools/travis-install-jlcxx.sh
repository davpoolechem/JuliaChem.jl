#/bin/bash

#== download jlcxx ==#
cd /home/travis 
git clone --quiet https://github.com/JuliaInterop/libcxxwrap-julia.git

#== install jlcxx ==#
cd libcxxwrap-julia 
mkdir build
cd build
cmake -DJulia_PREFIX=${JULIA_ROOT} -DCMAKE_INSTALL_PREFIX=/home/travis/jlcxx-install ../
make -j2 -s 
make -j2 -s install
