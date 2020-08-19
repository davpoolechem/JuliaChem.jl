#/bin/bash

#== download simint code generator ==#
cd /home/travis 
git clone https://github.com/simint-chem/simint-generator.git

#== run simint code generator ==#
cd simint-generator
mkdir build
cd build
cmake ../
make 
cd ..
./create.py -g build/generator/ostei -l 3 -p 3 /home/travis/simint 

#== build simint ==#
cd /home/travis/simint 
mkdir build
cd build
cmake -DSIMINT_VECTOR=scalar-sse -DCMAKE_INSTALL_PREFIX=/home/travis/simint-install -DCMAKE_C_FLAGS=-fPIC ../
make -j2
make -j2 install

#export SIMINT=/home/travis/simint-install 
