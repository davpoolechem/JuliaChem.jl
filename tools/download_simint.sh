#/bin/bash

#== download simint code generator ==#
cd /export/home/david/projects/Julia/JuliaChem.jl/tools/travis
git clone https://github.com/simint-chem/simint-generator.git

#== run simint code generator ==#
cd simint-generator
mkdir build
cd build
cmake ../
make 
cd ..
./create.py -g build/generator/ostei -l 2 -p 2 /export/home/david/projects/Julia/JuliaChem.jl/tools/travis/simint

#== build simint ==#
cd /export/home/david/projects/Julia/JuliaChem.jl/tools/travis/simint
mkdir build
cd build
cmake -DSIMINT_VECTOR=scalar-sse -DCMAKE_INSTALL_PREFIX=/export/home/david/projects/Julia/JuliaChem.jl/tools/travis/simint-install -DCMAKE_C_FLAGS=-fPIC ../
make
make install 
