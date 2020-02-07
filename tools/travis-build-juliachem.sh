#/bin/bash

cd ./deps
julia build.jl
cd ../
julia tools/travis-build-juliachem.jl
