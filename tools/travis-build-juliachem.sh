#/bin/bash

cd ./deps
julia build.jl
cd ../
julia tools/travis_build_juliachem.jl
