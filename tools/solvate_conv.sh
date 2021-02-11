#!/bin/bash

#== set up env variables ##=
export BABEL_ROOT=$HOME/programs/install/openbabel
export LD_LIBRARY_PATH=$BABEL_ROOT/lib64:$BABEL_ROOT/lib:$LD_LIBRARY_PATH

#== create pdf file of water cluster ==#
if [ ! -f water_$1\.pdb ] 
then
  ./solvate -t $1 water_$1
fi

#== convert pdb to xyz ==#
touch water_$1\.xyz
$BABEL_ROOT/bin/obabel water_$1\.pdb -O water_$1\.xyz

#== convert pdb to inp ==#
touch water_$1\.inp
$BABEL_ROOT/bin/obabel water_$1\.pdb -O water_$1\.inp
