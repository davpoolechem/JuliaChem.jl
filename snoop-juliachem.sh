#/bin/bash

#perform snoopcompile work
echo "----------------------------------------"
echo "    PERFORMING SNOOPCOMPILE ANALYSIS    "
echo "----------------------------------------"
if [ -d src/core/snoop ]; then
   rm -rf snoop
   rm snoop.csv
fi
JULIA_NUM_THREADS=1 /opt/pgi/linux86-64/2018/mpi/openmpi/bin/mpirun -np 1 ~/Local/Documents/Programs/julia-1.1.0/julia src/snoop.jl example_scripts/minimal_rhf_script.jl example_inputs/sto3g-water.json
