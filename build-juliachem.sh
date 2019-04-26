#/bin/bash

#perform snoopcompile work
echo "----------------------------------------"
echo "    PERFORMING SNOOPCOMPILE ANALYSIS    "
echo "----------------------------------------"
if [ -d src/core/snoop ]; then
   rm -rf snoop
   rm snoop.csv
fi
/opt/pgi/linux86-64/2018/mpi/openmpi/bin/mpirun -np 1 ~/Local/Documents/Programs/julia-1.1.0/julia snoop.jl

#do actual build
#echo ""
#echo "-----------------------------------------"
#echo "           BUILDING JULIACHEM            "
#echo "-----------------------------------------"
#if [ -d builddir ]; then
#   rm -rf builddir
#fi

MODULE_OPT="--compile=all --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O3 --inline=yes --check-bounds=no --math-mode=fast"
CC_COMP="--cc=gcc --startup-file=yes"
julia ~/.julia/packages/PackageCompiler/oT98U/juliac.jl -Rs $MODULE_OPT $CODE_OPT $CC_COMP example_scripts/minimal-rhf-script.jl
