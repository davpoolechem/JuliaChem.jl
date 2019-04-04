#/bin/bash

#perform snoopcompile work
echo "----------------------------------------" 
echo "    PERFORMING SNOOPCOMPILE ANALYSIS    "
echo "----------------------------------------" 
if [ -d snoop ]; then
   rm -rf snoop
   rm snoop.csv
fi
julia snoop.jl

#do actual build
echo ""
echo "----------------------------------------" 
echo "           BUILDING JULICHEM            "
echo "----------------------------------------" 
if [ -d builddir ]; then
   rm -rf builddir
fi

MODULE_OPT="--compile=all --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O3 --inline=yes --check-bounds=no --math-mode=fast"
CC_COMP="--cc=gcc" 
julia ~/.julia/packages/PackageCompiler/oT98U/juliac.jl -Re $MODULE_OPT $CODE_OPT $CC_COMP julichem.jl julichem.c 
