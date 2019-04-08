#/bin/bash

#perform snoopcompile work
echo "----------------------------------------"
echo "    PERFORMING SNOOPCOMPILE ANALYSIS    "
echo "----------------------------------------"
if [ -d src/core/snoop ]; then
   rm -rf src/core/snoop
   rm src/core/snoop.csv
fi
julia src/core/snoop.jl

#do actual build
echo ""
echo "-----------------------------------------"
echo "           BUILDING JULIACHEM            "
echo "-----------------------------------------"
if [ -d builddir ]; then
   rm -rf builddir
fi

MODULE_OPT="--compile=all --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O3 --inline=yes --check-bounds=no --math-mode=fast"
CC_COMP="--cc=gcc"
julia ~/.julia/packages/PackageCompiler/oT98U/juliac.jl -Re $MODULE_OPT $CODE_OPT $CC_COMP src/core/juliachem.jl src/core/juliachem.c
