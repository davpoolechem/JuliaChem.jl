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
julia ~/.julia/packages/PackageCompiler/oT98U/juliac.jl -Re julichem.jl
