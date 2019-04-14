MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""
DO_MPIRUN='/opt/pgi/linux86-64/2018/mpi/openmpi/bin/mpirun'
DO_JULIARUN='~/Local/Documents/Programs/julia-1.1.0/julia'

echo "Printing input file..."
~/Local/Documents/Programs/julia-1.1.0/julia -e "import JCInputFile; JCInputFile.assign(\"${2}\")"
echo "Performing calculation..."
JULIA_NUM_THREADS=${4} $DO_MPIRUN -np ${3} ~/Local/Documents/Programs/julia-1.1.0/julia $SYSIMG $MODULE_OPT $CODE_OPT ${1}