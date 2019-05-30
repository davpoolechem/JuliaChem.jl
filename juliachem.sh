MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

echo "Performing calculation..."
JULIA_NUM_THREADS=${4} mpirun -np ${3} julia $SYSIMG $MODULE_OPT $CODE_OPT ${1} ${2}
