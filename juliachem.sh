MODULE_OPT=" " 
CODE_OPT="-O2 --check-bounds=yes --math-mode=ieee"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

echo "Performing calculation..."
JULIA_NUM_THREADS=${4} mpirun -np ${3} julia $MODULE_OPT $CODE_OPT ${1} ${2}
