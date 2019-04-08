MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

DO_MPIRUN='mpirun'
DO_JULIARUN='~/Local/Documents/Programs/julia-1.1.0/julia'
JULIA_NUM_THREADS=${2} $DO_MPIRUN -np ${1} julia $SYSIMG $MODULE_OPT $CODE_OPT src/core/juliachem_shell.jl
