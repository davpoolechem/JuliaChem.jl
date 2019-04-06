MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

DO_MPIRUN='/opt/pgi/linux86-64/2018/mpi/openmpi/bin/mpirun'
DO_JULIARUN='~/Local/Documents/Programs/julia-1.1.0/julia'
JULIA_NUM_THREADS=${1} $DO_MPIRUN ~/Local/Documents/Programs/julia-1.1.0/julia $SYSIMG $MODULE_OPT $CODE_OPT src/core/julichem_shell.jl
