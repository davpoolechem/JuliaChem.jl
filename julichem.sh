MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

JULIA_NUM_THREADS=${2} ~/Local/Documents/Programs/julia-1.1.0/julia $SYSIMG $MODULE_OPT $CODE_OPT src/core/julichem_shell.jl ${1}
