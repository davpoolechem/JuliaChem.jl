MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=yes --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

JULIA_NUM_THREADS=${3} julia --procs=${2} $SYSIMG $MODULE_OPT $CODE_OPT julichem_shell.jl ${1} 
