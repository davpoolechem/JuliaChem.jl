MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=no --math-mode=fast --inline=yes"
#SYSIMG="-J /home/david/.julia/packages/PackageCompiler/oT98U/sysimg/backup/native/sys.ji"
SYSIMG=""

~/Local/Documents/Programs/julia-1.1.0/julia $SYSIMG $MODULE_OPT $CODE_OPT julichem_shell.jl  
