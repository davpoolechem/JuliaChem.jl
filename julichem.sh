MODULE_OPT=" --sysimage-native-code=yes --compiled-modules=yes"
CODE_OPT="-O0 --check-bounds=no --math-mode=fast"

julia $MODULE_OPT $CODE_OPT julichem_shell.jl $1 
