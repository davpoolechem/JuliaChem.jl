#=============================#
#== put needed modules here ==#
#=============================#
@everywhere import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
include("minimal-rhf-repl.jl")
#include("test.jl")

JuliaChem.initialize()
minimal_rhf(ARGS[1])
JuliaChem.finalize()

