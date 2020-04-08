#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
include("minimal-rhf-repl.jl")

JuliaChem.initialize()
minimal_rhf(ARGS[1])
JuliaChem.finalize()

