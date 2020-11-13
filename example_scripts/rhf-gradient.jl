#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
include("rhf-gradient-repl.jl")

JuliaChem.initialize()
rhf_gradient(ARGS[1])
JuliaChem.finalize()

