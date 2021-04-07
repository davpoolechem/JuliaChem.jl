#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
include("rhf-properties-xyz-repl.jl")

JuliaChem.initialize()
rhf_properties_xyz(ARGS[1])
JuliaChem.finalize()

