#=
"""
    Data
Summary
======
Core Hartree-Fock data structures

Fields
======
1. Fock::Array{Float64,2} = Fock Matrix
2. Density::Array{Float64,2} = Density Matrix
3. Coeff::Array{Float64,2} = Molecular Orbital Coefficient Matrix
4. Energy::Float64 = Electronic Energy
"""
=#
mutable struct Data{T<:Number}
    Fock::Array{T,2}
    Density::Array{T,2}
    Coeff::Array{T,2}
    Energy::T
end
export Data

mutable struct RHFRestartData
    H::Array{Float64,2}
    Ortho::Array{Float64,2}
    iter::UInt32

    Fock::Array{Float64,2}
    Density::Array{Float64,2}
    Coeff::Array{Float64,2}
    Energy::Float64
end
export RHFRestartData

"""
    HF_Flags

Flags which control the minutiae of a Hartree-Fock calculation. These flags
should be set in the input_hf_flags() function in the input file.

The flags are as follows:
1. NITER = Maximum number of SCF iterations (default = 50)
2. DELE = Change-in-energy convergence threshold (default = 1E-8)
3. RMSD = Change-in-root-mean-square-density convergence threshold (default = 1E-6)
"""
struct SCF_Flags
    NITER::UInt32
    DELE::Float64
    RMSD::Float64
    PREC::String
    DIRECT::Bool
    DEBUG::Bool
end
export SCF_Flags

#=
"""
    Flags
Summary
======
Conglomeration of other flag fields

Flag fields
======
1. CTRL::Ctrl_Flags = Core flags for controlling overall calculation
2. BASIS::Basis_Flags = Flags setting basis set information
3. HF::HF_Flags = Flags relevant to a Hartree-Fock calculation
"""
=#

struct RHF_Flags
    CTRL::Ctrl_Flags
    BASIS::Basis_Flags
    SCF::SCF_Flags
end
export RHF_Flags
