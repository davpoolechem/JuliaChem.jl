"""
    Ctrl_Flags
Summary
======
Core flags for controlling overall calculation

Flags
======
1. RUNTYP::String = Determines calculation type (no default setting)
"""
#struct Ctrl_Flags
#    RUNTYP::String
#end

"""
    Basis_Flags
Summary
======
Flags setting basis set information

Flags
======
1. NORB = Number of orbitals overall (no default value)
2. NOCC = Number of occupied orbitals (no default value)
"""
struct Basis_Flags
    NORB::Int64
    NOCC::Int64
end

"""
    HF_Flags
Summary
======
Flags relevant to a Hartree-Fock calculation

Flags
======
1. NITER = Maximum number of SCF iterations (default = 50)
2. DELE = Change-in-energy convergence threshold (default = 1E-8)
3. RMSD = Change-in-root-mean-square-density convergence threshold (default = 1E-6)
"""
struct HF_Flags
    NITER::Int64
    DELE::Float64
    RMSD::Float64
end

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
struct Flags
#    CTRL::Ctrl_Flags
    BASIS::Basis_Flags
    HF::HF_Flags
end
