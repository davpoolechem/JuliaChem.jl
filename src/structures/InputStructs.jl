"""
    Ctrl_Flags
Summary
======
Core flags for controlling overall calculation

Flags
======
1. RUNTYP::String = Determines calculation type (no default setting)
"""

struct Ctrl_Flags
    NAME::String
end
export Ctrl_Flags

"""
    struct Basis_Flags
Flags which control infomation about basis set information and orbital
occupancy. These flags should be set in the input_basis_flags() function
in the input file.

The flags are as follows:
1. NORB = Number of orbitals in the system overall (no default value)
2. NOCC = Number of doubly occupied orbitals (no default value)
"""
struct Basis_Flags
    NORB::Int64
    NOCC::Int64
end
export Basis_Flags
