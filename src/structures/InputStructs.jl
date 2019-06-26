"""
  struct Ctrl_Flags

Core flags for controlling overall calculation. These flags should be set via
an "Input":"Control Flags" section in the input file.

The flags are as follows:
1. name = Determines calculation label (no default setting)
"""
struct Ctrl_Flags
  NAME::String
end
export Ctrl_Flags

"""
  struct Basis_Flags
Flags which control infomation about basis set information and orbital
occupancy. These flags should be set via an "Input":"Basis Flags" section
in the input file.

The flags are as follows:
1. norb = Number of orbitals in the system overall (no default value)
2. nocc = Number of doubly occupied orbitals (no default value)
3. shells = Angular momentum of each shell in system (no default value)
"""
struct Basis_Flags
  NORB::Int64
  NOCC::Int64
end
export Basis_Flags
