mutable struct Ctrl_Flags
    RUNTYP::String
end

mutable struct Basis_Flags
    NORB::Int64
    NOCC::Int64
end

mutable struct HF_Flags
    NITER::Int64
    DELE::Float64
    RMSD::Float64
end

mutable struct Flags
    CTRL::Ctrl_Flags
    BASIS::Basis_Flags
    HF::HF_Flags
end
