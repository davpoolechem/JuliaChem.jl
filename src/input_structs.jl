struct Ctrl_Flags
    RUNTYP::String
end

struct Basis_Flags
    NORB::Int64
    NOCC::Int64
end

struct HF_Flags
    NITER::Int64
    DELE::Float64
    RMSD::Float64
end

struct Flags
    CTRL::Ctrl_Flags
    BASIS::Basis_Flags
    HF::HF_Flags
end
