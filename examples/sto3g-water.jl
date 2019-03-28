Base.include(@__MODULE__,"../src/input_structs.jl")

#------------------------------#
#      Control input flags     #
#------------------------------#
function input_ctrl_flags()
    #Set flags
    runtyp::String = "HF"

    #Form CTRL struct
    CTRL::Ctrl_Flags = Ctrl_Flags(runtyp)
    return CTRL
end

function input_basis_flags()
    #Set flags
    norb::Int64 = 7
    nocc::Int64 = 5

    #Form BASIS struct
    BASIS::Basis_Flags = Basis_Flags(norb, nocc)
    return BASIS
end

function input_hf_flags()
    #Set flags
    niter::Int64 = 50
    dele::Float64 = 1E-8
    rmsd::Float64 = 1E-6

    #Form HF struct
    HF::HF_Flags = HF_Flags(niter, dele, rmsd)
    return HF
end

#------------------------------#
#    Control input geometry    #
#------------------------------#
function input_geometry()
    #Set geometry
    geometry::Array{Float64,2} =
    [ 8.000000000000   0.000000000000  -0.143225816552   0.000000000000;
      1.000000000000   1.638036840407   1.136548822547  -0.000000000000;
      1.000000000000  -1.638036840407   1.136548822547  -0.000000000000 ]

    return geometry
end
