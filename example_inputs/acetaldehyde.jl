Base.include(@__MODULE__,"../src/input/input_structs.jl")

#------------------------------#
#      Control input flags     #
#------------------------------#
#function input_ctrl_flags()
    #Set flags
#    runtyp::String = "HF"

    #Form CTRL struct
#    CTRL::Ctrl_Flags = Ctrl_Flags(runtyp)
#    return CTRL
#end

function input_basis_flags()
    #Set flags
    norb::Int64 = 9
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
    [ 6.0  0.000000000000     0.000000000000     0.000000000000
      6.0  0.000000000000     0.000000000000     2.845112131228
      8.0  1.899115961744     0.000000000000     4.139062527233
      1.0 -1.894048308506     0.000000000000     3.747688672216
      1.0  1.942500819960     0.000000000000    -0.701145981971
      1.0 -1.007295466862    -1.669971842687    -0.705916966833
      1.0 -1.007295466862     1.669971842687    -0.705916966833 ]

    return geometry
end
