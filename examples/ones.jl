Base.include(@__MODULE__,"../src/input/input_structs.jl")

#------------------------------#
#      Control input flags     #
#------------------------------#

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
    norb::Int64 = 1
    nocc::Int64 = 1

    #Form BASIS struct
    BASIS::Basis_Flags = Basis_Flags(norb, nocc)
    return BASIS
end

function input_hf_flags()
    #Set flags
    niter::Int64 = 50
    dele::Float64 = 1E-6
    rmsd::Float64 = 1E-4

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
    [ 8.000000000000   0.000000000000  -0.143225816552   0.000000000000; ]

    return geometry
end

#------------------------------#
#   Control data input info    #
#------------------------------#
function input_enuc()
    enuc::Float64 = 8.002367061810450

    return enuc
end

function input_ovr()
    ovr::Array{Float64,2} =
    [ 1     1    1.000000000000000; ]

    return ovr
end

function input_kei()
    kei::Array{Float64,2} =
    [ 1     1   29.003199945539588 ]

    return kei
end

function input_nai()
    nai::Array{Float64,2} =
    [ 1     1  -61.580595358149914; ]

    return nai
end

function input_tei()
    tei::Array{Float64,2} =
    [ 1     1     1     1    4.785065404705506; ]

    return tei
end
