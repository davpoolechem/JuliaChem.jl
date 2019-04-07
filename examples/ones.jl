using InputStructs

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
export input_basis_flags

function input_hf_flags()
    #Set flags
    niter::Int64 = 50
    dele::Float64 = 1E-6
    rmsd::Float64 = 1E-4

    #Form HF struct
    HF::HF_Flags = HF_Flags(niter, dele, rmsd)
    return HF
end
export input_hf_flags

#------------------------------#
#    Control input geometry    #
#------------------------------#
function input_geometry()
    #Set geometry
    geometry::Array{Float64,2} =
    [ 8.000000000000   0.000000000000  -0.143225816552   0.000000000000; ]

    return geometry
end
export input_geometry

#------------------------------#
#   Control data input info    #
#------------------------------#
function input_enuc()
    enuc::Float64 = 8.002367061810450

    return enuc
end
export input_enuc

function input_ovr()
    ovr::Array{Float64,2} =
    [ 1     1    1.000000000000000; ]

    return ovr
end
export input_ovr

function input_kei()
    kei::Array{Float64,2} =
    [ 1     1   29.003199945539588 ]

    return kei
end
export input_kei

function input_nai()
    nai::Array{Float64,2} =
    [ 1     1  -61.580595358149914; ]

    return nai
end
export input_nai

function input_tei()
    tei::Array{Float64,2} =
    [ 1     1     1     1    4.785065404705506; ]

    return tei
end
export input_tei
