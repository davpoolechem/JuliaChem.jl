__precompile__(false)
module InputFunctions

using InputFile

using InputIntegrals
using InputStructs

function input_flags()
    #CTRL::Ctrl_Flags = input_ctrl_flags()
    BASIS::Basis_Flags = input_basis_flags()
    HF::HF_Flags = input_hf_flags()

    #FLAGS::Flags = Flags(CTRL,BASIS,HF)
    FLAGS::Flags = Flags(BASIS,HF)
    return FLAGS
end
export input_flags

function input_coord()
    #CTRL::Ctrl_Flags = input_ctrl_flags()
    coord::Array{Float64,2} = input_geometry()

    return coord
end
export input_coord

function read_in_enuc()
    enuc::Float64 = input_enuc()

    return enuc
end
export read_in_enuc

function read_in_ovr()
    ovr::Array{Float64,2} = input_ovr()

    ovr_matrix::Array{Float64,2} = get_oei_matrix(ovr)
    return ovr_matrix
end
export read_in_ovr

function read_in_kei()
    kei::Array{Float64,2} = input_kei()

    kei_matrix::Array{Float64,2} = get_oei_matrix(kei)
    return kei_matrix
end
export read_in_kei

function read_in_nai()
    nai::Array{Float64,2} = input_nai()

    nai_matrix::Array{Float64,2} = get_oei_matrix(nai)
    return nai_matrix
end
export read_in_nai

function read_in_tei()
    tei::Array{Float64,2} = input_tei()

    tei_array::Array{Float64,1} = get_tei_matrix(tei)
    return tei_array
end
export read_in_tei

end
