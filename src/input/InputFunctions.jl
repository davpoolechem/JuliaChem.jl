include("../../InputFile.jl")

include("InputIntegrals.jl")

using JCStructs

function input_flags()
    CTRL::Ctrl_Flags = input_ctrl_flags()
    BASIS::Basis_Flags = input_basis_flags()
    HF::HF_Flags = input_hf_flags()

    FLAGS::Flags = Flags(CTRL,BASIS,HF)
    #FLAGS::Flags = Flags(BASIS,HF)
    return FLAGS
end

function input_coord()
    #CTRL::Ctrl_Flags = input_ctrl_flags()
    coord::Array{Float64,2} = input_geometry()

    return coord
end

function read_in_enuc(type::T) where {T<:Number}
    enuc::T = input_enuc()

    return enuc
end

function read_in_ovr(type::T) where {T<:Number}
    ovr::Array{T,2} = input_ovr()

    ovr_matrix::Array{T,2} = get_oei_matrix(ovr)
    return ovr_matrix
end

function read_in_kei(type::T) where {T<:Number}
    kei::Array{T,2} = input_kei()

    kei_matrix::Array{T,2} = get_oei_matrix(kei)
    return kei_matrix
end

function read_in_nai(type::T) where {T<:Number}
    nai::Array{T,2} = input_nai()

    nai_matrix::Array{T,2} = get_oei_matrix(nai)
    return nai_matrix
end

function read_in_tei(type::T) where {T<:Number}
    tei::Array{T,2} = input_tei()

    tei_array::Array{T,1} = get_tei_matrix(tei)
    return tei_array
end
