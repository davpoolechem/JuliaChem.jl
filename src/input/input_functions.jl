Base.include(@__MODULE__, "../../input.jl")
Base.include(@__MODULE__, "integrals.jl")

function input_flags()
    #CTRL::Ctrl_Flags = input_ctrl_flags()
    BASIS::Basis_Flags = input_basis_flags()
    HF::HF_Flags = input_hf_flags()

    #FLAGS::Flags = Flags(CTRL,BASIS,HF)
    FLAGS::Flags = Flags(BASIS,HF)
    return FLAGS
end

function read_in_enuc()
    enuc::Float64 = input_enuc()

    return enuc
end

function read_in_ovr()
    ovr::Array{Float64,2} = input_ovr()

    ovr_matrix::Array{Float64,2} = read_in_oei(ovr)
    return ovr_matrix
end

function read_in_kei()
    kei::Array{Float64,2} = input_kei()

    kei_matrix::Array{Float64,2} = read_in_oei(kei)
    return kei_matrix
end

function read_in_nai()
    nai::Array{Float64,2} = input_nai()

    nai_matrix::Array{Float64,2} = read_in_oei(nai)
    return nai_matrix
end

function read_in_tei()
    tei::Array{Float64,2} = input_tei()

    tei_array::Array{Float64,1} = read_in_tei(tei)
    return tei_array
end
