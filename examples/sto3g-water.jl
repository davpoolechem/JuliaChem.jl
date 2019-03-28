Base.include(@__MODULE__,"../src/io_structs.jl")

function input_flags()
    #CTRL
    runtyp::String = "HF"

    CTRL::Ctrl_Flags = Ctrl_Flags(runtyp)

    #BASIS
    norb::Int64 = 7
    nocc::Int64 = 5

    BASIS::Basis_Flags = Basis_Flags(norb, nocc)

    #HF
    niter::Int64 = 50
    dele::Float64 = 1E-8
    rmsd::Float64 = 1E-6

    HF::HF_Flags = HF_Flags(niter, dele, rmsd)

    #overall flags
    FLAGS::Flags = Flags(CTRL,BASIS,HF)
    return FLAGS
end

function input_geometry()
    geometry::Array{Float64,2} = [ 8.000000000000   0.000000000000  -0.143225816552   0.000000000000;
                                   1.000000000000   1.638036840407   1.136548822547  -0.000000000000;
                                   1.000000000000  -1.638036840407   1.136548822547  -0.000000000000 ]
    return geometry
end
