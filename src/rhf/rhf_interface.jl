Base.include(@__MODULE__,"rhf_scf.jl")

"""
    do_rhf(dat::Array{String,1})
Summary
======
Execute the JuliChem RHF algorithm.

Arguments
======
dat = input data file object
"""
function do_rhf(dat::Array{String,1},flags::Flags)
    #determine and perform proper method
    GC.enable(false)
    scf::Data = rhf_energy(dat,flags)
    GC.enable(true)
    GC.gc()

    return scf
end
