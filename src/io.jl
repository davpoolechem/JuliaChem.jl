#------------------------------#
#             IO.jl            #
#------------------------------#
"""
    readin(input::String)
Summary
======
Read in input file and translate it to an input file object.

Arguments
======
input = name of input file to read in
"""
function readin(input::String)
    inp::IOStream = open(input)
    file::Array{String,1} = readlines(input)
    close(inp)

    return file
end

"""
    geomin(input::Array{String,1})
Summary
======
Extract geometric coordinates from input file.

Arguments
======
input = name of input file object to process
"""
function geomin(input::Array{String,1})
    geom::Int64 = findnext(input.=="\$\$GEOM",1)
    natoms::Int64 = parse(Int64, input[geom+1])

    coord::Matrix{Float64} = Matrix{Float64}(undef,natoms,4)
    for i::Int64 in [1,2...,natoms]
        coord[i,:] = parse.(Float64, split(input[i+geom+1]))
    end

    return coord
end

#extract coordinates from input file
#function flagsin(input::Array{String,1})
#end

"""
    enucin(input::Array{String,1})
Summary
======
Extract nuclear repulsion energy value from input data file.

Arguments
======
input = name of input file object to process
"""
function enucin(input::Array{String,1})
    enuc_loc::Int64 = findnext(input.=="\$\$ENUC",1)
    enuc::Float64 = parse(Float64, input[enuc_loc+1])
    return enuc
end

"""
    oeiin(input::Array{String,1}, type::String)
Summary
======
Extract one-electron integrals from input data file. Type dictates which
one-electron integrals are returned.

Arguments
======
input = name of input file object to process

type = tag for determining which type of one-electron integrals are returned:
1. type=OVR returns overlap integrals.
2. type=KEI returns kinetic energy integrals.
3. type=NAI returns nuclear-electron attraction integrals.
"""
function oeiin(input::Array{String,1}, type::String)
    loc::Int64 = findnext(input.=="\$\$$type",1)
    nbf::Int64 = 7
    nbf2::Int64 = nbf*(nbf+1)/2

    oei = Matrix{Float64}(undef,nbf,nbf)
    ij::Int64 = 1;
    for i::Int64 in 1:1:nbf, j::Int64 in 1:1:i
        ij = i*(i-1)/2 + j
        oei[i,j] = parse(Float64, split(input[ij+loc])[3])
        oei[j,i] = oei[i,j]
    end

    return oei
end

"""
    teiin(input::Array{String,1})
Summary
======
Extract two-electron integrals from input data file.

Arguments
======
input = name of input file object to process
"""
function teiin(input::Array{String,1})
    loc::Int64 = findnext(input.=="\$\$TEI",1)
    nint::Int64 = 228

    tei::Array{Float64,1} = zeros(2401)
    for index::Int64 in 1:1:nint
        i::Int64 = parse(Int64, split(input[index+loc])[1])
        j::Int64 = parse(Int64, split(input[index+loc])[2])
        k::Int64 = parse(Int64, split(input[index+loc])[3])
        l::Int64 = parse(Int64, split(input[index+loc])[4])

        ij::Int64 = (i > j) ? i*(i+1)/2 + j : j*(j+1)/2 + i
        kl::Int64 = (k > l) ? k*(k+1)/2 + l : l*(l+1)/2 + k
        ijkl::Int64 = (ij > kl) ? ij*(ij+1)/2 + kl : kl*(kl+1)/2 + ij

        tei[ijkl] = parse(Float64, split(input[index+loc])[5])
    end

    return tei
end

"""
    teiin(input::Array{String,1})
Summary
======
Re-print input file to display.

Arguments
======
input = name of input file object to process
"""
function iocat(input::String)
    inp::IOStream = open(input)

    file::Array{String,1} = readlines(input)

    for line in file
        println(line)
    end

    #flags::Int64 = findnext(file.=="\$FLAGS",1)
    #geom::Int64 = findnext(file.=="\$GEOM",flags)
    #vec::Int64 = findnext(file.=="\$VEC",geom)

    close(inp)
end

#function compile()
#    precompile(io.readin, tuple(String))
#    precompile(io.geomin, tuple(Array{String,1}))
#    precompile(io.enucin, tuple(Array{String,1}))
#    precompile(io.oeiin, tuple(Array{String,1}, String))
#    precompile(io.teiin, tuple(Array{String,1}))
#end

#we want to precompile all involved modules to reduce cold runs
#include("./snoop/precompile_Base.jl")
#_precompile_base()
#include("./snoop/precompile_Core.jl")
#_precompile_core()
#include("./snoop/precompile_hf.jl")
#_precompile_()
#include("./snoop/precompile_io.jl")
#_precompile_()
#include("./snoop/precompile_LinearAlgebra.jl")
#_precompile_()
#include("./snoop/precompile_openchem.jl")
#_precompile_()
#include("./snoop/precompile_unknown.jl")
#_precompile_()