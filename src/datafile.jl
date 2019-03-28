#------------------------------#
#             IO.jl            #
#------------------------------#
"""
    process_data_file(data::String)
Summary
======
Read in data file and translate it to an data file object.

Arguments
======
data = name of data file to read in
"""
function process_data_file(data::String)
    dat::IOStream = open(data)
    file::Array{String,1} = readlines(data)
    close(dat)

    return file
end

"""
    read_in_enuc(data::Array{String,1})
Summary
======
Extract nuclear repulsion energy value from data file object.

Arguments
======
data = name of data file object to process
"""
function read_in_enuc(data::Array{String,1})
    enuc_loc::Int64 = findnext(data.=="\$\$ENUC",1)
    enuc::Float64 = parse(Float64, data[enuc_loc+1])
    return enuc
end

"""
    read_in_oei(data::Array{String,1}, type::String)
Summary
======
Extract one-electron integrals from data file object. Type dictates which
one-electron integrals are returned.

Arguments
======
data = name of data file object to process

type = tag for determining which type of one-electron integrals are returned:
1. type=OVR returns overlap integrals.
2. type=KEI returns kinetic energy integrals.
3. type=NAI returns nuclear-electron attraction integrals.
"""
function read_in_oei(data::Array{String,1}, type::String)
    loc::Int64 = findnext(data.=="\$\$$type",1)
    nbf::Int64 = 7
    nbf2::Int64 = nbf*(nbf+1)/2

    oei = Matrix{Float64}(undef,nbf,nbf)
    ij::Int64 = 1;
    for i::Int64 in 1:1:nbf, j::Int64 in 1:1:i
        ij = i*(i-1)/2 + j
        oei[i,j] = parse(Float64, split(data[ij+loc])[3])
        oei[j,i] = oei[i,j]
    end

    return oei
end

"""
    read_in_tei(data::Array{String,1})
Summary
======
Extract two-electron integrals from data file object.

Arguments
======
data = name of data file object to process
"""
function read_in_tei(data::Array{String,1})
    loc::Int64 = findnext(data.=="\$\$TEI",1)
    nint::Int64 = 228

    tei::Array{Float64,1} = zeros(2401)
    for index::Int64 in 1:1:nint
        i::Int64 = parse(Int64, split(data[index+loc])[1])
        j::Int64 = parse(Int64, split(data[index+loc])[2])
        k::Int64 = parse(Int64, split(data[index+loc])[3])
        l::Int64 = parse(Int64, split(data[index+loc])[4])

        ij::Int64 = (i > j) ? i*(i+1)/2 + j : j*(j+1)/2 + i
        kl::Int64 = (k > l) ? k*(k+1)/2 + l : l*(l+1)/2 + k
        ijkl::Int64 = (ij > kl) ? ij*(ij+1)/2 + kl : kl*(kl+1)/2 + ij

        tei[ijkl] = parse(Float64, split(data[index+loc])[5])
    end

    return tei
end

"""
    display_data_file(data::Array{String,1})
Summary
======
Re-print data file to display.

Arguments
======
data = name of data file object to process
"""
function display_data_file(data::String)
    dat::IOStream = open(data)

    file::Array{String,1} = readlines(data)

    for line in file
        println(line)
    end

    close(dat)
end
