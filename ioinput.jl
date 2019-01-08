#module for io operations
module ioinput
    #individual atom coordinate
    mutable struct Atom
        number::Int64
        x::Float64
        y::Float64
        z::Float64

        Atom() = new(0,0,0,0)
        Atom(number,x,y,z) = new(number,x,y,z)
    end

    #array of atoms and their coordinates
    mutable struct Coords
        Atoms::Array{Atom,1}
    end

    #contains structures
    #mutable struct Flags
#        Ctrl
    #    Basis
#        HF
#    end

    #contains information of readin information
    mutable struct ReadIn
        TEI::Array{Float64,1}
        KEI::Array{Float64,1}
        NAI::Array{Float64,1}
        OVR::Array{Float64,1}
        Enuc::Float64

        ReadIn() = new(Array(undef,1),Array(undef,1),Array(undef,1),Array(undef,1), 0.0)
    end

    mutable struct Input
        Coords::Coords
        #Flags::Flags
        ReadIn::ReadIn
    end
end
