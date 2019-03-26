Base.include(Main,"julichem.jl")

Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    inp::String = ARGS[1]
    scf::Data = exe(inp)
    return 0
end
