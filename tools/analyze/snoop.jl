using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/colortypes_compiles.csv")

SnoopCompile.@snoopc "snoop.csv" begin
    using Pkg
    include("src/JuliaChem.jl")
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

data = SnoopCompile.read("snoop.csv")

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("snoop", pc)
