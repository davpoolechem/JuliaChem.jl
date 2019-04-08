using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/colortypes_compiles.csv")

SnoopCompile.@snoop "snoop.csv" begin
    using Pkg
    include("juliachem.jl")
end

### Parse the compiles and generate precompilation scripts
# This can be run repeatedly to tweak the scripts

data = SnoopCompile.read("snoop.csv")

pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("snoop", pc)
