include("./julichem.jl")
import Main.JuliChem

#test script
inp = JuliChem.input("test.inp")
scf = JuliChem.scf(inp)
#JuliChem.properties(scf)

#inp = io.readin("test.inp")
#coord = io.geomin(inp)
