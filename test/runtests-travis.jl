import Pkg; Pkg.add("Test")
import Test

include("../example_scripts/minimal-rhf-travis.jl")

#== select input files ==#
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
  scf = minimal_rhf(inputs[3])
  Test.@test scf["Energy"] ≈ -377.5889420312

  scf = minimal_rhf(inputs[10])
  Test.@test scf["Energy"] ≈ -270.9302743515 

  scf = minimal_rhf(inputs[18])
  Test.@test scf["Energy"] ≈ -286.9286565278 
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
