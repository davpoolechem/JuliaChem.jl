import Pkg; Pkg.add("Test")
import Test

include("../tools/travis-rhf.jl")

#== select input files ==#
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
  rhf_energy, properties = travis_rhf(inputs[3])                                
  Test.@test rhf_energy["Energy"] ≈ -377.5889420312                             
  Test.@test properties["Dipole"][:moment] <= 1.0E-10 #approximately 0.0 
 
  #scf = travis_rhf(inputs[10])
  #Test.@test scf["Energy"] ≈ -270.9302743515 

  #scf = travis_rhf(inputs[18])
  #Test.@test scf["Energy"] ≈ -286.9286565278 
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
