import Test

include("../tools/travis/travis-rhf.jl")
include("s22_gamess_values.jl")

#== select input files ==#
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

display(inputs)

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 calculations ==#
molecules = collect(2:3) 

s22_test_results = Dict([]) 
for imol in molecules 
  s22_test_results[imol] = travis_rhf(inputs[imol])
end

#== check energies ==#
Test.@testset "S22 Energy" begin
  for imol in molecules 
    Test.@test s22_test_results[imol][:Energy]["Energy"] ≈ S22_GAMESS[imol]["Energy"]
  end
end

#== check dipole moments ==#
Test.@testset "S22 Dipoles" begin
  for imol in molecules 
    if S22_GAMESS[imol]["Dipole"] == 1.0E-6
      Test.@test abs(s22_test_results[imol][:Properties]["Dipole"][:moment]) <= S22_GAMESS[imol]["Dipole"] #check if approximately zero 
    else
      Test.@test s22_test_results[imol][:Properties]["Dipole"][:moment] ≈ S22_GAMESS[imol]["Dipole"] atol=5.0E-5
    end
  end
end

#== check HOMO-LUMO gaps ==#
Test.@testset "S22 HOMO-LUMO Gaps" begin
  for imol in molecules 
    Test.@test s22_test_results[imol][:Properties]["MO Energies"][:homo_lumo] ≈ S22_GAMESS[imol]["HOMO-LUMO Gap"] atol=5.0E-4
  end
end

#== check Mulliken charges ==#
Test.@testset "S22 Mulliken Charges" begin
  for imol in molecules 
    Test.@test s22_test_results[imol][:Properties]["Mulliken Population"] ≈ 
      S22_GAMESS[imol]["Mulliken Population"] atol=5.0E-6
  end
end


#== finalize JuliaChem ==#
JuliaChem.finalize()       
