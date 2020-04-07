import Test

include("../example_scripts/minimal-rhf-repl.jl")

#== select input files ==#:w
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
  scf = minimal_rhf(inputs[1])
  Test.@test scf["Energy"] ≈ -112.3279768245 

  scf = minimal_rhf(inputs[2])
  Test.@test scf["Energy"] ≈ -151.9797610160

  scf = minimal_rhf(inputs[3])
  Test.@test scf["Energy"] ≈ -377.3533493344

  scf = minimal_rhf(inputs[4])
  Test.@test scf["Energy"] ≈ -337.7335100592

  scf = minimal_rhf(inputs[5])
  Test.@test scf["Energy"] ≈ -824.5903952242

  scf = minimal_rhf(inputs[6])
  Test.@test scf["Energy"] ≈ -623.0671382153

  scf = minimal_rhf(inputs[7])
  Test.@test scf["Energy"] ≈ -915.6187742465

  scf = minimal_rhf(inputs[8])
  Test.@test scf["Energy"] ≈ -80.3604114875

  scf = minimal_rhf(inputs[9])
  Test.@test scf["Energy"] ≈ -156.0078302833

  scf = minimal_rhf(inputs[10])
  Test.@test scf["Energy"] ≈ -270.8033880084

  scf = minimal_rhf(inputs[11])
  Test.@test scf["Energy"] ≈ -461.2416993739

  scf = minimal_rhf(inputs[12])
  Test.@test scf["Energy"] ≈ -525.1067034131

  scf = minimal_rhf(inputs[13])
  Test.@test scf["Energy"] ≈ -824.5593071962

  scf = minimal_rhf(inputs[14])
  Test.@test scf["Energy"] ≈ -591.9507402872

  scf = minimal_rhf(inputs[15])
  Test.@test scf["Energy"] ≈ -915.5895688619

  scf = minimal_rhf(inputs[16])
  Test.@test scf["Energy"] ≈ -154.7980797272

  scf = minimal_rhf(inputs[17])
  Test.@test scf["Energy"] ≈ -306.6112223872

  scf = minimal_rhf(inputs[18])
  Test.@test scf["Energy"] ≈ -286.7853098231

  scf = minimal_rhf(inputs[19])
  Test.@test scf["Energy"] ≈ -323.4546492465

  scf = minimal_rhf(inputs[20])
  Test.@test scf["Energy"] ≈ -461.2469197501

  scf = minimal_rhf(inputs[21])
  Test.@test scf["Energy"] ≈ -591.9609978993

  scf = minimal_rhf(inputs[22])
  Test.@test scf["Energy"] ≈ -610.8987019726
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
