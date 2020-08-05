import Test

include("../example_scripts/minimal-rhf-repl.jl")

#== select input files ==#:w
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

display(inputs)

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
  scf = minimal_rhf(inputs[1])
  Test.@test scf["Energy"] ≈ -112.3731088312

  scf = minimal_rhf(inputs[2])
  Test.@test scf["Energy"] ≈ -152.0298289870

  scf = minimal_rhf(inputs[3])
  Test.@test scf["Energy"] ≈ -377.5431167329

  scf = minimal_rhf(inputs[4])
  Test.@test scf["Energy"] ≈ -337.8788667317

  scf = minimal_rhf(inputs[5])
  Test.@test scf["Energy"] ≈ -824.9593611954

  scf = minimal_rhf(inputs[6])
  Test.@test scf["Energy"] ≈ -623.3152616759

  scf = minimal_rhf(inputs[7])
  Test.@test scf["Energy"] ≈ -916.0396657781

  scf = minimal_rhf(inputs[8])
  Test.@test scf["Energy"] ≈ -80.3897573362

  scf = minimal_rhf(inputs[9])
  Test.@test scf["Energy"] ≈ -156.0622806362

  scf = minimal_rhf(inputs[10])
  Test.@test scf["Energy"] ≈ -270.8967526555

  scf = minimal_rhf(inputs[11])
  Test.@test scf["Energy"] ≈ -461.3993110203

  scf = minimal_rhf(inputs[12])
  Test.@test scf["Energy"] ≈ -525.3578242397

  scf = minimal_rhf(inputs[13])
  Test.@test scf["Energy"] ≈ -824.9362675225

  scf = minimal_rhf(inputs[14])
  Test.@test scf["Energy"] ≈ -592.1602243851

  scf = minimal_rhf(inputs[15])
  Test.@test scf["Energy"] ≈ -916.0206299206

  scf = minimal_rhf(inputs[16])
  Test.@test scf["Energy"] ≈ -154.8498650045

  scf = minimal_rhf(inputs[17])
  Test.@test scf["Energy"] ≈ -306.7164812315

  scf = minimal_rhf(inputs[18])
  Test.@test scf["Energy"] ≈ -286.8877785284

  scf = minimal_rhf(inputs[19])
  Test.@test scf["Energy"] ≈ -323.5790904725

  scf = minimal_rhf(inputs[20])
  Test.@test scf["Energy"] ≈ -461.4043292962

  scf = minimal_rhf(inputs[21])
  Test.@test scf["Energy"] ≈ -592.1704932173

  scf = minimal_rhf(inputs[22])
  Test.@test scf["Energy"] ≈ -611.1204856578
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
