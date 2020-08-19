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
  Test.@test scf["Energy"] ≈ -112.4047194144

  scf = minimal_rhf(inputs[2])
  Test.@test scf["Energy"] ≈ -152.0671593526

  scf = minimal_rhf(inputs[3])
  Test.@test scf["Energy"] ≈ -377.5889420312

  scf = minimal_rhf(inputs[4])
  Test.@test scf["Energy"] ≈ -337.9245601137

  scf = minimal_rhf(inputs[5])
  Test.@test scf["Energy"] ≈ -825.0480089998

  scf = minimal_rhf(inputs[6])
  Test.@test scf["Energy"] ≈ -623.3855095909

  scf = minimal_rhf(inputs[7])
  Test.@test scf["Energy"] ≈ -916.1413738989

  scf = minimal_rhf(inputs[8])
  Test.@test scf["Energy"] ≈ -80.4073036256

  scf = minimal_rhf(inputs[9])
  Test.@test scf["Energy"] ≈ -156.0878533332

  scf = minimal_rhf(inputs[10])
  Test.@test scf["Energy"] ≈ -270.9302743515

  scf = minimal_rhf(inputs[11])
  Test.@test scf["Energy"] ≈ -461.4497600413

  scf = minimal_rhf(inputs[12])
  Test.@test scf["Energy"] ≈ -525.4098675881

  scf = minimal_rhf(inputs[13])
  Test.@test scf["Energy"] ≈ -825.0271971388

  scf = minimal_rhf(inputs[14])
  Test.@test scf["Energy"] ≈ -592.2263255700

  scf = minimal_rhf(inputs[15])
  Test.@test scf["Energy"] ≈ -916.1246231957

  scf = minimal_rhf(inputs[16])
  Test.@test scf["Energy"] ≈ -154.8750579586

  scf = minimal_rhf(inputs[17])
  Test.@test scf["Energy"] ≈ -306.7601557517

  scf = minimal_rhf(inputs[18])
  Test.@test scf["Energy"] ≈ -286.9286565278

  scf = minimal_rhf(inputs[19])
  Test.@test scf["Energy"] ≈ -323.6144109921

  scf = minimal_rhf(inputs[20])
  Test.@test scf["Energy"] ≈ -461.4537798340

  scf = minimal_rhf(inputs[21])
  Test.@test scf["Energy"] ≈ -592.2351527581

  scf = minimal_rhf(inputs[22])
  Test.@test scf["Energy"] ≈ -611.1928406081
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
