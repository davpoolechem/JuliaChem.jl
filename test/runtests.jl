import Test

include("../tools/travis-rhf.jl")

#== select input files ==#:w
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

display(inputs)

#== initialize JuliaChem ==#
JuliaChem.initialize()

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
  rhf_energy, properties = travis_rhf(inputs[1])
  Test.@test rhf_energy["Energy"] ≈ -112.4047194144
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[2])
  Test.@test rhf_energy["Energy"] ≈ -152.0671593526
  Test.@test properties["Dipole"][:moment] ≈ 2.696653

  rhf_energy, properties = travis_rhf(inputs[3])
  Test.@test rhf_energy["Energy"] ≈ -377.5889420312
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[4])
  Test.@test rhf_energy["Energy"] ≈ -337.9245601137
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[5])
  Test.@test rhf_energy["Energy"] ≈ -825.0480089998
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[6])
  Test.@test rhf_energy["Energy"] ≈ -623.3855095909
  Test.@test properties["Dipole"][:moment] ≈ 3.239972 

  rhf_energy, properties = travis_rhf(inputs[7])
  Test.@test rhf_energy["Energy"] ≈ -916.1413738989
  Test.@test properties["Dipole"][:moment] ≈ 1.990013 

  rhf_energy, properties = travis_rhf(inputs[8])
  Test.@test rhf_energy["Energy"] ≈ -80.4073036256
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[9])
  Test.@test rhf_energy["Energy"] ≈ -156.0878533332
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[10])
  Test.@test rhf_energy["Energy"] ≈ -270.9302743515
  Test.@test properties["Dipole"][:moment] ≈ 0.280546 

  rhf_energy, properties = travis_rhf(inputs[11])
  Test.@test rhf_energy["Energy"] ≈ -461.4497600413
  Test.@test abs(properties["Dipole"][:moment]) <= 1.0E-10 #approximately 0.0        

  rhf_energy, properties = travis_rhf(inputs[12])
  Test.@test rhf_energy["Energy"] ≈ -525.4098675881
  Test.@test properties["Dipole"][:moment] ≈ 0.280546 

  rhf_energy, properties = travis_rhf(inputs[13])
  Test.@test rhf_energy["Energy"] ≈ -825.0271971388
  Test.@test properties["Dipole"][:moment] ≈ 5.604885 

  rhf_energy, properties = travis_rhf(inputs[14])
  Test.@test rhf_energy["Energy"] ≈ -592.2263255700
  Test.@test properties["Dipole"][:moment] ≈ 1.860770 

  rhf_energy, properties = travis_rhf(inputs[15])
  Test.@test rhf_energy["Energy"] ≈ -916.1246231957
  Test.@test properties["Dipole"][:moment] ≈ 3.478761 

  rhf_energy, properties = travis_rhf(inputs[16])
  Test.@test rhf_energy["Energy"] ≈ -154.8750579586
  Test.@test properties["Dipole"][:moment] ≈ 0.493345 

  rhf_energy, properties = travis_rhf(inputs[17])
  Test.@test rhf_energy["Energy"] ≈ -306.7601557517
  Test.@test properties["Dipole"][:moment] ≈ 2.422111 

  rhf_energy, properties = travis_rhf(inputs[18])
  Test.@test rhf_energy["Energy"] ≈ -286.9286565278
  Test.@test properties["Dipole"][:moment] ≈ 1.770430 

  rhf_energy, properties = travis_rhf(inputs[19])
  Test.@test rhf_energy["Energy"] ≈ -323.6144109921
  Test.@test properties["Dipole"][:moment] ≈ 4.142949 

  rhf_energy, properties = travis_rhf(inputs[20])
  Test.@test rhf_energy["Energy"] ≈ -461.4537798340
  Test.@test properties["Dipole"][:moment] ≈ 0.651396 

  rhf_energy, properties = travis_rhf(inputs[21])
  Test.@test rhf_energy["Energy"] ≈ -592.2351527581
  Test.@test properties["Dipole"][:moment] ≈ 3.194881 

  rhf_energy, properties = travis_rhf(inputs[22])
  Test.@test rhf_energy["Energy"] ≈ -611.1928406081
  Test.@test properties["Dipole"][:moment] ≈ 3.486793 
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
