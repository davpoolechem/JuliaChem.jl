import Test

include("../tools/travis-rhf.jl")

#== select input files ==#
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

display(inputs)

#== initialize JuliaChem ==#
JuliaChem.initialize()


#== run S22 calculations ==#
s22_test_results = Dict([]) 
for imol in [ 2, 3 ] 
  s22_test_results[imol] = travis_rhf(inputs[imol])
end

#== check energies ==#
Test.@testset "S22 Energy" begin
  #Test.@test s22_test_results[1][:Energy]["Energy"] ≈ -112.4047194144
  Test.@test s22_test_results[2][:Energy]["Energy"] ≈ -152.0671593526
  Test.@test s22_test_results[3][:Energy]["Energy"] ≈ -377.5889420312
  #Test.@test s22_test_results[4][:Energy]["Energy"] ≈ -337.9245601137
  #Test.@test s22_test_results[5][:Energy]["Energy"] ≈ -825.0480089998
  #Test.@test s22_test_results[6][:Energy]["Energy"] ≈ -623.3855095909
  #Test.@test s22_test_results[7][:Energy]["Energy"] ≈ -916.1413738989
  #Test.@test s22_test_results[8][:Energy]["Energy"] ≈ -80.4073036256
  #Test.@test s22_test_results[9][:Energy]["Energy"] ≈ -156.0878533332
  #Test.@test s22_test_results[10][:Energy]["Energy"] ≈ -270.9302743515
  #Test.@test s22_test_results[11][:Energy]["Energy"] ≈ -461.4497600413
  #Test.@test s22_test_results[12][:Energy]["Energy"] ≈ -525.4098675881
  #Test.@test s22_test_results[13][:Energy]["Energy"] ≈ -825.0271971388
  #Test.@test s22_test_results[14][:Energy]["Energy"] ≈ -592.2263255700
  #Test.@test s22_test_results[15][:Energy]["Energy"] ≈ -916.1246231957
  #Test.@test s22_test_results[16][:Energy]["Energy"] ≈ -154.8750579586
  #Test.@test s22_test_results[17][:Energy]["Energy"] ≈ -306.7601557517
  #Test.@test s22_test_results[18][:Energy]["Energy"] ≈ -286.9286565278
  #Test.@test s22_test_results[19][:Energy]["Energy"] ≈ -323.6144109921
  #Test.@test s22_test_results[20][:Energy]["Energy"] ≈ -461.4537798340
  #Test.@test s22_test_results[21][:Energy]["Energy"] ≈ -592.2351527581
  #Test.@test s22_test_results[22][:Energy]["Energy"] ≈ -611.1928406081
end

#== check dipole moments ==#
Test.@testset "S22 Dipoles" begin
  #Test.@test abs(s22_test_results[1][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  Test.@test s22_test_results[2][:Properties]["Dipole"][:moment] ≈ 2.696653
  Test.@test abs(s22_test_results[3][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test abs(s22_test_results[4][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test abs(s22_test_results[5][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test s22_test_results[6][:Properties]["Dipole"][:moment] ≈ 3.239972 
  #Test.@test s22_test_results[7][:Properties]["Dipole"][:moment] ≈ 1.990013 
  #Test.@test abs(s22_test_results[8][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test abs(s22_test_results[9][:properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test s22_test_results[10][:Properties]["Dipole"][:moment] ≈ 0.280546 
  #Test.@test abs(s22_test_results[11][:Properties]["Dipole"][:moment]) <= 1.0E-6 #approximately 0.0        
  #Test.@test s22_test_results[12][:Properties]["Dipole"][:moment] ≈ 0.280546 
  #Test.@test s22_test_results[13][:Properties]["Dipole"][:moment] ≈ 5.604885 
  #Test.@test s22_test_results[14][:Properties]["Dipole"][:moment] ≈ 1.860770 
  #Test.@test s22_test_results[15][:Properties]["Dipole"][:moment] ≈ 3.478761 
  #Test.@test s22_test_results[16][:Properties]["Dipole"][:moment] ≈ 0.493345 
  #Test.@test s22_test_results[17][:Properties]["Dipole"][:moment] ≈ 2.422111 
  #Test.@test s22_test_results[18][:Properties]["Dipole"][:moment] ≈ 1.770430 
  #Test.@test s22_test_results[19][:Properties]["Dipole"][:moment] ≈ 4.142949 
  #Test.@test s22_test_results[20][:Properties]["Dipole"][:moment] ≈ 0.651396 
  #Test.@test s22_test_results[21][:Properties]["Dipole"][:moment] ≈ 3.194881 
  #Test.@test s22_test_results[22][:Properties]["Dipole"][:moment] ≈ 3.486793 
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
