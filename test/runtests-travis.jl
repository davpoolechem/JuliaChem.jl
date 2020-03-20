import Pkg; Pkg.add("Test")
import Test

include("../example_scripts/minimal-rhf-repl.jl")

#== select input files ==#
directory = joinpath(@__DIR__, "../example_inputs/S22/")
inputs = readdir(directory)
inputs .= directory .* inputs

#== run S22 tests ==#
Test.@testset "S22 accuracy test" begin
#  scf = minimal_rhf(inputs[1])
#  Test.@test scf[4] ≈ -112.0363811936 

#  scf = minimal_rhf(inputs[2])
#  Test.@test scf[4] ≈ -151.5637222519 

  scf = minimal_rhf(inputs[3])
  Test.@test scf[4] ≈ -376.3871095400

#  scf = minimal_rhf(inputs[4])
#  Test.@test scf[4] ≈ -336.8722288570 

#  scf = minimal_rhf(inputs[5])
#  Test.@test scf[4] ≈ -822.5278362138 

#  scf = minimal_rhf(inputs[6])
#  Test.@test scf[4] ≈ -621.5325650057

#  scf = minimal_rhf(inputs[7])
#  Test.@test scf[4] ≈ -913.3428352684

#  scf = minimal_rhf(inputs[8])
#  Test.@test scf[4] ≈ -80.1635874796

#  scf = minimal_rhf(inputs[9])
#  Test.@test scf[4] ≈ -155.6363922320

#  scf = minimal_rhf(inputs[10])
#  Test.@test scf[4] ≈ -270.1512256025

#  scf = minimal_rhf(inputs[11])
#  Test.@test scf[4] ≈ -460.1352685104

#  scf = minimal_rhf(inputs[12])
#  Test.@test scf[4] ≈ -523.8050189484

#  scf = minimal_rhf(inputs[13])
#  Test.@test scf[4] ≈ -822.4956611365

#  scf = minimal_rhf(inputs[14])
#  Test.@test scf[4] ≈ -590.5292651211

#  scf = minimal_rhf(inputs[15])
#  Test.@test scf[4] ≈ -913.3126739534

#  scf = minimal_rhf(inputs[16])
#  Test.@test scf[4] ≈ -154.4360760563

#  scf = minimal_rhf(inputs[17])
#  Test.@test scf[4] ≈ -305.8480434104

#  scf = minimal_rhf(inputs[18])
#  Test.@test scf[4] ≈ -286.0846164433

#  scf = minimal_rhf(inputs[19])
#  Test.@test scf[4] ≈ -322.6781928490

#  scf = minimal_rhf(inputs[20])
#  Test.@test scf[4] ≈ -460.1400279398

#  scf = minimal_rhf(inputs[21])
#  Test.@test scf[4] ≈ -590.5383556631

#  scf = minimal_rhf(inputs[22])
#  Test.@test scf[4] ≈ -609.3928162070
end

#== finalize JuliaChem ==#
JuliaChem.finalize()       
