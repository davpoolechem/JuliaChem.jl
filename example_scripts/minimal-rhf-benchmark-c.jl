#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import Statistics
import HypothesisTests

#================================#
#== JuliaChem execution script ==#
#================================#
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== read in input file ==#
  input_file = ARGS[1]
  input_time1 = @elapsed JuliaChem.JCInput.run(input_file)
  input_time2 = @elapsed JuliaChem.JCInput.run(input_file)
  input_jit = input_time1 - input_time2
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  basis_time1 = @elapsed JuliaChem.JCBasis.run(molecule, model)
  basis_time2 = @elapsed JuliaChem.JCBasis.run(molecule, model)
  basis_jit = basis_time1 - basis_time2
  basis = JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  timeof = Vector{Float64}(undef,0) 
  scf_time1 = 0.0
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf_time1 = @elapsed JuliaChem.JCRHF.run(basis, molecule, keywords)
      JuliaChem.reset()
      
      for index in 1:100
        molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)
        basis = JuliaChem.JCBasis.run(molecule, model)

        push!(timeof, @elapsed JuliaChem.JCRHF.run(basis, molecule, keywords))
        JuliaChem.reset()
      end
    end
  end
  scf_jit = scf_time1 - Statistics.median(timeof)
 
  #== output relevant information ==# 
  println("Median: ", Statistics.median(timeof))
  println("Standard Deviation: ", Statistics.std(timeof))
  println("")
  println("Input JIT: ", input_jit)
  println("Basis JIT: ", basis_jit)
  println("SCF JIT: ", scf_jit)
  println("SCF JIT: ", input_jit + basis_jit + scf_jit)
  println("")

  #== perform t-test to compare to GAMESS ==#
  p = HypothesisTests.OneSampleTTest(timeof,0.39+0.17)
  println(p)

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()

  #== we are done ==#
  return 0
end
