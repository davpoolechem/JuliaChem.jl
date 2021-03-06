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

  input_time1_t1 = time_ns()/1e9
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)
  input_time1_t2 = time_ns()/1e9
  input_time1 = input_time1_t2 - input_time1_t1

  input_time2_t1 = time_ns()/1e9
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)
  input_time2_t2 = time_ns()/1e9
  input_time2 = input_time2_t2 - input_time2_t1

  input_jit = input_time1 - input_time2
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  basis_time1_t1 = time_ns()/1e9
  basis = JuliaChem.JCBasis.run(molecule, model)
  basis_time1_t2 = time_ns()/1e9
  basis_time1 = basis_time1_t2 - basis_time1_t1

  basis_time2_t1 = time_ns()/1e9
  basis = JuliaChem.JCBasis.run(molecule, model)
  basis_time2_t2 = time_ns()/1e9
  basis_time2 = basis_time2_t2 - basis_time2_t1

  basis_jit = basis_time1 - basis_time2
  basis = JuliaChem.JCBasis.run(molecule, model)

  #== perform scf benchmark ==#
  timeof = Vector{Float64}(undef,0)
  scf_time1 = 0.0
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf_time1_t1 = time_ns()/1e9
      scf = JuliaChem.JCRHF.run(basis, molecule, keywords)
      scf_time1_t2 = time_ns()/1e9
      scf_time1 = scf_time1_t2 - scf_time1_t1

      JuliaChem.reset()

      for index in 1:10
        molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)
        basis = JuliaChem.JCBasis.run(molecule, model)

        scf_timeof_t1 = time_ns()/1e9
        scf = JuliaChem.JCRHF.run(basis, molecule, keywords) #initial run
        scf_timeof_t2 = time_ns()/1e9
        push!(timeof, scf_timeof_t2 - scf_timeof_t1)

        JuliaChem.reset()
      end
    end
  end
  scf_jit = scf_time1 - Statistics.mean(timeof)
 
  #== output relevant information ==# 
  println("Input JIT: ", input_jit)
  println("Basis JIT: ", basis_jit)
  println("SCF JIT: ", scf_jit)
  println("Total JIT: ", input_jit + basis_jit + scf_jit)
  println("")

  #== perform t-test to compare to GAMESS ==#
  p = HypothesisTests.OneSampleTTest(timeof,2.94+0.19)
  println(p)

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()

  #== we are done ==#
  return 0
end
