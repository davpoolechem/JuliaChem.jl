#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import Statistics
#import HypothesisTests 
import MPI
using BenchmarkTools

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file)
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  try
    #== read in input file ==#
    input_time1_t1 = time_ns()/1e9
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file; 
      output=0)
    input_time1_t2 = time_ns()/1e9
    input_time1 = input_time1_t2 - input_time1_t1 

    input_time2_t1 = time_ns()/1e9
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file; 
      output=0)
    input_time2_t2 = time_ns()/1e9
    input_time2 = input_time2_t2 - input_time2_t1
  
    input_jit = input_time1 - input_time2
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file; 
      output=0)

    #== generate basis set ==#
    basis_time1_t1 = time_ns()/1e9
    basis = JuliaChem.JCBasis.run(molecule, model; output=0)
    basis_time1_t2 = time_ns()/1e9
    basis_time1 = basis_time1_t2 - basis_time1_t1 

    basis_time2_t1 = time_ns()/1e9
    basis = JuliaChem.JCBasis.run(molecule, model; output=0)
    basis_time2_t2 = time_ns()/1e9
    basis_time2 = basis_time2_t2 - basis_time2_t1
  
    basis_jit = basis_time1 - basis_time2
    mol, basis = JuliaChem.JCBasis.run(molecule, model; output=0)

    #== perform scf benchmark ==#
    timeof = Vector{Float64}(undef,0)
    scf_time1 = 0.0
    if (driver == "energy")
      if (model["method"] == "RHF")
        #scf_time1_t1 = time_ns()/1e9
        #scf = JuliaChem.JCRHF.run(mol, basis, keywords["scf"]; 
        #  output=2) 
        #scf_time1_t2 = time_ns()/1e9
        #scf_time1 = scf_time1_t2 - scf_time1_t1 
      
        molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file; 
          output=0)
        mol, basis = JuliaChem.JCBasis.run(molecule, model; output=0)

        scf_timeof_t1 = time_ns()/1e9
        #scf = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; 
        #  output=2)
        scf = BenchmarkTools.@benchmark begin
          JuliaChem.JCRHF.Energy.run($mol, $basis, $(keywords["scf"]); 
            output=2) #initial run
        end
        display(scf)
        scf_timeof_t2 = time_ns()/1e9
        push!(timeof, scf_timeof_t2 - scf_timeof_t1) 
      end
    end
    scf_jit = scf_time1 - Statistics.mean(timeof)
  
    MPI.Barrier(MPI.COMM_WORLD)
    if (MPI.Comm_rank(MPI.COMM_WORLD) == 0) display(timeof) end
    
    #if (MPI.Comm_rank(MPI.COMM_WORLD) == 0)  
      #== output relevant information ==# 
    #  println("Input JIT: ", input_jit)
    #  println("Basis JIT: ", basis_jit)
    #  println("SCF JIT: ", scf_jit)
    #  println("Total JIT: ", input_jit + basis_jit + scf_jit)
    #  println("")
  
      #== perform t-test to compare to GAMESS ==#
    #  p = HypothesisTests.OneSampleTTest(timeof,2.94+0.19)
    #  println(p)
    #end
    MPI.Barrier(MPI.COMM_WORLD)
    
    #== finalize JuliaChem runtime ==#
    JuliaChem.finalize()
  
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end    
end

script(ARGS[1])
