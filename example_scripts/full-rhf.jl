#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function full_rhf(input_file)
  try
    #== initialize JuliaChem ==#
    JuliaChem.initialize()

    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="verbose")       
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; 
      output="verbose")          

    #== calculation driver ==# 
    if driver == "energy"
      if model["method"] == "RHF"
        #== perform scf calculation ==#
        rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; 
          output="verbose") 
     
        #== compute molecular properties such as dipole moment ==#
        JuliaChem.JCProperties.run(mol, basis, rhf_energy, keywords, 
          output="verbose")
      else
        throw("Exception: Methods other than RHF are not supported yet!")
      end  
    else
      throw("Exception: Only energy calculations are currently supported!")
    end
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end

full_rhf(ARGS[1])
