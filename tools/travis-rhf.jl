#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function travis_rhf(input_file)
  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="verbose")       
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; output="verbose")          

    #== compute molecular inform ation ==#
    JuliaChem.JCMolecule.run(mol)
 
    #== perform scf calculation ==#
    rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; 
      output="verbose") 

    #== compute molecular properties ==# 
    rhf_properties = JuliaChem.JCProperties.run(mol, basis, rhf_energy, keywords["prop"],
      output="verbose")  

    return (Energy = rhf_energy, Properties = rhf_properties) 
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end
