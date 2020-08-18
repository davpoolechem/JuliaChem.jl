#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function minimal_rhf(input_file)
  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="verbose")       
    
    #== generate basis set ==#
    mol, basis, jeri_engine = JuliaChem.JCBasis.run(molecule, model; 
      output="none")          

    JuliaChem.JCMolecule.run(mol)

    #== perform scf calculation ==#
    rhf_energy = JuliaChem.JCRHF.run(mol, basis, jeri_engine, keywords["scf"]; 
      output="verbose") 
 
    #== perform gradient ==#
    rhf_gradient = JuliaChem.JCGrad.run(mol, basis; output="verbose")

    #== reset JuliaChem runtime ==#
    JuliaChem.reset()
    return rhf_energy, rhf_gradient
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end
