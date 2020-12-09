#=============================#
#== put needed modules here ==#
#=============================#
import sys

import julia
from julia import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
def minimal_rhf(input_file):
  try:  
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file, 
      output="none")

    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model, output="none")

    #== perform scf calculation ==#
    rhf_energy = {}
    if "scf" in keywords:
      rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"],
        output="verbose")
    else:
      rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis,
        output="verbose")

    return rhf_energy 
  except Exception as e:                                                        
    bt = Base.catch_backtrace()                                                 
    msg = Base.sprint(Base.showerror, e, bt)
    print(msg)                                                                  
                                                                                
    JuliaChem.finalize()                                                        
    exit()    
