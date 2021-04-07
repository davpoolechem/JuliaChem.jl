#=============================#
#== put needed modules here ==#
#=============================#
import sys

import julia
from julia import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
def rhf_properties_xyz(input_file):
  try:
    #== read in xyz input file ==#
    mol_charge = 0
    molecule = JuliaChem.JCInput.xyz_to_molecule(input_file, mol_charge) 

    #== define other input parameters ==#
    driver = "energy"      
    
    model = { 
      "method": "RHF",
      "basis": "6-31G*"
    }

    keywords = {
      "scf": {},
      "prop": {
        "formation": True,
        "mo energies": True,
        "mulliken": True,
        "multipole": "dipole"
      } 
    }
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model, 
      output=0)          

    #JuliaChem.JCMolecule.run(mol)

    #== perform scf calculation ==#
    rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"],
      output=1)

    #== perform property calculations ==#
    rhf_properties = JuliaChem.JCRHF.Properties.run(mol, basis, rhf_energy,
        keywords["prop"], output=2)  

    return rhf_energy, rhf_properties
  except Exception as e:                                                            
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    print(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
