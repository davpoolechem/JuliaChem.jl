#=============================#
#== put needed modules here ==#
#=============================#
import qcportal as ptl

import julia
from julia import JuliaChem
from julia import JSON

#================================#
#== JuliaChem execution script ==#
#================================#
def script():
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== get molecule information from QCArchive ==# 
  client = ptl.FractalClient()
  s22_database = client.get_collection("ReactionDataset", "S22")
  
  mol_index = s22_database.get_index()[0]
  mol_value = s22_data.get_molecules(mol_index)

  print(mol_value.to_json())

  #== create input system ==#
  molecule = JSON.parse(mol_value.to_json())

  driver = "energy"

  model = { 
    "method": "RHF",
    "basis": "6-31G"
  }

  keywords = { 
    "scf": { 
      "niter": 50,
      "ndiis": 3,
      "dele": 1E-8,
      "rmsd": 1E-6,
      "prec": "Float64",
      "direct": False,
      "debug": False
    }
  }

  #== generate basis set ==#
  basis = JuliaChem.JCBasis.run(molecule, model)
  
  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()

script()
