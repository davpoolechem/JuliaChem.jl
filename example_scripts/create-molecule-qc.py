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
  mol = client.query_molecules(1234)[0]

  #== create input system ==#
  molecule = JSON.parse(mol.json())

  driver = "energy"

  model = { 
    "method": "RHF",
    "basis": "6-31G(d,p)"
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
