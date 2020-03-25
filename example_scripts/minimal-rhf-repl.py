#=============================#
#== put needed modules here ==#
#=============================#
import sys

import julia
from julia import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
#== initialize JuliaChem ==#
JuliaChem.initialize()

#== read in input file ==#
molecule, driver, model, keywords = JuliaChem.JCInput.run(sys.argv[0])

#== generate basis set ==#
mol, basis = JuliaChem.JCBasis.run(molecule, model)

#== perform scf calculation ==#
if (driver == "energy"):
  if (model["method"] == "RHF"):
    scf = JuliaChem.JCRHF.run(mol, basis, keywords)

JuliaChem.reset()
