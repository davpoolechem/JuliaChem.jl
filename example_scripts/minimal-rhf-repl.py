#=============================#
#== put needed modules here ==#
#=============================#
import julia
from julia import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
def script(input_file):
  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  basis = JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy"):
    if (model["method"] == "RHF"):
      scf = JuliaChem.JCRHF.run(basis, molecule, keywords)

JuliaChem.initialize()
