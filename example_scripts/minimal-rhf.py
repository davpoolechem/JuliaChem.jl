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
  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file,
    output="none")

  #== generate basis set ==#
  mol, basis = JuliaChem.JCBasis.run(molecule, model, output="none")

  #== perform scf calculation ==#
  scf = JuliaChem.JCRHF.run(mol, basis, keywords["scf"], output="minimal")
  
  #== reset JuliaChem runtime ==#
  JuliaChem.reset()
  return scf

JuliaChem.initialize()
scf = minimal_rhf(sys.argv[1])
JuliaChem.finalize()
