#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function minimal_rhf(input_file)
  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  mol, basis = JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  scf = JuliaChem.JCRHF.run(mol, basis, keywords)

  #== reset JuliaChem runtime ==#
  JuliaChem.reset()

  return scf
end
