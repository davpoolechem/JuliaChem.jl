#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
#== initialize JuliaChem ==#
JuliaChem.initialize()

#== read in input file ==#
molecule, driver, model, keywords = JuliaChem.JCInput.run(ARGS[1])

#== generate basis set ==#
mol, basis = JuliaChem.JCBasis.run(molecule, model)

#== perform scf calculation ==#
if (driver == "energy")
  if (model["method"] == "RHF")
    scf = JuliaChem.JCRHF.run(mol, basis, keywords)
  end
end

#== reset JuliaChem runtime ==#
JuliaChem.reset()

