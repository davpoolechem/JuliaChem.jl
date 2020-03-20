#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import InteractiveUtils

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file)
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  basis = JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf = JuliaChem.JCRHF.run(basis, molecule, keywords)
    end
  end

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()
end

script(ARGS[1])
