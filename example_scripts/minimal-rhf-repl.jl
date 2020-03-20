#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import InteractiveUtils

#================================#
#== JuliaChem execution script ==#
#================================#
function minimal_rhf(input_file)
  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

  #== generate basis set ==#
  InteractiveUtils.@code_warntype JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf = JuliaChem.JCRHF.run(basis, molecule, keywords)
    end
  end

  #== reset JuliaChem runtime ==#
  JuliaChem.reset()
  return scf
end

JuliaChem.initialize()
