#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

using JuliaChem.JCInput
using JuliaChem.JCBasis
using JuliaChem.JCRHF

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file::String)
  #== read in input file ==#
  molecule, driver, model, keywords = JCInput.run(input_file)

  #== generate basis set ==#
  basis = JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf = JCRHF.run(basis, molecule, keywords)
      write("output.json",JSON.json(scf[5]))
    end
  end
end

JuliaChem.initialize()
