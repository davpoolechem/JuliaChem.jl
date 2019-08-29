#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

using JuliaChem.JCInput
using JuliaChem.JCBasis
using JuliaChem.JCRHF

using BenchmarkTools

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file::String)
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== read in input file ==#
  molecule, driver, model, keywords = JCInput.run(input_file)

  #== generate basis set ==#
  basis = JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf = JCRHF.run(basis, molecule, keywords)
    end
  end

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()
end

script(ARGS[1])
