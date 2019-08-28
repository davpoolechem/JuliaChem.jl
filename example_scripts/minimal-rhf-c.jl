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
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== read in input file ==#
  input_file::String = ARGS[1]
  @time molecule, driver, model, keywords = JCInput.run(input_file)

  #== generate basis set ==#
  @time basis = JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      @time scf = JCRHF.run(basis, molecule, keywords)
    end
  end

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()

  #== we are done ==#
  return 0
end
