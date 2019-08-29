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
function script()
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== create list of geometries to analyze ==#
  geometries = map(x -> [x,0.0,0.0,0.0,0.0,0.0], LinRange(1.0,2.0,11))

  for geometry in geometries
    #== create input system ==#
    molecule = Dict(
      "geometry" => geometry,
      "symbols" => ["H", "H"],
      "charge" => 0,
      "enuc" => 0.5291772492,
    )

    driver = "energy"
 
    model = Dict(
      "method" => "RHF",
      "basis" => "PCSeg-0"
    )

    keywords = Dict(
      "scf" => Dict(
        "niter" => 50,
        "ndiis" => 3,
        "dele" => 1E-8,
        "rmsd" => 1E-6,
        "prec" => "Float64",
        "direct" => false,
        "debug" => false
      )
    )

    #== generate basis set ==#
    basis = JCBasis.run(molecule, model)
  end

  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()
end

script()
