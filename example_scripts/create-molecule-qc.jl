#=============================#
#== put needed modules here ==#
#=============================#
import PyCall
ptl = PyCall.pyimport("qcportal")

import JuliaChem

import JSON

#================================#
#== JuliaChem execution script ==#
#================================#
function script()
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  #== get molecule information from QCArchive ==# 
  client = ptl.FractalClient()
  mol = client.query_molecules(6)[1]

  #== create input system ==#
  molecule = JSON.parse(mol.json())

  driver = "energy"

  model = Dict(
    "method" => "RHF",
    "basis" => "6-31G(d,p)"
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
  basis = JuliaChem.JCBasis.run(molecule, model)

  display(basis[1])
  
  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()
end

script()
