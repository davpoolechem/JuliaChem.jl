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
  #== get molecule information from QCArchive ==#
  client = ptl.FractalClient()
  s22_database = client.get_collection("ReactionDataset", "S22")

  mol_index = s22_database.get_index()[1]
  mol_object_full = s22_database.get_molecules(first) #doesn't seem to work
  mol_object = get(mol_object_full.molecule, 0)
  
  #== create input system ==#
  molecule = JSON.parse(mol_object.json())
  display(molecule)
  
  driver = "energy"

  #== create input system ==#
  molecule = JSON.parse(mol.json())

  display(molecule)
  
  driver = "energy"

  model = Dict(
    "method" => "RHF",
    "basis" => "STO-3G"
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
end

JuliaChem.initialize()
script()
JuliaChem.finalize()
