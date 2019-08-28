import julia

julia.install()

from julia import JuliaChem

from julia import JSON
from julia import MPI

def script(input_file):
    MPI.Init()

    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

    basis = JuliaChem.JCBasis.run(molecule, model)

    if (driver == "energy"):
      if (model["method"] == "RHF"):
        scf = JuliaChem.JCRHF.run(basis, molecule, keywords)

    MPI.Finalize()

script("example_inputs/631Gdp-H2O.json")
