using Documenter, Input, Molecule, RHF, Properties

makedocs(
    sitename="JuliChem Documentation",
    modules = [Input,
               Molecule,
               RHF,
               Properties],
)
