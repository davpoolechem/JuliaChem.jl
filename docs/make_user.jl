using Documenter, Input, Molecule, RHF, Properties
using DocumenterLaTeX

makedocs(
    authors="David Poole",
    sitename="JuliChem Documentation",
    modules = [Input,
               Molecule,
               RHF,
               Properties],
    #format = LaTeX()
)
