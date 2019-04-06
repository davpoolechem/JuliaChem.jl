using Documenter, Input, Molecule, RHF, Properties, InputStructs
using DocumenterLaTeX

makedocs(
    authors="David Poole",
    sitename="JuliChem Documentation",
    modules = [Input,
               Molecule,
               RHF,
               Properties,
               InputStructs],
    #format = LaTeX()
)
