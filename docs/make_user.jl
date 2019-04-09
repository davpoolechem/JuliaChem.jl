using Documenter, JCInput, JCMolecule, JCRHF, JCProperties, InputStructs
using DocumenterLaTeX

makedocs(
    authors="David Poole",
    sitename="JuliChem Documentation",
    modules = [JCInput,
               JCMolecule,
               JCRHF,
               JCProperties,
               InputStructs],
    #format = LaTeX()
)
