using Documenter, JCInput, JCMolecule, JCRHF, JCProperties, JCStructs
using DocumenterLaTeX

makedocs(
    authors="David Poole",
    sitename="JuliChem Documentation",
    modules = [JCInput,
               JCMolecule,
               JCRHF,
               JCProperties,
               JCStructs],
    #format = LaTeX()
)
