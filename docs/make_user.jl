using JuliaChem

using Documenter, JuliaChem.JCInput, JuliaChem.JCBasis, JuliaChem.JCRHF 

makedocs(
    authors="David Poole",
    sitename="JuliChem Documentation",
    modules = [JCInput,
               JCBasis,
               JCRHF],
    #format = LaTeX()
)
