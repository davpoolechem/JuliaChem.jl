#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file)
  #== initialize JuliaChem runtime ==#
  JuliaChem.initialize()

  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file)

    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model)
  
    #== print molecular information ==#
    JuliaChem.JCMolecule.run(mol)

    #== perform scf calculation ==#
    #if (driver == "energy")
    #  if (model["method"] == "RHF")
    #    scf = JuliaChem.JCRHF.run(mol, basis, keywords)
    #  end
    #end
  catch e
    bt = catch_backtrace()
    msg = sprint(showerror, e, bt)
    println(msg)

    JuliaChem.finalize()
    exit()
  end
  
  #== finalize JuliaChem runtime ==#
  JuliaChem.finalize()
end

script(ARGS[1])
