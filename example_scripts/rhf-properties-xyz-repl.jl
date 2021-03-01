#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function rhf_properties_xyz(input_file)
  try
    #== read in xyz input file ==#
    mol_charge = 0
    molecule = JuliaChem.JCInput.xyz_to_molecule(input_file, mol_charge) 

    #== define other input parameters ==#
    driver = "energy"      
    
    model = Dict(
      "method" => "RHF",
      "basis" => "6-31G*"
    )

    keywords = Dict(
      "scf" => Dict(),
      "prop" => Dict(
        "formation" => true,
        "mo energies" => true,
        "mulliken" => true,
        "multipole" => "dipole"
      ) 
    )
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; 
      output=0)          

    #JuliaChem.JCMolecule.run(mol)

    #== perform scf calculation ==#
    rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"];
      output=1)

    #== perform property calculations ==#
    rhf_properties = JuliaChem.JCRHF.Properties.run(mol, basis, rhf_energy,
        keywords["prop"]; output=2)  

    return rhf_energy, rhf_properties
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end
