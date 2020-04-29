using HDF5

function test()
  try
    h5open(joinpath(@__DIR__,"../records/bsed.h5"), "r") do file
      #== set up some calculation parameters ==#
      driver = "energy"

      keywords = Dict(
        "scf" => Dict(
          "niter" => 50,
          "ndiis" => 8,
          "dele" => 1E-8,
          "rmsd" => 1E-6,
          "prec" => "Float64",
          "direct" => true,
          "debug" => false,
          "load" => "static"
        )
      )

      #== loop over basis sets in bsed ==# 
      basis_sets = read(file)
      for basis_set in keys(basis_sets)
        #== set up more calculation parameters ==#
        model = Dict(
          "method" => "RHF",
          "basis" => basis_set 
        )

        #== loop over atoms described by basis set ==# 
        for atom in keys(basis_sets[basis_set])
          #== create molecule parameters ==#
          molecule = Dict("geometry" => [ 0.0, 0.0, 0.0 ],
                          "symbols" => [ atom ],
                          "molecular_charge" => 0 
                         )       
          open(joinpath(@__DIR__,"sad_inputs/$basis_set-$atom.inp"), 
            "w") do gamess_input
            
            write(gamess_input, "$basis_set + $atom") 
          end 
          #== compute density matrix ==#   
          #mol, basis = JuliaChem.JCBasis.run(molecule, model; output="none")

          #scf = JuliaChem.JCRHF.run(mol, basis, keywords["scf"]; 
          #  output="none")
  
          #display(scf["Density"])
        end
      end
    end
  catch e
    bt = catch_backtrace()
    msg = sprint(showerror, e, bt)
    println(msg)

    JuliaChem.finalize()
    exit()
  end
end

test()
