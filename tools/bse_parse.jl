module BSEParse

using PyCall
using HDF5
using JSON

bse = pyimport("basis_set_exchange")
const am_to_shell_mapping = Dict(
    0 => "S", 
    1 => "P",
    2 => "D", 
    3 => "F", 
    4 => "G", 
    5 => "H", 
    6 => "I" 
  )

#===========================================================#
#= perform parsing operation on single atom/basis set pair =#
#===========================================================#
function parse_individual(atom::Dict{String,Any}, atomid::String, basis::String, 
  bsed::HDF5File)
  
  for (shell_num, shell) in enumerate(atom["electron_shells"])
    ang_mom = shell["angular_momentum"]
    shell_type_string = length(ang_mom) > 1 ? "L" : am_to_shell_mapping[ang_mom[1]]  
    
    exponents::Array{Float64,1} = parse.(Float64,shell["exponents"])
   
    coeff::Array{Float64} = begin
      if shell_type_string == "L"
        s_coeff = parse.(Float64,shell["coefficients"][1])
        p_coeff = parse.(Float64,shell["coefficients"][2])
        vcat(s_coeff,p_coeff)
      else
        parse.(Float64,shell["coefficients"][1])
      end
    end

    #println("Writing to $atomid/$basis/$shell_num...")
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Shell Type", shell_type_string)
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Exponents", exponents)
    h5write("bsed.h5",
      "$atomid/$basis/$shell_num/Coefficients", coeff)
    
    #=
    println("Exponents:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$shell_num/Exponents"))
    println("")
    println("Contraction Coefficients:")
    display(h5read("bsed.h5",
      "$atomid/$basis/$shell_num/Coefficients"))
    println("")
    =#
  end
end

#===========================================#
#= parse all selected atom-basis set pairs =#
#===========================================#
function parse_all()
    #== open HDF5 file for writing ==#
    hdf5name::String = "bsed.h5"
    h5open(hdf5name, "w") do bsed

        #== parse "sto" basis family ==#
        atoms::Array{String,1} = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe",
            "Cs"
            ] #H-Cs
        #=
        atoms::Array{String,1} = [
            "1", "2",
            "3", "4", "5", "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15", "16", "17", "18",
            "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
            "30", "31", "32", "33", "34", "35", "36",
            "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48",
            "49", "50", "51", "52", "53", "54",
            "55"
            ] #H-Cs
        =#
        basis_sets::Array{String,1} = ["STO-2G", "STO-3G", "STO-4G", "STO-5G",
            "STO-6G" ] #STO family
        
        for basis::String in basis_sets
            println("Basis: $basis")
            bs_dict = bse.get_basis(basis,fmt="json", header=false)
            bs_dict_json = JSON.parse(bs_dict)

            for atom in bs_dict_json["elements"] 
                parse_individual(atom.second, atoms[parse(Int64,atom.first)], 
                  basis, bsed)
            end          
        end

        #== parse "pople" basis family ==#
        basis_sets = ["3-21G"          , "4-31G"             , "5-21G"        , "6-21G"       , 
                      "6-31G"          , "6-31G*"            , "6-31G**"      , 
                      "6-31+G"         , "6-31+G*"           , "6-31+G**"     , "6-31+G*-J"   ,
                      "6-31++G"        , "6-31++G*"          , "6-31++G**"    , "6-31++G**-J" ,
                      "6-31G(2df,p)"   , "6-31G(3df,3pd)"    , "6-31G(d,p)"   , "6-31G-J"     ,
                      "6-311G"         , "6-311G*"           , "6-311G**"     , 
                      "6-311+G"        , "6-311+G*"          , "6-311+G**"    , "6-311+G*-J"  ,
                      "6-311++G"       , "6-311++G*"         , "6-311++G**"   , "6-311++G**-J",
                      "6-311G(d,p)"    , "6-311G(2df,2pd)"   , "6-311+G(2d,p)", 
                      "6-311++G(2d,2p)", "6-311++G(3df,3pd)" ] #pople family

        for basis::String in basis_sets
            println("Basis: $basis")
            bs_dict = bse.get_basis(basis,fmt="json", header=false)
            bs_dict_json = JSON.parse(bs_dict)

            for atom in bs_dict_json["elements"] 
                parse_individual(atom.second, atoms[parse(Int64,atom.first)], 
                  basis, bsed)
            end
        end

        #== parse polarization-consistent basis family ==#
        basis_sets = ["PCSeg-0"] #polarization-consistent family
        for basis::String in basis_sets
            println("Basis: $basis")
            bs_dict = bse.get_basis(basis,fmt="json", header=false)
            bs_dict_json = JSON.parse(bs_dict)

            for atom in bs_dict_json["elements"] 
                parse_individual(atom.second, atoms[parse(Int64,atom.first)], 
                  basis, bsed)
            end
        end

        #== parse correlation-consistent basis family ==#
        #=
        basis_sets = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"] 
        for basis::String in basis_sets
            println("Basis: $basis")
            bs_dict = bse.get_basis(basis,fmt="json", header=false)
            bs_dict_json = JSON.parse(bs_dict)

            for atom in bs_dict_json["elements"] 
              parse_individual(atom.second, atoms[parse(Int64,atom.first)], 
                basis, bsed)
            end
        end
        =#
    end
end

end

BSEParse.parse_all()
