function xyz_to_qc(ARGS)
  #== read in xyz file ==#
  xyz = Vector{String}([])
  open(ARGS[1]) do file
    input_file = readlines(file)
    xyz = input_file[3:end]
  end
  
  #== write file ==#
  filename = splitext(ARGS[1])[1]
  output = open("$filename.json", "w") do file
    #== print initial stuff ==#
    write(file, "{\n")
    write(file, "  \"molecule\": {\n")

    #== print geometry ==# 
    write(file, "    \"geometry\": [ \n") 
    for atom in xyz
      atom_coords_match = eachmatch(r"([-]?[\d]{1,}\.[\d]{1,})",atom)
      atom_coords = collect(atom_coords_match)
      atom_x = atom_coords[1].match
      atom_y = atom_coords[2].match
      atom_z = atom_coords[3].match 
  
      if atom == xyz[end] 
        write(file, "                  $atom_x, $atom_y, $atom_z \n") 
      else
        write(file, "                  $atom_x, $atom_y, $atom_z, \n") 
      end
    end
    write(file, "                ],\n") 

    #== print symbols ==#
    write(file, "    \"symbols\": [ \n")
    for atom in xyz
      atom_symbol_match=match(r"[A-Za-z]{1,2}",atom).match
      if atom == xyz[end] 
        write(file, "                  \"$atom_symbol_match\" \n") 
      else
        write(file, "                  \"$atom_symbol_match\", \n") 
      end
    end
    write(file, "                ],\n") 

    #== write rest of molecule subsection ==#
    write(file, "    \"molecular_charge\":0\n")
    write(file, "  },\n")

    #== write driver subsection ==#
    write(file, "  \"driver\": \"energy\",\n")

    #== write model subsection ==#
    write(file, "  \"model\": {\n")
    write(file, "    \"method\": \"RHF\",\n")
    write(file, "    \"basis\": \"6-31+G(d)\n")
    write(file, "  },\n")
    
    #== write keywords subsection ==#
    write(file, "  \"keywords\": {\n")
    write(file, "    \"scf\":{\n")
    write(file, "      \"niter\":100,\n")
    write(file, "      \"ndiis\":8,\n")
    write(file, "      \"dele\":1E-8,\n")
    write(file, "      \"rmsd\":1E-5,\n")
    write(file, "      \"prec\":\"Float64\",\n")
    write(file, "      \"direct\":true,\n")
    write(file, "      \"debug\":false,\n")
    write(file, "      \"load\":\"static\",\n")
    write(file, "      \"fdiff\":true\n")
    write(file, "    }\n")
    write(file, "  }\n")
    write(file, "}\n")
  end
end

xyz_to_qc(ARGS)
