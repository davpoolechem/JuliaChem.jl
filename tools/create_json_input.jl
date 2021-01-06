function xyz_to_geometry(xyzfile)
  input_file = []
  open(joinpath(@__DIR__, xyzfile)) do file
    input_file = readlines(file)
  end

  natoms = parse(Int, input_file[1])
  #println(natoms)
  geometry = input_file[3:end]
  #display(geometry); println()

  coords_vector = Vector{Float64}([])
  symbols = Vector{String}([])
  for geometry_line in geometry
    coords_match = eachmatch(r"([-]?[\d]{1,}\.[\d]{1,})",geometry_line)
    atom_coords = collect(coords_match)
    atom_coords_match = map(x -> x.match, atom_coords) 
    append!(coords_vector, parse.(Float64,atom_coords_match))
    
    symbol_match = match(r"[A-Z]{1,2}",geometry_line)
    push!(symbols, symbol_match.match)
  end
  #coords = reshape(coords_vector, (3,length(geometry)))
  #coords = transpose(coords)
  coords = coords_vector

  #@assert natoms == size(coords)[1]
  #@assert natoms == length(symbols)

  #display(coords); println()
  #display(symbols); println()

  charge = 0

  return coords, symbols, charge
end

function create_input(input)
  #== create json file if it does not exist ==#
  filename = input[1:(end-4)] * ".json"

  #== write input json file ==#
  output = open(joinpath(@__DIR__, filename), write=true) do file
    write(file, "{\n")
    #== write in molecule information ==#
    write(file, "  \"molecule\": {\n")
    coords, symbols, charge = xyz_to_geometry(input)
    #== write fragment geometry ==#
    write(file, "    \"geometry\" : [\n") 
    for icoord in 1:length(coords)
      write(file, "       ") 
        
      value = coords[icoord]
      #println(value)
      if icoord == length(coords)
        write(file, "$value\n")
      elseif icoord%3 == 0 
        write(file, "$value,\n")
      else
        write(file, "$value,")
      end
    end
    write(file, "    ],\n") 
 
    #== write fragment symbols ==#
    write(file, "    \"symbols\" : [\n") 
    for iatom in 1:length(symbols)
      write(file, "      ") 
      symbol = symbols[iatom]
      if iatom == length(symbols)
        write(file, "\"$symbol\"\n")
      elseif iatom%5 == 0 
        write(file, "\"$symbol\",\n")
      else
        write(file, "\"$symbol\",")
      end
    end
    write(file, "    ],\n") 
    #== write fragment charge ==#
    write(file, "    \"molecular_charge\": $charge\n") 
    write(file, "  },\n")
 
    #== write calculation driver and model information ==# 
    write(file,"  \"driver\": \"energy\",\n")
    write(file,"  \"model\": {\n")
    write(file,"    \"method\": \"RHF\",\n")
    write(file,"    \"basis\": \"6-31G*\"\n")
    write(file,"  },\n")
    
    #== write keywords ==# 
    write(file,"  \"keywords\": {\n")
    write(file,"    \"scf\": {\n")
    write(file,"      \"dele\":1.0,\n")
    write(file,"      \"rmsd\":2E-5\n")
    write(file,"    }\n")
    write(file,"  }\n")
    write(file,"}")
  end 
end
