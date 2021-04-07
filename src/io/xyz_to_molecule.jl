using JSON

function xyz_to_geometry(xyzfile)
  input_file = []
  #println(joinpath(@__DIR__, xyzfile))
  open(xyzfile) do file
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
    append!(coords_vector, parse.(Float64,atom_coords_match[1:3]))

    symbol_match = match(r"[A-Za-z]{1,2}",geometry_line)
    push!(symbols, symbol_match.match)
  end
  #coords = reshape(coords_vector, (3,length(geometry)))
  #coords = transpose(coords)
  coords = coords_vector

  #@assert natoms == size(coords)[1]
  #@assert natoms == length(symbols)

  #display(coords); println()
  #display(symbols); println()

  return coords, symbols
end

function xyz_to_molecule(input, charge = 0)
  geom_json = "{\n"
  
  #== write in molecule geometry ==#
  coords, symbols = xyz_to_geometry(input)
  geom_json *= "    \"geometry\" : [\n" 
  for icoord in 1:length(coords)
    geom_json *= "      " 
      
    value = coords[icoord]
    #println(value)
    if icoord == length(coords)
      geom_json *= "$value\n"
    elseif icoord%3 == 0 
      geom_json *= "$value,\n"
    else
      geom_json *= "$value,"
    end
  end
  geom_json *= "    ],\n" 

  #== write molecule symbols ==#
  geom_json *= "    \"symbols\" : [\n" 
  for iatom in 1:length(symbols)
    geom_json *= "      " 
    symbol = symbols[iatom]
    if iatom == length(symbols)
      geom_json *= "\"$symbol\"\n"
    elseif iatom%5 == 0 
      geom_json *= "\"$symbol\",\n"
    else
      geom_json *= "\"$symbol\","
    end
  end
  geom_json *= "    ],\n" 

  #== write molecule charge ==#
  geom_json *= "    \"molecular_charge\": $charge\n" 
  geom_json *= "}\n"

  molecule = JSON.parse(geom_json)
  return molecule
end 
export xyz_to_molecule
