#== extract geometry and symbols ==# 
function get_geometry_and_symbols(input_file_string::String)
  #== read in input file ==#
  input_file = [] 
  open(input_file_string) do file
    input_file = readlines(file)
  end

  #== find DATA line ==#
  data_line = 0
  for input_line in 1:length(input_file)
    if occursin("\$DATA", input_file[input_line])
      data_line = input_line
      break
    end
  end

  #== extract geometry from input ==# 
  starting = data_line + 3 
  ending = length(input_file) - 1

  geometry = input_file[starting:ending]
  symbols::Array{String,1} = []

  #== write output string array ==#
  output_file_array::Array{String,1} = []
  push!(output_file_array,"    \"geometry\": [")
  
  for geometry_line in geometry
    #coords_match = eachmatch(r"([-]?[\d]{1,}\.[\d]{1,}\s{0,})",geometry_line)
    coords_match = eachmatch(r"([-]?[\d]{1,}\.[\d]{1,})",geometry_line)
    coords = collect(coords_match)
    coords_x = coords[2].match 
    coords_y = coords[3].match
    coords_z = coords[4].match

    if geometry_line == geometry[end]
      push!(output_file_array, "      $coords_x, $coords_y, $coords_z")
    else
      push!(output_file_array, "      $coords_x, $coords_y, $coords_z,")
    end
 
    symbol_match = match(r"[A-Z]",geometry_line)
    push!(symbols, symbol_match.match)
  end

  push!(output_file_array,"    ],") 
  push!(output_file_array,"    \"symbols\": $symbols,") 

  #== write to output file ==#
  open("geometry.json", "w") do file
    for line in output_file_array
      write(file, "$line\n")
    end
  end
end

function create_input(input_file_string::String) 
  #== read in geometry ==# 
  geometry::Array{String,1} = []
  open("geometry.json") do file
    geometry = readlines(file)
  end

  #== generate output file array ==#
  output_file_array::Array{String,1} = [] 
  push!(output_file_array,"{")
  push!(output_file_array,"  \"molecule\": {")
  push!(output_file_array, geometry...)
  push!(output_file_array,"    \"molecular_charge\":0")
  push!(output_file_array,"  },")
  push!(output_file_array,"  \"driver\": \"energy\",")
  push!(output_file_array,"  \"model\": {")
  push!(output_file_array,"    \"method\": \"RHF\",")
  push!(output_file_array,"    \"basis\": \"PCSeg-0\"")
  push!(output_file_array,"  },")
  push!(output_file_array,"  \"keywords\": {")
  push!(output_file_array,"    \"scf\":{")
  push!(output_file_array,"      \"niter\":100,")
  push!(output_file_array,"      \"ndiis\":8,")
  push!(output_file_array,"      \"dele\":1E-10,")
  push!(output_file_array,"      \"rmsd\":1E-8,")
  push!(output_file_array,"      \"prec\":\"Float64\",")
  push!(output_file_array,"      \"direct\":true,")
  push!(output_file_array,"      \"debug\":false,")
  push!(output_file_array,"      \"load\":\"static\"")
  push!(output_file_array,"    }")
  push!(output_file_array,"  }")
  push!(output_file_array,"}")
  
  #== write to output file ==#
  output_file = input_file_string*".json" 
  open(output_file, "w") do file
    for line::String in output_file_array
      write(file, "$line\n")
    end
  end
end

#== create_all S22 inputs ==#
#function get_all_S22(directory)
#  for input_number in 1:22
input_file_inp_string = ARGS[1]
get_geometry_and_symbols(input_file_inp_string*".inp") 
create_input(input_file_inp_string)
#  end
#end

#get_all_S22(ARGS[1])


