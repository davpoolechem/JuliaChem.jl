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
  
#== extract nuclear energy ==# 
function get_nuclear_energy(input_file_string::String)
  input_file::Array{String,1} = []
  open(input_file_string) do file
    input_file = readlines(file)
  end
  
  #== find nuclear energy line ==#
  nuclear_energy_line = "" 
  for input_line in input_file
    if occursin("NUCLEAR ENERGY =", input_line)
      nuclear_energy_line = input_line
      break
    end
  end

  #== get nuclear energy from line ==#
  nuclear_energy_match = match(r"=(.*)", nuclear_energy_line)

  nuclear_energy = parse(Float64, nuclear_energy_match.match[2:end])
  return nuclear_energy 
end


#=
#== extract hamiltonian integrals from file ==#
function get_hamiltonian_integrals(input_file_string::String)
    #== read in input file ==#
    input_file::Array{String,1} = []
    open(input_file_string) do file
        input_file = readlines(file)
    end

    #== extract integrals from input file ==#
    starting::Int64 = 1
    ending::Int64 = 1

    for line_num in 1:length(input_file)
        if occursin("START HAMILTONIAN INTEGRALS",input_file[line_num])
            starting = line_num
        elseif occursin("END HAMILTONIAN INTEGRALS",input_file[line_num])
            ending = line_num
        end
    end

    integrals::Array{String,1} = input_file[starting+1:ending-1]

    #== write output string array ==#
    output_file_array::Array{String,1} = []
    push!(output_file_array,"\"hcore\":[ "*integrals[1]*",")
    for integral_num in 2:length(integrals)
        integral::String = integrals[integral_num]
        if integral_num == length(integrals)
            push!(output_file_array,"$integral"*" ]")
        else
            push!(output_file_array,"$integral"*",")
        end
    end

    #== write to output file ==#
    open("hcore.json", "w") do file
        for line::String in output_file_array
            write(file, "$line\n")
        end
    end
end

#== extract overlap integrals from file ==#
function get_overlap_integrals(input_file_string::String)
    #== read in input file ==#
    input_file::Array{String,1} = []
    open(input_file_string) do file
        input_file = readlines(file)
    end

    #== extract integrals from input file ==#
    starting::Int64 = 1
    ending::Int64 = 1

    for line_num in 1:length(input_file)
        if occursin("START OVERLAP INTEGRALS",input_file[line_num])
            starting = line_num
        elseif occursin("END OVERLAP INTEGRALS",input_file[line_num])
            ending = line_num
        end
    end

    integrals::Array{String,1} = input_file[starting+1:ending-1]

    #== write output string array ==#
    output_file_array::Array{String,1} = []
    push!(output_file_array,"\"ovr\":[ "*integrals[1]*",")
    for integral_num in 2:length(integrals)
        integral::String = integrals[integral_num]
        if integral_num == length(integrals)
            push!(output_file_array,"$integral"*" ],")
        else
            push!(output_file_array,"$integral"*",")
        end
    end

    #== write to output file ==#
    open("overlap.json", "w") do file
        for line::String in output_file_array
            write(file, "$line\n")
        end
    end
end

#== extract overlap integrals from file ==#
function get_overlap_integrals(input_file_string::String)
    #== read in input file ==#
    input_file::Array{String,1} = []
    open(input_file_string) do file
        input_file = readlines(file)
    end

    #== extract integrals from input file ==#
    starting::Int64 = 1
    ending::Int64 = 1

    for line_num in 1:length(input_file)
        if occursin("START OVERLAP INTEGRALS",input_file[line_num])
            starting = line_num
        elseif occursin("END OVERLAP INTEGRALS",input_file[line_num])
            ending = line_num
        end
    end

    integrals::Array{String,1} = input_file[starting+1:ending-1]

    #== write output string array ==#
    output_file_array::Array{String,1} = []
    push!(output_file_array,"\"ovr\":[ "*integrals[1]*",")
    for integral_num in 2:length(integrals)
        integral::String = integrals[integral_num]
        if integral_num == length(integrals)
            push!(output_file_array,"$integral"*" ],")
        else
            push!(output_file_array,"$integral"*",")
        end
    end

    #== write to output file ==#
    open("overlap.json", "w") do file
        for line::String in output_file_array
            write(file, "$line\n")
        end
    end
end

#== extract two-electron integrals from file ==#
function get_two_electron_integrals(input_file_string::String)
    #== read in input file ==#
    input_file::Array{String,1} = []
    open(input_file_string) do file
        input_file = readlines(file)
    end

    #== extract integrals from input file ==#
    starting::Int64 = 1
    ending::Int64 = 1

    for line_num in 1:length(input_file)
        if occursin("START TWO-ELECTRON INTEGRALS",input_file[line_num])
            starting = line_num
        elseif occursin("END TWO-ELECTRON INTEGRALS",input_file[line_num])
            ending = line_num
        end
    end

    integrals::Array{String,1} = input_file[starting+1:ending-1]

    #== write output string array ==#
    output_file_array::Array{String,1} = []
    push!(output_file_array,"\"tei\":[ "*integrals[1]*",")
    for integral_num in 2:length(integrals)
        integral::String = integrals[integral_num]
        if integral_num == length(integrals)
            push!(output_file_array,"$integral"*" ],")
        else
            push!(output_file_array,"$integral"*",")
        end
    end

    #== write to output file ==#
    open("tei.json", "w") do file
        for line::String in output_file_array
            write(file, "$line\n")
        end
    end
end
=#
function create_input(input_file_string::String)
  #== read in geometry ==# 
  geometry::Array{String,1} = []
  open("geometry.json") do file
    geometry = readlines(file)
  end

  #== read in overlap integrals ==# 
  overlap::Array{String,1} = []
  open("overlap.json") do file
    overlap = readlines(file)
  end

  #== read in hcore integrals ==# 
  hcore::Array{String,1} = []
  open("hcore.json") do file
    hcore = readlines(file)
  end

  #== generate output file array ==#
  output_file_array::Array{String,1} = [] 
  push!(output_file_array,"{")
  push!(output_file_array,"  \"molecule\": {")
  push!(output_file_array, geometry...)
  push!(output_file_array,"    \"molecular_charge\":0,")
  
  nuclear_energy = get_nuclear_energy(input_file_string)
  push!(output_file_array,"    \"enuc\":$nuclear_energy ")

  #push!(output_file_array, overlap...)
  #push!(output_file_array, hcore...)

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
  push!(output_file_array,"      \"load\":static")
  push!(output_file_array,"    }")
  push!(output_file_array,"  }")
  push!(output_file_array,"}")
  #== write to output file ==#
  output_file_match = match(r"^(.*?)\.",input_file_string)
  output_file = output_file_match.match*"json" 
  open(output_file, "w") do file
    for line::String in output_file_array
      write(file, "$line\n")
    end
  end
end

#== create_all S22 inputs ==#
function get_all_S22()
  for input_number in 1:22
    input_number_string = input_number < 10 ? "0"*"$input_number" : "$input_number"
    input_file_log_string = "S22/"*input_number_string*"_MP2.log" 
    input_file_inp_string = "S22/"*input_number_string*"_MP2.inp" 
    println("Working on input file "*input_file_log_string)
   
    get_geometry_and_symbols(input_file_inp_string) 
    get_hamiltonian_integrals(input_file_log_string)
    get_overlap_integrals(input_file_log_string)
    create_input(input_file_log_string)
  end
end

get_all_S22()


