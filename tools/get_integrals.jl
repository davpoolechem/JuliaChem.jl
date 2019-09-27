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
            push!(output_file_array,"$integral"*" ],")
        else
            push!(output_file_array,"$integral"*",")
        end
    end

    #== write to output file ==#
    open("tools/hcore.json", "w") do file
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
    open("tools/overlap.json", "w") do file
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
    open("tools/overlap.json", "w") do file
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
    open("tools/tei.json", "w") do file
        for line::String in output_file_array
            write(file, "$line\n")
        end
    end
end

function create_input()
  #== read in overlap integrals ==# 
  overlap::Array{String,1} = []
  open("tools/overlap.json") do file
    overlap = readlines(file)
  end

  #== read in hcore integrals ==# 
  hcore::Array{String,1} = []
  open("tools/hcore.json") do file
    hcore = readlines(file)
  end

  #== read in two-electron integrals ==# 
  tei::Array{String,1} = []
  open("tools/tei.json") do file
    tei = readlines(file)
  end

  #== generate output file array ==#
  output_file_array::Array{String,1} = [] 
  push!(output_file_array,"{")
  push!(output_file_array,"  \"molecule\": {")
  push!(output_file_array,"    \"geometry\": [")
  push!(output_file_array,"      0.000000000000,  -0.143225816552,   0.000000000000,")
  push!(output_file_array,"      1.638036840407,   1.136548822547,  -0.000000000000,")
  push!(output_file_array,"     -1.638036840407,   1.136548822547,  -0.000000000000
    ],")
  push!(output_file_array,"    \"symbols\": [\"O\", \"H\", \"H\"],")
  push!(output_file_array,"    \"molecular_charge\":0,")
  push!(output_file_array,"    \"enuc\":4.2346705892, ")

  push!(output_file_array, overlap...)
  push!(output_file_array, hcore...)
  push!(output_file_array, tei...)

  push!(output_file_array,"  },")
  push!(output_file_array,"  \"driver\": \"energy\",")
  push!(output_file_array,"  \"model\": {")
  push!(output_file_array,"    \"method\": \"RHF\",")
  push!(output_file_array,"    \"basis\": \"6-31G\"")
  push!(output_file_array,"  },")
  push!(output_file_array,"  \"keywords\": {")
  push!(output_file_array,"    \"scf\":{")
  push!(output_file_array,"      \"niter\":100,")
  push!(output_file_array,"      \"ndiis\":8,")
  push!(output_file_array,"      \"dele\":1E-8,")
  push!(output_file_array,"      \"rmsd\":1E-6,")
  push!(output_file_array,"      \"prec\":\"Float64\",")
  push!(output_file_array,"      \"direct\":false,")
  push!(output_file_array,"      \"debug\":false")
  push!(output_file_array,"    }")
  push!(output_file_array,"  }")
  push!(output_file_array,"}")
  #== write to output file ==#
  open("tools/input.json", "w") do file
    for line::String in output_file_array
      write(file, "$line\n")
    end
  end
end

get_hamiltonian_integrals(ARGS[1])
get_overlap_integrals(ARGS[1])
#get_two_electron_integrals(ARGS[1])
create_input()
