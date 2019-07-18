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

get_hamiltonian_integrals("example_inputs/pc0-water.log")
get_overlap_integrals("example_inputs/pc0-water.log")
get_two_electron_integrals("example_inputs/pc0-water.log")
