using HDF5
using JSON

function get_integrals(input::String, nbf::Int64)
    input_file = open(input)
        input_string = readlines(input_file)
    close(input_file)

    input_info = Dict([])
    i = 1
    while (i <= length(input_string))
        if (input_string[i] != "{")
            i += 1
        else
            j = i
            json_subsection = ""
            input_name = ""
            while (input_string[j] != "}")
                json_subsection *= input_string[j]
                if (occursin(r"\\\"Input\\\"\:(.*)",input_string[j]))
                    input_name = match(r"\\\"Input\\\"\:(.*)",input_string[j])[1]
                    input_name = input_name[2:end-2]
                end
                j += 1
            end
            json_subsection *= "}"
            json_parse = JSON.parse(json_subsection)
            merge!(input_info,Dict([(input_name,json_parse)]))

            i = j
        end
    end

    nbf2::Int64 = nbf*(nbf+1)/2
    h5open("tei.h5", "w") do file
        for i::Int64 in 1:nbf2
            array::Array{Array{Float64,1},1} = input_info["Two-Electron-$i"]["tei"]
            matrix::Array{Float64,2} = fill(0.0,(length(array),5))
            for irow::Int64 in 1:length(array), icol::Int64 in 1:5
                matrix[irow,icol] = array[irow][icol]
            end
            write(file, "tei-$i", matrix)  # alternatively, say "@write file A"
        end
    end

    c = h5open("tei.h5", "r") do file
        test = read(file, "tei-5")
        display(test[1,5])
    end
end

function test(input::String, nbf::Int64)
    get_integrals(input,nbf)
end

test("sto3g-water-gms.json",7)
