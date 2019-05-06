function get_tei_integrals(input::String, nbf::Int64)
    #--read in input file--#
    f_int::IOStream = open(input)
        integral_strings::Array{String,1} = readlines(f_int)
    close(f_int)

    #--process read-in file to create integral list--#
    integrals_temp::Array{Array{Real,1}} = []
	for iline::String in integral_strings
      if (!occursin("II",iline) && !occursin("TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS",iline))
	    integral_line_1::Array{Real,1} = []
		push!(integral_line_1,parse(Int64, iline[1:4]))
		push!(integral_line_1,parse(Int64, iline[5:8]))
		push!(integral_line_1,parse(Int64, iline[9:12]))
		push!(integral_line_1,parse(Int64, iline[13:16]))
		push!(integral_line_1,parse(Float64, iline[27:38]))
	   	push!(integrals_temp,integral_line_1)

        integral_line_2::Array{Real,1} = []
	    push!(integral_line_2,parse(Int64, iline[39:43]))
		push!(integral_line_2,parse(Int64, iline[44:47]))
		push!(integral_line_2,parse(Int64, iline[48:51]))
		push!(integral_line_2,parse(Int64, iline[52:55]))
		push!(integral_line_2,parse(Float64, iline[66:77]))
	   	push!(integrals_temp,integral_line_2)
      end
    end

	sort!(integrals_temp)

    #--format integrals for json format--#
    integrals::Array{String,1} = []
	nbf2::Int64 = nbf*(nbf+1)/2

	for bf in 1:nbf2
		ibf::Int64 = ceil(((-1+sqrt(1+8*bf))/2))
		jbf::Int64 = bf%ibf + 1
		nadd::Int64 = 0

		push!(integrals,"{")
		push!(integrals,"   \"Input\":\"Two-Electron-$bf\",")

	    for int::Int64 in 1:length(integrals_temp)
	        ish::Int64 = integrals_temp[int][1]
	        jsh::Int64 = integrals_temp[int][2]
	        ksh::Int64 = integrals_temp[int][3]
	        lsh::Int64 = integrals_temp[int][4]
	        integral::Float64 = integrals_temp[int][5]
			if (ish == ibf && jsh == jbf)
				if (nadd == 0 && ish == ksh && jsh == lsh)
					push!(integrals,"   \"tei\":[ [$ish, $jsh, $ksh, $lsh, $integral] ]")
		        elseif (nadd == 0)
		            push!(integrals,"   \"tei\":[ [$ish, $jsh, $ksh, $lsh, $integral],")
		        elseif (ish == ksh && jsh == lsh)
		            push!(integrals,"           [$ish, $jsh, $ksh, $lsh, $integral] ]")
		        else
		            push!(integrals,"           [$ish, $jsh, $ksh, $lsh, $integral],")
		        end
				nadd += 1
			end
	    end
	    push!(integrals,"}")
	end

    #--write json block to output file--#
    f_tei_w::IOStream = open("tei_integrals.log","w")
        for iline in integrals
            write(f_tei_w,"$iline\n")
        end
    close(f_tei_w)
end

#=
function get_oei_integrals(input::String, nbf::Int64)
    #--read in input file--#
    f_int::IOStream = open(input)
        integral_strings::Array{String,1} = readlines(f_int)
    close(f_int)

    #--determine starting locations of OEI blocks--#
    ovr_start = 0
    ovr_end = 0

    hcore_start = 0
    hcore_end = 0

    kei_start = 0
    kei_end = length(integral_strings)

    for i in 1:length(integral_strings)
        if (occursin("OVERLAP MATRIX",integral_strings[i]))
            ovr_start = i
        elseif (occursin("BARE NUCLEUS HAMILTONIAN INTEGRALS",integral_strings[i]))
            hcore_start = i
            ovr_end = i-1
        elseif (occursin("KINETIC ENERGY INTEGRALS",integral_strings[i]))
            kei_start = i
            hcore_end = i-1
        end
    end

    ovr = integral_strings[ovr_start:ovr_end]
    hcore = integral_strings[hcore_start:hcore_end]
    kei = integral_strings[kei_start:kei_end]

    #display(ovr[3+(1-1)*(nbf+3)][17:25])

    ovr::Array{Real,1} = []
    iblock::Int64 = 0
    while (3+iblock*(nbf+3) < length(ovr))
        columns::String = ovr[3+iblock*(nbf+3)]
        icols::Array{Int64,1} = [ parse(Int64,columns[17:min(25,length(columns))]) ]

        if (length(columns) > 25)
            push!(icols,parse(Int64,columns[28:min(36,length(columns))]))
        end
        if (length(columns)) > 37
            push!(icols,parse(Int64,columns[39:min(47,length(columns))]))
        end
        if (length(columns)) > 48
            push!(icols,parse(Int64,columns[50:min(58,length(columns))]))
        end
        if (length(columns)) > 59
            push!(icols,parse(Int64,columns[61:min(69,length(columns))]))
        end

        first = 5+iblock*(nbf+3)
        last = min(4+nbf+iblock*(nbf+3),length(ovr))
        display(ovr[first:last])



        iblock += 1
    end
end
=#

get_tei_integrals("sto3g-water.tei",7)
