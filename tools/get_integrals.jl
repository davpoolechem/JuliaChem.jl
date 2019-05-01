function get_integrals(input::String, nshell::Int64)
    f_int::IOStream = open(input)
        integral_strings::Array{String,1} = readlines(f_int)
    close(f_int)

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

    integrals::Array{String,1} = []
    push!(integrals,"{")
    push!(integrals,"   \"Input\":\"Two-Electron\",")
    for int::Int64 in 1:length(integrals_temp)
        ish::Int64 = integrals_temp[int][1]
        jsh::Int64 = integrals_temp[int][2]
        ksh::Int64 = integrals_temp[int][3]
        lsh::Int64 = integrals_temp[int][4]
        integral::Float64 = integrals_temp[int][5]
        if (int == 1)
            push!(integrals,"   \"tei\":[ [$ish, $jsh, $ksh, $lsh, $integral],")
        elseif (int == length(integrals_temp))
            push!(integrals,"           [$ish, $jsh, $ksh, $lsh, $integral] ]")
        else
            push!(integrals,"           [$ish, $jsh, $ksh, $lsh, $integral],")
        end
    end
    push!(integrals,"}")

    f_int_w::IOStream = open("integrals.log","w")
        for iline in integrals
            write(f_int_w,"$iline\n")
        end
    close(f_int_w)

#=
	integrals_temp_b::Array{Array{Array{Array{Float64,1}}}} = []

	for ish::Int64 in 1:nshell
	  push!(integrals_temp_b,[])

      for jsh::Int64 in 1:nshell
	    push!(integrals_temp_b[ish],[])

        for ksh::Int64 in 1:nshell
	      push!(integrals_temp_b[ish][jsh],[])

          for lsh::Int64 in 1:nshell
	        push!(integrals_temp_b[ish][jsh][ksh], zero(Float64))
		  end
		end
	  end
    end

	for int in 1:length(integrals_temp)
	  ish = integrals_temp[int][1]
	  jsh = integrals_temp[int][2]
	  ksh = integrals_temp[int][3]
	  lsh = integrals_temp[int][4]
      integrals_temp_b[ish][jsh][ksh][lsh] = integrals_temp[int][5]
    end

	integrals::Array{String,1} = []
    for ish::Int64 in 1:nshell, jsh::Int64 in 1:nshell
      for ksh::Int64 in 1:nshell, lsh::Int64 in 1:nshell
	    integral::Float64 = integrals_temp_b[ish][jsh][ksh][lsh]
	    push!(integrals,"[$ish, $jsh, $ksh, $lsh, $integral]")
	  end
    end

	#display(integrals)
	=#
end

get_integrals("sto3g-water.int",7)
