function get_tei_integrals(input::String, nbf::Int64)
    #--read in input file--#
    f_int::IOStream = open(input)
        integral_strings::Array{String,1} = readlines(f_int)
    close(f_int)

    integrals::Array{String,1} = []
<<<<<<< HEAD
	nbf2::Int64 = nbf*(nbf+1)/2

	for ibf in 1:nbf, jbf in 1:ibf
		bf::Int64 = (ibf-1)*ibf/2 + jbf
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
    f_tei_w::IOStream = open("tei_integrals.json","w")
=======

    #--write tei block--#
    push!(integrals,"{")
    push!(integrals,"   \"Input\":\"Two-Electron\",")
    for iline in integral_strings
        display(iline)
        i::Int64 = parse(Int64,iline[1:5])
        j::Int64 = parse(Int64,iline[20:25])
        k::Int64 = parse(Int64,iline[43:48])
        l::Int64 = parse(Int64,iline[60:66])

        integral::Float64 = parse(Float64,iline[67:end])

        if (i == 1 && j == 1 && k == 1 && l == 1 && length(integral_strings) == 1)
            push!(integrals,"   \"tei\":[ [$i, $j, $k, $l, $integral] ]")
        elseif (i == 1 && j == 1 && k == 1 && l == 1)
            push!(integrals,"   \"tei\":[ [$i, $j, $k, $l, $integral],")
        elseif (i == nbf && j == nbf && k == nbf && l == nbf)
            push!(integrals,"           [$i, $j, $k, $l, $integral] ]")
        else
            push!(integrals,"           [$i, $j, $k, $l, $integral],")
        end
    end
    push!(integrals,"}")
    push!(integrals,"")

    #--write json block to output file--#
    f_tei_w::IOStream = open("tools/tei_integrals.log","w")
>>>>>>> development
        for iline in integrals
            write(f_tei_w,"$iline\n")
        end
    close(f_tei_w)
end

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
            ovr_start = i+1
        elseif (occursin("BARE NUCLEUS HAMILTONIAN INTEGRALS",integral_strings[i]))
            hcore_start = i+1
            ovr_end = i-1
        elseif (occursin("KINETIC ENERGY INTEGRALS",integral_strings[i]))
            kei_start = i+1
            hcore_end = i-1
        end
    end

    ovr_string = integral_strings[ovr_start:ovr_end]
    hcore_string = integral_strings[hcore_start:hcore_end]
    kei_string = integral_strings[kei_start:kei_end]
    integrals::Array{String,1} = []

    #--write oei block--#
    push!(integrals,"{")
    push!(integrals,"   \"Input\":\"Overlap\",")
    for iline in ovr_string
        display(iline)
        i::Int64 = parse(Int64,iline[1:5])
        j::Int64 = parse(Int64,iline[20:23])
        integral::Float64 = parse(Float64,iline[26:end])

        if (i == 1 && j == 1 && length(ovr_string) == 1)
            push!(integrals,"   \"ovr\":[ [$i, $j, $integral] ]")
        elseif (i == 1 && j == 1)
            push!(integrals,"   \"ovr\":[ [$i, $j, $integral],")
        elseif (i == nbf && j == nbf)
            push!(integrals,"           [$i, $j, $integral] ]")
        else
            push!(integrals,"           [$i, $j, $integral],")
        end
    end
    push!(integrals,"}")
    push!(integrals,"")

    #--write hcore block--#
    push!(integrals,"{")
    push!(integrals,"   \"Input\":\"One-Electron Hamiltonian\",")
    for iline in hcore_string
        display(iline)
        i::Int64 = parse(Int64,iline[1:5])
        j::Int64 = parse(Int64,iline[20:23])
        integral::Float64 = parse(Float64,iline[26:end])

        if (i == 1 && j == 1 && length(hcore_string) == 1)
            push!(integrals,"   \"hcore\":[ [$i, $j, $integral] ]")
        elseif (i == 1 && j == 1)
            push!(integrals,"   \"hcore\":[ [$i, $j, $integral],")
        elseif (i == nbf && j == nbf)
            push!(integrals,"           [$i, $j, $integral] ]")
        else
            push!(integrals,"           [$i, $j, $integral],")
        end
    end
    push!(integrals,"}")
    push!(integrals,"")

    #--write json block to output file--#
    f_oei_w::IOStream = open("tools/oei_integrals.log","w")
        for iline in integrals
            write(f_oei_w,"$iline\n")
        end
    close(f_oei_w)
end

get_tei_integrals("tools/sto3g-water.tei",7)
get_oei_integrals("tools/sto3g-water.oei",7)
