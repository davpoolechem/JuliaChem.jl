module BSEParse

using PyCall
using HDF5

bse = pyimport("basis_set_exchange")

#=================#
#= parse a shell =#
#=================#
function parse_shell(bs_dict::String,multishell::Bool,
    inner_regex_index::Int64,inner_loop_index::Int64)

    #== initialize variables ==#
    exponents::Array{Float64,1} = [ ]
    coeff::Array{Float64} = [ ]

    number_counter::Int64 = 0

    #== use this parsing scheme if not multishell... ==#
    if multishell == false
        #== loop over numbers in shell ==#
        for i::Int64 in 1:inner_loop_index
            #== get current number ==#
            regexmatch::RegexMatch = match(r"0","0")
            non_scientific::Bool = false
            if (typeof(findnext(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                bs_dict,inner_regex_index)) != Nothing) && #there is a next scientific notation number
                (typeof(findnext(r"([-\s][\d].[\d]{1,}[\s])",
                bs_dict,inner_regex_index)) != Nothing) #there is a next non-scientific notation number

                if (findnext(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                    bs_dict,inner_regex_index)[1] > findnext(
                    r"([-\s][\d].[\d]{1,}[\s])",bs_dict,
                    inner_regex_index)[1]) #the non-scientific number is next on the list
                  non_scientific = true
                end
            end

            if (typeof(match(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                bs_dict,inner_regex_index)) == Nothing) || non_scientific

                regexmatch = match(r"([-\s][\d].[\d]{1,}[\s])",
                bs_dict,inner_regex_index)
            else
                regexmatch= match(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                bs_dict,inner_regex_index)
            end

            #== add exponent number to exponent array... ==#
            if number_counter%2 == 0 push!(exponents,
                parse(Float64,regexmatch[1]))

            #== ...or add contraction coefficient to coeff array ==#
            else push!(coeff,parse(Float64,regexmatch[1]))
            end

            #== increment variables for next number in string ==#
            number_counter += 1
            inner_regex_index = findnext(regexmatch[1],bs_dict,
                inner_regex_index)[end]
        end
    #== ...else use this scheme ==#
    else
        #== loop over numbers in shell ==#
        for i::Int64 = 1:inner_loop_index
            #== get current number ==#
            regexmatch::RegexMatch = match(r"0","0")
            non_scientific::Bool = false
            if (typeof(findnext(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                bs_dict,inner_regex_index)) != Nothing) && #there is a next scientific notation number
                (typeof(findnext(r"([-\s][\d].[\d]{1,}[\s])",
                bs_dict,inner_regex_index)) != Nothing) #there is a next non-scientific notation number

                if (findnext(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                    bs_dict,inner_regex_index)[1] > findnext(
                    r"([-\s][\d].[\d]{1,}[\s])",bs_dict,
                    inner_regex_index)[1]) #the non-scientific number is next on the list
                  non_scientific = true
                end
            end

            if (typeof(match(r"([\s()-][\d].[\d]{1,}[E][+()-][\d]{2}[\s])",
                bs_dict,inner_regex_index)) == Nothing) || non_scientific

                regexmatch= match(r"([-\s][\d].[\d]{1,}[\s])",
                bs_dict,inner_regex_index) #match numbers not in scientific notation
            else
                regexmatch= match(r"([-\s][\d].[\d]{1,}[E][+-][\d]{2}[\s])",
                bs_dict,inner_regex_index) #match numbers in scientific notation
            end

            #== add exponent number to exponent array... ==#
            if number_counter%3 == 0 push!(exponents,
                parse(Float64,regexmatch[1]))

            #== ...or add contraction coefficient to coeff array ==#
            else push!(coeff,parse(Float64,regexmatch[1]))
            end

            #== increment variables for next number in string ==#
            number_counter += 1
            inner_regex_index = findnext(regexmatch[1],bs_dict,
                inner_regex_index)[end]
        end

        #== convert coeff array into 2-column matrix if multishell ==#
        coeff2::Array{Float64,2} = reshape(coeff,(2,Int64(length(coeff)/2)))
        coeff = transpose(coeff2)
    end

    return exponents, coeff, inner_regex_index
end

#===========================================================#
#= perform parsing operation on single atom/basis set pair =#
#===========================================================#
function parse_individual(atom::String, basis::String, bsed::HDF5File)
    #== get basis set ==#
    bs_dict::String = ""
    try
        bs_dict = bse.get_basis(basis,fmt="gamess_us",
            elements=[atom], header=false)
    #== if invalid atom/basis set combination, return ==#
    catch
        println("WARNING. The basis set $basis does not exist for atom $atom.")
        return
    end

    println(bs_dict)

    #== initialize some variables ==#
    outer_counter::Int64 = 1

    outer_regex_index::Int64 = 1
    inner_regex_index::Int64 = 1

    #== loop over shells in basis set ==#
    while true
        #== if no more shells, we are done ==#
        if typeof(match(r"([A-Z]\s+[0-9])",
            bs_dict,outer_regex_index)) == Nothing
            break
        end

        #== determine if next shell is a multishell==#
        multishell::Bool = false

        regexmatch::RegexMatch = match(r"([A-Z]\s+[0-9])",
            bs_dict,outer_regex_index)
        if regexmatch[1][1] == 'L' multishell = true end

        shell_type_char::Char = regexmatch[1][1]
        shell_type_string::String = "$shell_type_char"

        #== determine number of elements to loop over in shell ==#
        inner_loop_index::Int64 = multishell == true ?
            3 * parse(Int64, regexmatch[1][end]) :
            2 * parse(Int64, regexmatch[1][end])

        #== parse the shell elements ==#
        exponents::Array{Float64,1}, coeff::Array{Float64},
            inner_regex_index = parse_shell(bs_dict,
            multishell, inner_regex_index, inner_loop_index)

        #== write shell info to hdf5 database ==#
        shell_num::String = string(outer_counter)

        h5write("bsed.h5",
            "$atom/$basis/$shell_num/Shell Type", shell_type_string)
        h5write("bsed.h5",
            "$atom/$basis/$shell_num/Exponents", exponents)
        h5write("bsed.h5",
            "$atom/$basis/$shell_num/Coefficients", coeff)


        println("Exponents:")
        display(h5read("bsed.h5",
            "$atom/$basis/$shell_num/Exponents"))
        println("")
        println("Contraction Coefficients:")
        display(h5read("bsed.h5",
            "$atom/$basis/$shell_num/Coefficients"))
        println("")


        #== set index of next shell to loop over ==#
        outer_regex_index = findnext(regexmatch[1], bs_dict,
            outer_regex_index)[end]
        outer_counter += 1
    end
end

#===========================================#
#= parse all selected atom-basis set pairs =#
#===========================================#
function parse_all()
    #== open HDF5 file for writing ==#
    hdf5name::String = "bsed.h5"
    h5open(hdf5name, "w") do bsed

        #== parse "sto" basis family ==#
        atoms::Array{String,1} = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe"
            ] #H-Xe
        #=
        atoms::Array{String,1} = [
            "1", "2",
            "3", "4", "5", "6", "7", "8", "9", "10",
            "11", "12", "13", "14", "15", "16", "17", "18",
            "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
            "30", "31", "32", "33", "34", "35", "36",
            "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47",
            "48", "49", "50", "51", "52", "53"
            ] #H-Xe
        =#
        basis_sets::Array{String,1} = ["STO-2G", "STO-3G", "STO-4G", "STO-5G",
            "STO-6G" ] #STO family

        for atom::String in atoms
            for basis::String in basis_sets
                parse_individual(atom, basis, bsed)
            end
        end

        #== parse "pople" basis family ==#
        atoms = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe"
            ] #H-Xe
        basis_sets = ["3-21G", "4-31G", "5-21G", "5-31G", "6-21G",
            "6-31G","6-31G(d,p)" ] #pople family

        for atom::String in atoms
            for basis::String in basis_sets
                parse_individual(atom, basis, bsed)
            end
        end

        #== parse polarization-consistent basis family ==#
        atoms = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe"
            ] #H-Xe
        basis_sets = ["PCSeg-0"] #polarization-consistent family

        for atom::String in atoms
            for basis::String in basis_sets
                parse_individual(atom, basis, bsed)
            end
        end
    
        #== parse correlation-consistent basis family ==#
        atoms = [
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe"
            ] #H-Xe
        basis_sets = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"] 

        for atom::String in atoms
            for basis::String in basis_sets
                parse_individual(atom, basis, bsed)
            end
        end

    end
end

end

BSEParse.parse_all()
