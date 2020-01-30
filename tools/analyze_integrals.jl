function get_two_electron_integrals_julia(input_file_string::String)
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
  
  #== set up integral dictionary ==#
  integral_dict = Dict([])
  num_duplicates = 0
  for integral in integrals
    #==get ijkl numbers ==#
    ijkl_match = eachmatch(r"[0-9]+", integral)
    ijkl = collect(ijkl_match)
    i = ijkl[1].match
    j = ijkl[2].match
    k = ijkl[3].match
    l = ijkl[4].match
  
    #== get actual integral ==#
    integral_match = match(r"0","0")
    if typeof(match(r"(-?[\d].[\d]{1,}[e][+()-][\d]{1,})",integral)) != Nothing
      integral_match = match(r"(-?[\d].[\d]{1,}[e()E][+()-][\d]{1,})",integral)
    else
      integral_match = match(r"(-?[\d].[\d]{1,})",integral)
    end
   
    integral_value = parse(Float64,integral_match.match)
    if haskey(integral_dict,"$i,$j,$k,$l")
      println("Duplicate Julia key at: "*"$i,$j,$k,$l")
      num_duplicates += 1
    else
      integral_dict["$i,$j,$k,$l"] = integral_value
    end
  end

  println("Number of duplicate keys: ", num_duplicates)

  return integral_dict
end

function get_two_electron_integrals_gamess(input_file_string::String)
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
 
  #== set up integral dictionary ==#
  integral_dict = Dict([])
  for integral in 1:2:length(integrals)
    #==get ijk numbers ==#
    ijk_match = eachmatch(r"[0-9]+", integrals[integral])
    ijk = collect(ijk_match)
    
    i = ijk[1].match
    j = ijk[2].match
    k = ijk[3].match
   
    l_match = eachmatch(r"[0-9]+", integrals[integral+1])
    l = collect(l_match)[1].match

    #== get actual integral ==#
    integral_match = match(r"0","0")
    if typeof(match(r"(-?[\d].[\d]{1,}[E][+()-][\d]{1,})",integrals[integral+1])) != Nothing
      integral_match = match(r"(-?[\d].[\d]{1,}[E][+()-][\d]{1,})",integrals[integral+1])
    else
      integral_match = match(r"(-?[\d].[\d]{1,})",integrals[integral+1])
    end
    
    integral_value = parse(Float64,integral_match.match)
    integral_dict["$i,$j,$k,$l"] = integral_value
  end

  return integral_dict
end

julia_dict = get_two_electron_integrals_julia(ARGS[1])
gamess_dict = get_two_electron_integrals_gamess(ARGS[2])

println(length(keys(julia_dict)))
println(length(keys(gamess_dict)))

for key in keys(gamess_dict)
  has_julia_key = haskey(julia_dict,key)
  if !has_julia_key
    println("Missing Julia key:", key)
  end 
end

for key in keys(julia_dict)
  has_julia_key = haskey(gamess_dict,key)
  if !has_julia_key
    println("Missing GAMESS key:", key,": ",julia_dict[key])
  end 
  
  #if has_julia_key
  #  same = isapprox(gamess_dict[key], julia_dict[key], rtol = 1e-10)
  #  if !same
  #    println("Not same!: ", gamess_dict[key], ", ", julia_dict[key])
  #  end
  #end
end
