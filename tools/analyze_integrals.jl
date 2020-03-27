function get_two_electron_integrals(input_file_string::String)
  #== read in input file ==#
  input_file::Array{String,1} = []
  open(input_file_string) do file
    input_file = readlines(file)
  end

  #== extract integrals from input file ==#
  starting::Int64 = 1
  ending::Int64 = 1

  integrals = Vector{String}([])
  for line in input_file
    if occursin("ERI(",line)
      push!(integrals,line)
    end
  end

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
    if typeof(match(r"(-?[\d]{1,}.[\d]{1,}[e()E][+()-][\d]{1,})",integral)) != Nothing
      integral_match = match(r"(-?[\d]{1,}.[\d]{1,}[e()E][+()-][\d]{1,})",integral)
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

julia_dict = get_two_electron_integrals(ARGS[1])
gamess_dict = get_two_electron_integrals(ARGS[2])

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
    println("Missing GAMESS key:", key)
  end 
  
  if has_julia_key
    same = gamess_dict[key] < 1.0E-10 && julia_dict[key] < 1.0E-10 ? true : 
      isapprox(gamess_dict[key], julia_dict[key], rtol = 1e-9)
    if !same
      println("Not same!: ", key, "; ", gamess_dict[key], ", ", julia_dict[key])
    end
  end
end
