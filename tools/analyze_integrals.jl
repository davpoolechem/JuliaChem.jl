function gamess_to_julia_key(gamess_key)
  gamess_to_julia_mapping = Dict([ 1 => 1,
                                   2 => 2,
                                   3 => 3,
                                   4 => 4,
                                   5 => 5,
                                   6 => 6,
                                   7 => 7,
                                   8 => 8,
                                   9 => 9,
                                   10 => 10,
                                   11 => 13,
                                   12 => 15,
                                   13 => 11,
                                   14 => 12,
                                   15 => 14,
                                   16 => 16,
                                   17 => 17,
                                   18 => 18,
                                   19 => 19,
                                   20 => 20 
                                 ])

  ijkl_match = eachmatch(r"[0-9]+", gamess_key)
  ijkl = collect(ijkl_match)
  i = parse(Int,ijkl[1].match)
  j = parse(Int,ijkl[2].match)
  k = parse(Int,ijkl[3].match)
  l = parse(Int,ijkl[4].match)

  i_new = gamess_to_julia_mapping[i]
  j_new = gamess_to_julia_mapping[j]
  k_new = gamess_to_julia_mapping[k]
  l_new = gamess_to_julia_mapping[l]

  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  if k_new <= l_new
    tmp = k_new
    k_new = l_new
    l_new = tmp
  end

  ij = i_new > j_new ? i_new*(i_new-1)/2 + j_new : j_new*(j_new-1)/2 + i_new
  kl = k_new > l_new ? k_new*(k_new-1)/2 + l_new : l_new*(l_new-1)/2 + k_new

  if ij <= kl
    tmp = i_new
    i_new = k_new
    k_new = tmp

    tmp = j_new
    j_new = l_new
    l_new = tmp
  end

  return "$i_new,$j_new,$k_new,$l_new"
  #return "$i,$j,$k,$l"
end

function julia_to_gamess_key(julia_key)
  julia_to_gamess_mapping = Dict([ 1 => 1,
                                   2 => 2,
                                   3 => 3,
                                   4 => 4,
                                   5 => 5,
                                   6 => 6,
                                   7 => 7,
                                   8 => 8,
                                   9 => 9,
                                   10 => 10,
                                   13 => 11,
                                   15 => 12,
                                   11 => 13,
                                   12 => 14,
                                   14 => 15,
                                   16 => 16,
                                   17 => 17,
                                   18 => 18,
                                   19 => 19,
                                   20 => 20 
                                 ])
  
  ijkl_match = eachmatch(r"[0-9]+", julia_key)
  ijkl = collect(ijkl_match)
  i = parse(Int,ijkl[1].match)
  j = parse(Int,ijkl[2].match)
  k = parse(Int,ijkl[3].match)
  l = parse(Int,ijkl[4].match)

  i_new = julia_to_gamess_mapping[i]
  j_new = julia_to_gamess_mapping[j]
  k_new = julia_to_gamess_mapping[k]
  l_new = julia_to_gamess_mapping[l]

  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  if k_new <= l_new
    tmp = k_new
    k_new = l_new
    l_new = tmp
  end

  ij = i_new > j_new ? i_new*(i_new-1)/2 + j_new : j_new*(j_new-1)/2 + i_new
  kl = k_new > l_new ? k_new*(k_new-1)/2 + l_new : l_new*(l_new-1)/2 + k_new

  if ij <= kl
    tmp = i_new
    i_new = k_new
    k_new = tmp

    tmp = j_new
    j_new = l_new
    l_new = tmp
  end

  return "$i_new,$j_new,$k_new,$l_new"
  #return "$i,$j,$k,$l"
end

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
    if typeof(match(r"(-?[\d]{1,}\.[\d]{1,}[e()E][+()-][\d]{1,})",integral)) != Nothing
      integral_match = match(r"(-?[\d]{1,}\.[\d]{1,}[e()E][+()-][\d]{1,})",integral)
    elseif typeof(match(r"(-\.[\d]{1,})", integral)) != Nothing
      integral_match = match(r"(-\.[\d]{1,})",integral)
    else
      integral_match = match(r"(-?[\d]\.[\d]{1,})",integral)
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

for gamess_key in keys(gamess_dict)
  julia_key = gamess_to_julia_key(gamess_key)
  has_julia_key = haskey(julia_dict,julia_key)
  if !has_julia_key
    println("Missing Julia key:", julia_key)
  else
    same = gamess_dict[gamess_key] < 1.0E-10 && julia_dict[julia_key] < 1.0E-10 ? true : 
      isapprox(gamess_dict[gamess_key], julia_dict[julia_key], atol = 1e-9)
    if !same
      println("Not same!: $gamess_key => $julia_key: ", gamess_dict[gamess_key], ", ", julia_dict[julia_key])
    end
  end 
end

for julia_key in keys(julia_dict)
  gamess_key = julia_to_gamess_key(julia_key) 
  has_gamess_key = haskey(gamess_dict,gamess_key)
  if !has_gamess_key
    println("Missing GAMESS key:", gamess_key)
  else
    same = gamess_dict[gamess_key] < 1.0E-10 && julia_dict[julia_key] < 1.0E-10 ? true : 
      isapprox(gamess_dict[gamess_key], julia_dict[julia_key], atol = 1e-9)
    if !same
      println("Not same!: $julia_key => $gamess_key: ", julia_dict[julia_key], ", ", gamess_dict[gamess_key])
    end
  end
end
