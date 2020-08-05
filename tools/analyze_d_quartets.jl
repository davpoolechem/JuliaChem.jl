function julia_to_gamess_key(julia_key)
  julia_to_gamess_mapping = Dict( 7 => 1,
                                  3 => 2,
                                  1 => 3,
                                  9 => 4,
                                  5 => 5,
                                  8 => 6,
                                  4 => 7,
                                  2 => 8,
                                  6 => 9
                                )
 
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
end

function gamess_to_julia_key(gamess_key)
  gamess_to_julia_mapping = Dict( 1 => 7,
                                  2 => 3,
                                  3 => 1,
                                  4 => 9,
                                  5 => 5,
                                  6 => 8,
                                  7 => 4,
                                  8 => 2,
                                  9 => 6
                                )
 
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
end

function get_d_quartet_integrals(input_file_string::String)
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
    if occursin("QUARTET(",line)
      push!(integrals,line)
    end
  end

  #== set up integral dictionary ==#
  quartet_dict = Dict([])
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
    if !haskey(quartet_dict, "$i,$j,$k,$l")
      quartet_dict["$i,$j,$k,$l"] = [ integral_value ]
    else 
      push!(quartet_dict["$i,$j,$k,$l"], integral_value)
    end
  end
  
  return quartet_dict
end

function check_lists(julia_dict, gamess_dict)
  julia_key_amount = length(keys(julia_dict))
  gamess_key_amount = length(keys(gamess_dict))

  same_size = true
  if julia_key_amount != gamess_key_amount
    println("The lists are of different size!: $julia_key_amount, $gamess_key_amount")
    same_size = false
  else
    println("Both lists have $julia_key_amount quartets. You are good.")
  end

  return same_size
end

function check_gamess_keys_integrals(julia_dict, gamess_dict)
  for gamess_quartet in keys(gamess_dict)
    julia_quartet = gamess_to_julia_key(gamess_quartet)
    if !haskey(julia_dict, julia_quartet)
      println("Missing Julia key! $gamess_quartet => $julia_quartet")
      continue
    end

    julia_value = julia_dict[julia_quartet]
    gamess_value = gamess_dict[gamess_quartet]
     
    julia_array_sorted = sort(julia_value)
    gamess_array_sorted = sort(gamess_value)
    
    @assert length(julia_value) == length(gamess_value)
    same = true
    for idx in 1:length(julia_value)
      same = same && isapprox(julia_array_sorted[idx], 
        gamess_array_sorted[idx]; atol=1e-6) 
    end
    if !same
      println("Different quartet! $julia_quartet and $gamess_quartet")
    end
  end
end

function check_julia_keys_integrals(julia_dict, gamess_dict)
  for julia_quartet in keys(julia_dict)
    gamess_quartet = julia_to_gamess_key(julia_quartet)
    if !haskey(gamess_dict, gamess_quartet)
      println("Missing GAMESS key! $julia_quartet => $gamess_quartet")
      continue  
    end

    julia_value = julia_dict[julia_quartet]
    gamess_value = gamess_dict[gamess_quartet]
     
    julia_array_sorted = sort(julia_value)
    gamess_array_sorted = sort(gamess_value)
    
    @assert length(julia_value) == length(gamess_value)
    same = true
    for idx in 1:length(julia_value)
      same = same && isapprox(julia_array_sorted[idx], 
        gamess_array_sorted[idx]; atol=1e-6) 
    end
    if !same
      println("Different quartet! $julia_quartet and $gamess_quartet")
    end
  end
end

julia_dict = get_d_quartet_integrals(ARGS[1])
gamess_dict = get_d_quartet_integrals(ARGS[2])

@assert check_lists(julia_dict, gamess_dict)
check_gamess_keys_integrals(julia_dict, gamess_dict)
check_julia_keys_integrals(julia_dict, gamess_dict)
