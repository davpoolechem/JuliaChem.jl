function gamess_to_julia_key(gamess_key)
  gamess_to_julia_mapping = Dict([ 1 => 1,
                                   2 => 2,
                                   3 => 4,
                                   4 => 7,
                                   5 => 3,
                                   6 => 5,
                                   7 => 8,
                                   8 => 6,
                                   9 => 9,
                                   10 => 10,
                                   11 => 11,
                                   12 => 12
                                 ])

  ij_match = eachmatch(r"[0-9]+", gamess_key)
  ij = collect(ij_match)
  i = parse(Int,ij[1].match)
  j = parse(Int,ij[2].match)

  i_new = gamess_to_julia_mapping[i]
  j_new = gamess_to_julia_mapping[j]

  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  return "$i_new,$j_new"
end

function julia_to_gamess_key(julia_key)
  julia_to_gamess_mapping = Dict([ 1 => 1,
                                   2 => 2,
                                   3 => 5,
                                   4 => 3,
                                   5 => 6,
                                   6 => 8,
                                   7 => 4,
                                   8 => 7,
                                   9 => 9,
                                   10 => 10,
                                   11 => 11,
                                   12 => 12
                                 ])
  
  ij_match = eachmatch(r"[0-9]+", julia_key)
  ij = collect(ij_match)
  i = parse(Int,ij[1].match)
  j = parse(Int,ij[2].match)

  i_new = julia_to_gamess_mapping[i]
  j_new = julia_to_gamess_mapping[j]

  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  return "$i_new,$j_new"
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
    if occursin("SCHWARZ(",line)
      push!(integrals,line)
    end
  end

  #== set up integral dictionary ==#
  integral_dict = Dict([])
  num_duplicates = 0
  for integral in integrals
    #==get ij numbers ==#
    ij_match = eachmatch(r"[0-9]+", integral)

    ij = collect(ij_match)
    i = ij[1].match
    j = ij[2].match

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
    #if haskey(integral_dict,"$i,$j")
    if false 
      println("Duplicate Julia key at: "*"$i,$j")
      num_duplicates += 1
    else
      integral_dict["$i,$j"] = integral_value
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
      isapprox(gamess_dict[gamess_key], julia_dict[julia_key], atol = 1e-7)
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
      isapprox(gamess_dict[gamess_key], julia_dict[julia_key], atol = 1e-7)
    if !same
      println("Not same!: $julia_key => $gamess_key: ", julia_dict[julia_key], ", ", gamess_dict[gamess_key])
    end
  end
end
