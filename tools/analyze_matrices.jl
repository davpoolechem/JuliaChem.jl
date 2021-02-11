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

  ij_match = eachmatch(r"[0-9]+", gamess_key)
  ij = collect(ij_match)
  i = parse(Int,ij[1].match)
  j = parse(Int,ij[2].match)

  #i_new = gamess_to_julia_mapping[i]
  #j_new = gamess_to_julia_mapping[j]

  i_new = i 
  j_new = j 


  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  return "$i_new,$j_new"
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
  
  ij_match = eachmatch(r"[0-9]+", julia_key)
  ij = collect(ij_match)
  i = parse(Int,ij[1].match)
  j = parse(Int,ij[2].match)

  #i_new = julia_to_gamess_mapping[i]
  #j_new = julia_to_gamess_mapping[j]

  i_new = i 
  j_new = j 

  if i_new <= j_new
    tmp = i_new
    i_new = j_new
    j_new = tmp
  end

  return "$i_new,$j_new"
end

function get_matrix(input_file_string::String, matrix_string::String)
  #== read in input file ==#
  input_file::Array{String,1} = []
  open(input_file_string) do file
    input_file = readlines(file)
  end

  #== extract matrixs from input file ==#
  starting::Int64 = 1
  ending::Int64 = 1

  matrix = Vector{String}([])
  for line in input_file
    if occursin("$matrix_string(",line)
      push!(matrix,line)
    end
  end

  #== set up matrix dictionary ==#
  matrix_dict = Dict([])
  num_duplicates = 0
  for element in matrix
    #==get ij numbers ==#
    ij_match = eachmatch(r"[0-9]+", element)
    ij = collect(ij_match)
    i = ij[1].match
    j = ij[2].match

    #== get actual matrix ==#
    matrix_match = match(r"0","0")
    if typeof(match(r"(-?[\d]{1,}\.[\d]{1,}[e()E][+()-][\d]{1,})",element)) != Nothing
      matrix_match = match(r"(-?[\d]{1,}\.[\d]{1,}[e()E][+()-][\d]{1,})",element)
    elseif typeof(match(r"(-\.[\d]{1,})", element)) != Nothing
      matrix_match = match(r"(-\.[\d]{1,})",element)
    else
      matrix_match = match(r"(-?[\d]\.[\d]{1,})",element)
    end
  
    matrix_value = parse(Float64,matrix_match.match)
    if haskey(matrix_dict,"$i,$j")
      println("Duplicate key at: "*"$i,$j")
      num_duplicates += 1
    else
      matrix_dict["$i,$j"] = matrix_value
    end
  end

  println("Number of duplicate keys: ", num_duplicates)
  
  return matrix_dict
end

julia_dict = get_matrix(ARGS[1], ARGS[3])
gamess_dict = get_matrix(ARGS[2], ARGS[3])

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
