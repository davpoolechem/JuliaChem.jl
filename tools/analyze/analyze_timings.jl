using Statistics

function get_timings_for_file(input_file_string, mean_timings,
  median_timings) 
  #== read in input file ==#
  input_file::Array{String,1} = []
  open(input_file_string) do file
    input_file = readlines(file)
  end

  for time_to_parse in collect(keys(mean_timings))
    #== extract timings from input file ==#
    timings_string::Array{String,1} = [] 
    for line in input_file
      if occursin("$time_to_parse TIME",line)
        push!(timings_string, line)
      end
    end
    #display(timings_string)
   
    timings::Array{Float64,1} = [] 
    for line in timings_string
      timing_match = eachmatch(r"[0-9]+.[0-9]+", line)
      timing = parse(Float64, collect(timing_match)[1].match)
      push!(timings, timing)
    end
    #display(timings)
    #println()

    avg = mean(timings)
    med = median(timings)
    stdev = std(timings)
    reldev = 100*stdev/avg
 
    push!(mean_timings["$time_to_parse"], avg) 
    push!(median_timings["$time_to_parse"], med) 
  end
end
 
#=
function analyze_timing(input_file_string::String) 
  #== read in input file ==#
  input_file::Array{String,1} = []
  open(input_file_string) do file
    input_file = readlines(file)
  end

  for time_to_parse in ["DECOMPOSE", "ASSIGN", "INTEGRALS", "FOCK", 
    "KERNEL SUM", "TWOEI"]
    #== extract timings from input file ==#
    timings_string::Array{String,1} = [] 
    for line in input_file
      if occursin("$time_to_parse TIME",line)
        push!(timings_string, line)
      end
    end
    #display(timings_string)
   
    timings::Array{Float64,1} = [] 
    for line in timings_string
      timing_match = eachmatch(r"[0-9]+.[0-9]+", line)
      timing = parse(Float64, collect(timing_match)[1].match)
      push!(timings, timing)
    end
    #display(timings)
    #println()

    avg = mean(timings)
    med = median(timings)
    stdev = std(timings)
    reldev = 100*stdev/avg
  
    println("Mean $time_to_parse time: $avg")
    println("Median $time_to_parse time: $med")
    println("Standard Deviation $time_to_parse time: $stdev")
    println("Relative Deviation $time_to_parse time: $reldev")
    println()
  end
end
=#

files = [ "01_MP2.log",
          "02_MP2.log",
          "03_MP2.log",
          "04_MP2.log",
          "05_MP2.log",
          "06_MP2.log",
          "07_MP2.log",
          "08_MP2.log",
          "09_MP2.log",
          "10_MP2.log",
          "11_MP2.log",
          "12_MP2.log",
          "13_MP2.log",
          "14_MP2.log",
          "15_MP2.log",
          "16_MP2.log",
          "17_MP2.log",
          "18_MP2.log",
          "19_MP2.log",
          "20_MP2.log",
          "21_MP2.log",
          "22_MP2.log"
        ] 

mean_timings = Dict("DECOMPOSE" => Vector{Float64}([]),
                    "ASSIGN" => Vector{Float64}([]),
                    "INTEGRALS" => Vector{Float64}([]),
                    "FOCK" => Vector{Float64}([]),
                    "KERNEL SUM" => Vector{Float64}([]),
                    "TWOEI" => Vector{Float64}([])
                   )
  
median_timings = Dict("DECOMPOSE" => Vector{Float64}([]),
                      "ASSIGN" => Vector{Float64}([]),
                      "INTEGRALS" => Vector{Float64}([]),
                      "FOCK" => Vector{Float64}([]),
                      "KERNEL SUM" => Vector{Float64}([]),
                      "TWOEI" =>Vector{Float64}([])
                     )

for file in files
  file_location = joinpath(@__DIR__, "example_inputs/S22-profile/$file") 
  get_timings_for_file(file_location, mean_timings, median_timings)
end

for timing_key in collect(keys(mean_timings))
  println("MEAN $timing_key TIMINGS:")
  for timing in mean_timings[timing_key]
    println(timing)
  end
    println()
end 
  
