function build_jeri()
  #cd(joinpath(@__DIR__,"src"))
  run(`make`)
end

build_jeri()
