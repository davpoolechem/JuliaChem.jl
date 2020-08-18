function build()
  cd(joinpath(@__DIR__,"src"))
  run(`make`)
end

build()
