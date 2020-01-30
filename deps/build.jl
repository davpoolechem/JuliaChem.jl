function build(ARGS)
  if (ARGS == [])
    SIMINT = ENV["SIMINT"] 
    run(`cmake -DSIMINT_PATH=$SIMINT .`)
    run(`make`)
  else
    SIMINT = ARGS[1] 
    run(`cmake -DSIMINT_PATH=$SIMINT .`)
    run(`make`)
  end
end

build(ARGS)
