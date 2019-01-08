function bigfact(n::Int64)
    bigint::BigInt = 1
    for i in 1:n
        bigint *= i
    end

    #return bigint
end

inp = 100
bigfact(inp)
