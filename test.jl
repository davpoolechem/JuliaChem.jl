function swap(x,y)
        tmp = x
        x = y
        y = tmp
end

quartet = 0
for i in 1:2, j in 1:i, k in 1:2, l in 1:k
        ij = i*(i+1)/2 + j
        kl = k*(k+1)/2 + l
        if (ij < kl) continue end

        F = [i, j]
        D = [k, l]
        eri = [i, j, k, l]

        #println(F,D,eri)
        println(eri)
        println(i," ",j," <- ",k," ",l)
        println(k," ",l," <- ",i," ",j)
        println(i," ",k," <- ",j," ",l)
        println(i," ",l," <- ",j," ",k)
        println(j," ",k," <- ",i," ",l)
        println(j," ",l," <- ",i," ",j)
end

function disassemble(x,y)
        return x+y
end
