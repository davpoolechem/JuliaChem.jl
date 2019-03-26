#------------------------------#
#             IO.jl            #
#------------------------------#
"""
    readin(input::String)
Summary
======
Read in input file and translate it to an input file object.

Arguments
======
input = name of input file to read in
"""
function readin(input::String)
    inp::IOStream = open(input)
    file::Array{String,1} = readlines(input)
    close(inp)

    return file
end

"""
    geomin(input::Array{String,1})
Summary
======
Extract geometric coordinates from input file.

Arguments
======
input = name of input file object to process
"""
function geomin(input::Array{String,1})
    geom::Int64 = findnext(input.=="\$\$GEOM",1)
    natoms::Int64 = parse(Int64, input[geom+1])

    coord::Matrix{Float64} = Matrix{Float64}(undef,natoms,4)
    for i::Int64 in [1,2...,natoms]
        coord[i,:] = parse.(Float64, split(input[i+geom+1]))
    end

    return coord
end

#extract coordinates from input file
#function flagsin(input::Array{String,1})
#end

"""
    enucin(input::Array{String,1})
Summary
======
Extract nuclear repulsion energy value from input data file.

Arguments
======
input = name of input file object to process
"""
function enucin(input::Array{String,1})
    enuc_loc::Int64 = findnext(input.=="\$\$ENUC",1)
    enuc::Float64 = parse(Float64, input[enuc_loc+1])
    return enuc
end

"""
    oeiin(input::Array{String,1}, type::String)
Summary
======
Extract one-electron integrals from input data file. Type dictates which
one-electron integrals are returned.

Arguments
======
input = name of input file object to process

type = tag for determining which type of one-electron integrals are returned:
1. type=OVR returns overlap integrals.
2. type=KEI returns kinetic energy integrals.
3. type=NAI returns nuclear-electron attraction integrals.
"""
function oeiin(input::Array{String,1}, type::String)
    loc::Int64 = findnext(input.=="\$\$$type",1)
    nbf::Int64 = 7
    nbf2::Int64 = nbf*(nbf+1)/2

    oei = Matrix{Float64}(undef,nbf,nbf)
    ij::Int64 = 1;
    for i::Int64 in 1:1:nbf, j::Int64 in 1:1:i
        ij = i*(i-1)/2 + j
        oei[i,j] = parse(Float64, split(input[ij+loc])[3])
        oei[j,i] = oei[i,j]
    end

    return oei
end

"""
    teiin(input::Array{String,1})
Summary
======
Extract two-electron integrals from input data file.

Arguments
======
input = name of input file object to process
"""
function teiin(input::Array{String,1})
    loc::Int64 = findnext(input.=="\$\$TEI",1)
    nint::Int64 = 228

    tei::Array{Float64,1} = zeros(2401)
    for index::Int64 in 1:1:nint
        i::Int64 = parse(Int64, split(input[index+loc])[1])
        j::Int64 = parse(Int64, split(input[index+loc])[2])
        k::Int64 = parse(Int64, split(input[index+loc])[3])
        l::Int64 = parse(Int64, split(input[index+loc])[4])

        ij::Int64 = (i > j) ? i*(i+1)/2 + j : j*(j+1)/2 + i
        kl::Int64 = (k > l) ? k*(k+1)/2 + l : l*(l+1)/2 + k
        ijkl::Int64 = (ij > kl) ? ij*(ij+1)/2 + kl : kl*(kl+1)/2 + ij

        tei[ijkl] = parse(Float64, split(input[index+loc])[5])
    end

    return tei
end

"""
    teiin(input::Array{String,1})
Summary
======
Re-print input file to display.

Arguments
======
input = name of input file object to process
"""
function iocat(input::String)
    inp::IOStream = open(input)

    file::Array{String,1} = readlines(input)

    for line in file
        println(line)
    end

    #flags::Int64 = findnext(file.=="\$FLAGS",1)
    #geom::Int64 = findnext(file.=="\$GEOM",flags)
    #vec::Int64 = findnext(file.=="\$VEC",geom)

    close(inp)
end

#function compile()
#    precompile(io.readin, tuple(String))
#    precompile(io.geomin, tuple(Array{String,1}))
#    precompile(io.enucin, tuple(Array{String,1}))
#    precompile(io.oeiin, tuple(Array{String,1}, String))
#    precompile(io.teiin, tuple(Array{String,1}))
#end

#we want to precompile all involved modules to reduce cold runs
#include("./snoop/precompile_Base.jl")
#_precompile_base()
#include("./snoop/precompile_Core.jl")
#_precompile_core()
#include("./snoop/precompile_hf.jl")
#_precompile_()
#include("./snoop/precompile_io.jl")
#_precompile_()
#include("./snoop/precompile_LinearAlgebra.jl")
#_precompile_()
#include("./snoop/precompile_openchem.jl")
#_precompile_()
#include("./snoop/precompile_unknown.jl")
#_precompile_()

#------------------------------#
#            Math.jl           #
#------------------------------#
import LinearAlgebra

"""
    Method #1 of ∑:
    ∑(array::Array{Float64})
Summary
======
Return sum of elements in an array.

Arguments
======
array = array to be summed across
"""
function ∑(array::Array{Float64})
    return sum(array)
end

"""
    Method #2 of ∑:
    ∑(array_1::Array{Float64},array_2::Array{Float64})
Summary
======
Return dot product of two arrays.

Arguments
======
array_1 = first array in the dot product

array_2 = second array in the dot product

"""
function ∑(array_1::Array{Float64},array_2::Array{Float64})
    return LinearAlgebra.dot(array_1,array_2)
end

#------------------------------#
#             HF.jl            #
#------------------------------#
"""
    Data
Summary
======
Core Hartree-Fock data structures

Fields
======
1. Fock::Array{Float64,2} = Fock Matrix
2. Density::Array{Float64,2} = Density Matrix
3. Coeff::Array{Float64,2} = Molecular Orbital Coefficient Matrix
4. Energy::Float64 = Electronic Energy
"""
mutable struct Data
    Fock::Array{Float64,2}
    Density::Array{Float64,2}
    Coeff::Array{Float64,2}
    Energy::Float64
end

ioff = map((x) -> x*(x+1)/2, collect(1:7*8))

"""
     energy(dat::Array{String,1})
Summary
======
Perform the core RHF SCF algorithm.

Arguments
======
dat = Input data file object
"""
function energy(dat::Array{String,1})
    println("----------------------------------------")
    println("  RESTRICTED CLOSED-SHELL HARTREE-FOCK"  )
    println("----------------------------------------")
    println("")

    scf::Data = Data(zeros(7,7), zeros(7,7), zeros(7,7), 0)

    #Step #1: Nuclear Repulsion Energy
    E_nuc::Float64 = enucin(dat)

    #Step #2: One-Electron Integrals
    S::Matrix{Float64} = oeiin(dat,"OVR")
    T::Matrix{Float64} = oeiin(dat,"KEI")
    V::Matrix{Float64} = oeiin(dat,"NAI")
    H::Matrix{Float64} = T+V

    #Step #3: Two-Electron Integrals
    tei::Array{Float64} = teiin(dat)

    #Step #4: Build the Orthogonalization Matrix
    #println("Initial S matrix:")
    #display(S)
    #println("")
    S_evec::Array{Float64} = LinearAlgebra.eigvecs(S)

    S_eval_diag::Array{Float64} = LinearAlgebra.eigvals(S)
    #println("Initial S_evec matrix:")
    #display(S_evec)
    #println("")
    S_eval::Array{Float64} = zeros(7,7)
    for i::Int64 in 1:7
        S_eval[i,i] = S_eval_diag[i]
    end

    ortho::Array{Float64} = S_evec*(S_eval^-0.5)*transpose(S_evec)

    #Step #5: Build the Initial (Guess) Density
    F::Array{Float64} = transpose(ortho)*H*ortho
    D::Array{Float64} = zeros(7,7)
    C::Array{Float64} = zeros(7,7)

    println("----------------------------------------")
    println("        Starting SCF iterations...")
    println("----------------------------------------")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")

    F, D, C, E_elec = iteration(F, D, H, ortho)
    E::Float64 = E_elec + E_nuc

    println(0,"     ", E)

    #start scf cycles: #7-10
    converged::Bool = false
    iter::Int64 = 1
    while(!converged)

        #Step #7: Compute the New Fock Matrix
        F = twoei(F, D, tei, H)
        #println("Initial Fock matrix:")
        #display(F)
        #println("")

        #Step #8: Build the New Density Matrix
        D_old::Array{Float64,2} = deepcopy(D)
        E_old::Float64 = E

        F = transpose(ortho)*F*ortho
        F, D, C, E_elec = iteration(F, D, H, ortho)
        E = E_elec+E_nuc

        #Step #10: Test for Convergence
        ΔE::Float64 = E - E_old

        ΔD::Array{Float64,2} = D - D_old
        D_rms::Float64 = √(∑(ΔD,ΔD))

        println(iter,"     ", E,"     ", ΔE,"     ", D_rms)

        converged = (ΔE <= 1E-6) && (D_rms <= 1E-4)
        iter += 1
        if (iter > 50) break end
    end

    println(" ")
    println("----------------------------------------")
    println("The SCF calculation has converged!")
    println("----------------------------------------")
    println("Total SCF Energy: ",E," h")
    println(" ")
    println("----------------------------------------")
    println("      END RESTRICTED CLOSED-SHELL"       )
    println("              HARTREE-FOCK"              )
    println("----------------------------------------")
    println(" ")

    scf.Fock = F
    scf.Density = D
    scf.Coeff = C
    scf.Energy = E

    return scf
end

"""
     properties(scf::Data)
Summary
======
Compute properties for RHF wave function.

Arguments
======
scf = Core HF data structures
"""
function properties(scf::Data)
    println("----------------------------------------")
    println("            SYSTEM PROPERTIES"           )
    println("----------------------------------------")
    println("")
    println("----------------------------------------")
    println("          END SYSTEM PROPERTIES"         )
    println("----------------------------------------")
end

"""
     iteration(F::Array{Float64,2}, D::Array{Float64,2}, H::Array{Float64,2}, ortho::Array{Float64,2})
Summary
======
Perform single SCF iteration.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

H = One-electron Hamiltonian Matrix

ortho = Symmetric Orthogonalization Matrix
"""
function iteration(F::Array{Float64,2}, D::Array{Float64,2}, H::Array{Float64,2}, ortho::Array{Float64,2})

    #Step #8: Build the New Density Matrix
    F_eval::Array{Float64,1} = LinearAlgebra.eigvals(F)

    F_evec::Array{Float64,2} = LinearAlgebra.eigvecs(F)
    F_evec = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

    C::Array{Float64,2} = ortho*F_evec

    for i::Int64 in 1:7, j::Int64 in 1:i
        D[i,j] = ∑(C[i,1:5],C[j,1:5])
        D[j,i] = D[i,j]
    end

    #Step #9: Compute the New SCF Energy
    E_elec::Float64 = ∑(D,H + F)

    return (F, D, C, E_elec)
end

"""
     index(a::Int64,b::Int64)
Summary
======
Triangular indexing determination.

Arguments
======
a = row index

b = column index
"""
function index(a::Int64,b::Int64)
    index::Int64 = (a > b) ? ioff[a] + b : ioff[b] + a
    return index
end
#=
function twoei_base(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
    F = deepcopy(H)
    J::Array{Float64,2} = zeros(7,7)
    K::Array{Float64,2} = zeros(7,7)

    for i::Int64 in 1:7
        for j::Int64 in 1:7
            for k::Int64 in 1:7
                for l::Int64 in 1:7
                    ij::Int64 = index(i,j)
                    kl::Int64 = index(k,l)
                    ijkl::Int64 = index(ij,kl)
                    #J_idx = 7*(μ-1)+ν
                    #F[i,j] += D[k,l] * (2*tei[ijkl] - tei[ikjl])
                    #println(ijkl,"        ",tei[ijkl])
                    J[i,j] += D[k,l] * tei[ijkl]
                    K[i,k] += D[j,l] * tei[ijkl]
                end
            end
        end
    end
    F += 2*J - K
    return F
end

function twoei_ij(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
    F = deepcopy(H)
    J::Array{Float64,2} = zeros(7,7)
    K::Array{Float64,2} = zeros(7,7)
    for i::Int64 in 1:7
        for j::Int64 in 1:i
            for k::Int64 in 1:7
                for l::Int64 in 1:7
                    ij::Int64 = index(i,j)
                    kl::Int64 = index(k,l)
                    ijkl::Int64 = index(ij,kl)
                    #J_idx = 7*(μ-1)+ν
                    #F[i,j] += D[k,l] * (2*tei[ijkl] - tei[ikjl])
                    #F[i,j] += 2 * D[k,l] * tei[ijkl]
                    #F[i,k] += -1 * D[j,l] * tei[ijkl]
                    val = (i == j) ? 0.5 : 1
                    J[i,j] += val * D[k,l] * tei[ijkl]
                    J[j,i] += val * D[k,l] * tei[ijkl]
                    #J[k,l] += val * D[i,j] * tei[ijkl]
                    #J[l,k] += val * D[i,j] * tei[ijkl]
                    if (i == j)
                        K[i,k] += D[j,l] * tei[ijkl]
                    else
                        K[i,k] += val * D[j,l] * tei[ijkl]
                        K[j,k] += val * D[i,l] * tei[ijkl]
                    end
                end
            end
            #F[j,i] = F[i,j]
        end
    end
    F += 2*J - K
    return F
end

function twoei_kl(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
    F = deepcopy(H)
    #F = zeros(7,7)
    J::Array{Float64,2} = zeros(7,7)
    K::Array{Float64,2} = zeros(7,7)
    for i::Int64 in 1:7
        for j::Int64 in 1:i
            for k::Int64 in 1:7
                for l::Int64 in 1:k
                    ij::Int64 = index(i,j)
                    kl::Int64 = index(k,l)
                    ijkl::Int64 = index(ij,kl)
                    #val = (k == l) ? 0.5 : 1
                    #F[i,j] += val * D[k,l] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
                    val = (i == j) ? 0.5 : 1
                    val *= (k == l) ? 0.5 : 1
                    #J[i,j] += val * D[k,l] * tei[ijkl]
                    #J[j,i] += val * D[k,l] * tei[ijkl]
                    #J[i,j] += val * D[l,k] * tei[ijkl]
                    #J[j,i] += val * D[l,k] * tei[ijkl]
                    J[k,l] += val * D[i,j] * tei[ijkl]
                    J[k,l] += val * D[j,i] * tei[ijkl]
                    J[l,k] += val * D[i,j] * tei[ijkl]
                    J[l,k] += val * D[j,i] * tei[ijkl]
                    K[i,k] += val * D[j,l] * tei[ijkl]
                    K[i,l] += val * D[j,k] * tei[ijkl]
                    K[j,k] += val * D[i,l] * tei[ijkl]
                    K[j,l] += val * D[i,k] * tei[ijkl]
                end
            end
            #F[j,i] = F[i,j]
        end
    end
    #F += H
    F += 2*J - K
    return F
end

function twoei_kl_opt(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
    #F = deepcopy(H)
    F = zeros(7,7)
    J::Array{Float64,2} = zeros(7,7)
    K::Array{Float64,2} = zeros(7,7)
    for i::Int64 in 1:7
        for j::Int64 in 1:i
            ij::Int64 = index(i,j)
            for k::Int64 in 1:7
                for l::Int64 in 1:k
                    kl::Int64 = index(k,l)
                    ijkl::Int64 = index(ij,kl)
                    val = (i == j) ? 0.5 : 1
                    val *= (k == l) ? 0.5 : 1
                    J[k,l] += 2 * val * D[i,j] * tei[ijkl]
                    J[l,k] += 2 * val * D[i,j] * tei[ijkl]
                    K[i,k] += val * D[j,l] * tei[ijkl]
                    K[i,l] += val * D[j,k] * tei[ijkl]
                    K[j,k] += val * D[i,l] * tei[ijkl]
                    K[j,l] += val * D[i,k] * tei[ijkl]
            #        if (i == j && k == l)
        #                K[i,k] += D[j,l] * tei[ijkl]
    #                elseif (i == j)
#                        K[i,k] += D[j,l] * tei[ijkl]
                     #    K[i,l] += D[j,k] * tei[ijkl]
                     #elseif (k == l)
            #            K[i,k] += D[j,l] * tei[ijkl]
            #            K[j,k] += D[i,l] * tei[ijkl]
            #        else
            #        end
                end
            end
            #F[j,i] = F[i,j]
        end
    end
    #F += H
    F = H + 2*J - K
    return F
end

function twoei_fock2a(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
    #F = deepcopy(H)
    F = zeros(7,7)
    J::Array{Float64,2} = zeros(7,7)
    Ki::Array{Float64,2} = zeros(7,7)
    Kj::Array{Float64,2} = zeros(7,7)
    for i::Int64 in 1:7
        for j::Int64 in 1:i
            ij::Int64 = index(i,j)
            for k::Int64 in 1:7
                for l::Int64 in 1:k
                    kl::Int64 = index(k,l)
                    ijkl::Int64 = index(ij,kl)
                    val = (i == j) ? 0.5 : 1
                    val *= (k == l) ? 0.5 : 1
                    J[k,l] += 2 * val * D[i,j] * tei[ijkl]
                    J[l,k] += 2 * val * D[i,j] * tei[ijkl]
                    if (k == l)
                        Ki[i,k] += 2 * val * D[j,l] * tei[ijkl]
                        Kj[j,k] += 2 * val * D[i,l] * tei[ijkl]
                    else
                        Ki[i,k] += val * D[j,l] * tei[ijkl]
                        Ki[i,l] += val * D[j,k] * tei[ijkl]
                        Kj[j,k] += val * D[i,l] * tei[ijkl]
                        Kj[j,l] += val * D[i,k] * tei[ijkl]
                    end
            #        if (i == j && k == l)
        #                K[i,k] += D[j,l] * tei[ijkl]
    #                elseif (i == j)
#                        K[i,k] += D[j,l] * tei[ijkl]
                       #  K[i,l] += D[j,k] * tei[ijkl]
                     #elseif (k == l)
            #            K[i,k] += D[j,l] * tei[ijkl]
            #            K[j,k] += D[i,l] * tei[ijkl]
            #        else
            #        end
                end
            end
            #F[j,i] = F[i,j]
        end
    end
    #F += H
    F = H + 2*J - Ki - Kj
    return F
end
=#

"""
     twoei(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""
function twoei(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
    F = deepcopy(H)
    #F = zeros(7,7)
    coulomb::Array{Float64} = zeros(7,7)
    exchange::Array{Float64} = zeros(7,7)
    for μ::Int64 in 1:7
        μμ = ioff[μ]
        for ν::Int64 in 1:μ
            νν = ioff[ν]
            μν::Int64 = μμ + ν
            for λ::Int64 in 1:7
                λλ = ioff[λ]
                μλ::Int64 = (μ > λ) ? μμ + λ : μ + λλ
                νλ::Int64 = (ν > λ) ? νν + λ : ν + λλ
                for σ::Int64 in 1:λ
                    σσ = ioff[σ]
                    μσ::Int64 = (μ > σ) ? μμ + σ : μ + σσ
                    νσ::Int64 = (ν > σ) ? νν + σ : ν + σσ
                    λσ::Int64 = (λ > σ) ? λλ + σ : λ + σσ
                    μνλσ::Int64 = index(μν,λσ)
                    μλνσ::Int64 = index(μλ,νσ)
                    μσνλ::Int64 = index(μσ,νλ)
                    #eri = tei[ijkl]
                    #eri1 = 2*tei[ijkl] - tei[ikjl]
                    #eri2 = 2*tei[ijkl] - tei[iljk]
                    #eri = eri1 + eri2
                    val = (λ == σ) ? 0.5 : 1
                    #if (k == l)
                #        eri1 *= 0.5
            #            eri2 *= 0.5
            #        end
                    #F[i,j] += D[k,l] * eri1
                    #F[i,j] += D[l,k] * eri2
                    coulomb[μ,ν] += val*D[λ,σ]*tei[μνλσ]
                    exchange[μ,ν] += val*D[λ,σ]*tei[μλνσ]
                    coulomb[μ,ν] += val*D[σ,λ]*tei[μνλσ]
                    exchange[μ,ν] += val*D[σ,λ]*tei[μσνλ]
                end
            end
            coulomb[ν,μ] = coulomb[μ,ν]
            exchange[ν,μ] = exchange[μ,ν]
            #F[j,i] = F[i,j]
        end
    end
    F += 2*coulomb - exchange
    return F
end

#=
function twoei_edit(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
    F = deepcopy(H)
    #F = zeros(7,7)
    for i::Int64 in 1:7, j::Int64 in 1:i
        for k::Int64 in 1:7, l::Int64 in 1:k
            ij::Int64 = index(i,j)
            kl::Int64 = index(k,l)
            ijkl::Int64 = index(ij,kl)
            ik::Int64 = index(i,k)
            jl::Int64 = index(j,l)
            ikjl::Int64 = index(ik,jl)
            il::Int64 = index(i,l)
            jk::Int64 = index(j,k)
            iljk::Int64 = index(il,jk)
            eri = tei[ijkl]
            val = (k == l) ? 0.5 : 1
            val = (j == l) ? 0.5 : 1
            val3 = (j == k) ? 0.5 : 1
            F[i,j] += val * D[k,l] * (2*tei[ijkl] - tei[ikjl])
            F[i,j] += val * D[l,k] * (2*tei[ijkl] - tei[iljk])
            #F[i,j] += 4 * val * D[k,l] * (tei[ijkl])
            #F[i,k] -= val * D[j,l] * (tei[ijkl])
            #F[i,l] -= val * D[j,k] * (tei[ijkl])
            #F[i,j] += 2*D[k,l] * eri
            #F[i,k] -= D[j,l] * eri
            #F[k,i] = F[i,k]
            #F[l,i] = F[i,l]
        end
        F[j,i] = F[i,j]
        #F[i,j] += H[i,j]
    end
    return F
end

function twoei_six(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
    #F = deepcopy(H)
    F = zeros(7,7)
    J::Array{Float64,2} = zeros(7,7)
    K::Array{Float64,2} = zeros(7,7)
    for i::Int64 in 1:7
        for j::Int64 in 1:i
            for k::Int64 in 1:7
                for l::Int64 in 1:k
                    ij::Int64 = index(i,j)
                    kl::Int64 = index(k,l)
                    if (ij < kl) continue end
                    ijkl::Int64 = index(ij,kl)
                    ik::Int64 = index(i,k)
                    jl::Int64 = index(j,l)
                    ikjl::Int64 = index(ik,jl)
                    il::Int64 = index(i,k)
                    jk::Int64 = index(j,l)
                    iljk::Int64 = index(il,jk)
                    val = (k == l) ? 0.5 : 1
                    F[i,j] += val * D[k,l] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
                    F[k,l] += val * D[i,j] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
                    #val = (i == j) ? 0.5 : 1
                    #val *= (k == l) ? 0.5 : 1
                    #val *= (i == k && j == l) ? 0.5 : 1
                    #J[i,j] += val * D[k,l] * tei[ijkl]
                    #J[j,i] += val * D[k,l] * tei[ijkl]
                    #J[i,j] += val * D[l,k] * tei[ijkl]
                    #J[j,i] += val * D[l,k] * tei[ijkl]
                    #J[k,l] += val * D[i,j] * tei[ijkl]
                    #J[l,k] += val * D[i,j] * tei[ijkl]
                    #J[k,l] += val * D[j,i] * tei[ijkl]
                    #J[l,k] += val * D[j,i] * tei[ijkl]
                    #K[i,k] += val * D[j,l] * tei[ijkl]
                    #K[i,l] += val * D[j,k] * tei[ijkl]
                    #K[j,k] += val * D[i,l] * tei[ijkl]
                    #K[j,l] += val * D[i,k] * tei[ijkl]
                end
            end
            F[j,i] = F[i,j]
        end
    end
    F += H
    #F += 2*J - K
    return F
end
=#

#function compile()
#    precompile(hf.E, tuple(Array{String,1}))
#    precompile(hf.iteration, tuple(Array{Float64},Array{Float64},Array{Float64},Array{Float64}))
#    precompile(hf.twoei_base, tuple(Array{Float64},Array{Float64},Array{Float64},Array{Float64}))
#end

#we want to precompile all involved Cdules to reduce cold runs
#include("./snoop/precompile_Base.jl")
#_precompile_base()
#include("./snoop/precompile_Core.jl")
#_precompile_core()
#include("./snoop/precompile_hf.jl")
#_precompile_()
#include("./snoop/precompile_io.jl")
#_precompile_()
#include("./snoop/precompile_LinearAlgebra.jl")
#_precompile_()
#include("./snoop/precompile_openchem.jl")
#_precompile_()
#include("./snoop/precompile_unknown.jl")
#_precompile_()

#------------------------------#
#          JuliChem.jl         #
#------------------------------#
"""
     do_input(file::String)
Summary
======
Read in input file.

Arguments
======
file = name of input file to read in
"""
function do_input(file::String)
    #read in .inp and .dat files
    print("Reading in input files...")
    inpname::String = file
    inp::Array{String,1} = readin(inpname)
    datname::String = replace(inpname, ".inp"=>".dat")
    dat::Array{String,1} = readin(datname)
    println("Done!")
    println("")

    #do file processing
    #println("Processing input file...")
    coord::Matrix{Float64} = geomin(inp)

    return dat
end

"""
    do_scf(dat::Array{String,1})
Summary
======
Execute the JuliChem SCF algorithm.

Arguments
======
dat = input data file object
"""
function do_scf(dat::Array{String,1})
    #determine and perform proper method
    GC.enable(false)
    scf = energy(dat)
    GC.enable(true)
    GC.gc()

    return scf
end

#function properties(scf::Data)
#    properties(scf)
#end

"""
    do_exe(file::String)
Summary
======
Combined input file read-in + SCF execution.

Arguments
======
file = name of input file to read in
"""
function do_exe(file::String)
    println("----------------------------------------")
    println("Welcome to JuliChem!")
    println("JuliChem is a software package written")
    println("in Julia for the purpose of quantum")
    println("chemical calculations.")
    println("Let's get this party started!")
    println("----------------------------------------")
    println(" ")

    #read in input file
    inp = do_input(file)

    #perform scf calculation
    scf = do_scf(inp)

    #we have run to completion! :)
    println("----------------------------------------")
    println("The calculation has run to completion!")
    println("Sayonara!")
    println("----------------------------------------")

    return scf
end

#we want to precompile all involved modules to reduce cold runs
#include("./snoop/precompile_Base.jl")
#_precompile_base()
#include("./snoop/precompile_Core.jl")
#_precompile_core()
#include("./snoop/precompile_hf.jl")
#_precompile_()
#include("./snoop/precompile_io.jl")
#_precompile_()
#include("./snoop/precompile_LinearAlgebra.jl")
#_precompile_()
#include("./snoop/precompile_openchem.jl")
#_precompile_()
#include("./snoop/precompile_unknown.jl")
#_precompile_()

#function compile()
#    print("Precompiling OpenChem...")
#
#
#
#

#always do module recompilation upon running module block
#openchem.compile()

#------------------------------#
#           Script.jl          #
#------------------------------#
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    inp::String = ARGS[1]
    scf = exe(inp)
    return 0
end
