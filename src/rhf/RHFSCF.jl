module RHFSCF

Base.include(@__MODULE__,"../math/math.jl")

using RHFStructs

using InputFunctions
using InputStructs

import MPI
import Base.Threads
import Distributed
import LinearAlgebra
import LinearAlgebra.eigvecs
import LinearAlgebra.eigvals
import SparseArrays

#=
"""
     rhf_energy(dat::Array{String,1})
Summary
======
Perform the core RHF SCF algorithm.

Arguments
======
dat = Input data file object
"""
=#
function rhf_energy(FLAGS::Flags)
    norb::Int64 = FLAGS.BASIS.NORB
    scf::Data = Data(zeros(norb,norb), zeros(norb,norb), zeros(norb,norb), 0)
    comm=MPI.COMM_WORLD

    #Step #1: Nuclear Repulsion Energy
    E_nuc::Float64 = read_in_enuc()
    #println(E_nuc)

    #Step #2: One-Electron Integrals
    S::Array{Float64,2} = read_in_ovr()
    T::Array{Float64,2} = read_in_kei()
    V::Array{Float64,2} = read_in_nai()
    H::Array{Float64,2} = T+V

    #Step #3: Two-Electron Integrals
    tei::Array{Float64,1} = read_in_tei()

    #Step #4: Build the Orthogonalization Matrix
    #println("Initial S matrix:")
    #display(S)
    #println("")
    S_evec::Array{Float64,2} = eigvecs(LinearAlgebra.Hermitian(S))

    S_eval_diag::Array{Float64,1} = eigvals(LinearAlgebra.Hermitian(S))
    #println("Initial S_evec matrix:")
    #display(S_evec)
    #println("")
    S_eval::Array{Float64,2} = zeros(norb,norb)
    for i::Int64 in 1:norb
        S_eval[i,i] = S_eval_diag[i]
    end

    ortho::Array{Float64,2} = S_evec*(LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)
    #ortho::Array{Float64,2} = S_evec*(LinearAlgebra.sqrt(LinearAlgebra.inv(S_eval)))*transpose(S_evec)

    #Step #5: Build the Initial (Guess) Density
    F::Array{Float64,2} = transpose(ortho)*H*ortho
    D::Array{Float64,2} = zeros(norb,norb)
    C::Array{Float64,2} = zeros(norb,norb)

    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("       Starting RHF iterations...                 ")
        println("----------------------------------------          ")
        println(" ")
        println("Iter      Energy                   ΔE                   Drms")
    end

    F, D, C, E_elec = iteration(F, D, H, ortho, FLAGS)
    E::Float64 = E_elec + E_nuc

    if (MPI.Comm_rank(comm) == 0)
        println(0,"     ", E)
    end

    #start scf cycles: #7-10
    converged::Bool = false
    iter::Int64 = 1
    while(!converged)
        #multilevel MPI+threads parallel algorithm
        F_local = zeros(norb,norb)
        for μν::Int64 in 1:((norb*(norb+1))/2)
            if(MPI.Comm_rank(comm) == μν%MPI.Comm_size(comm))
                F_local += twoei(F, D, tei, H, FLAGS, μν)
            end
        end
        MPI.Barrier(comm)
            F = MPI.Allreduce(F_local,MPI.SUM,comm)
        MPI.Barrier(comm)

        F += deepcopy(H)

        #println("Initial Fock matrix:")
        #display(F)
        #println("")

        #Step #8: Build the New Density Matrix
        D_old::Array{Float64,2} = deepcopy(D)
        E_old::Float64 = E

        F = transpose(ortho)*F*ortho
        F, D, C, E_elec = iteration(F, D, H, ortho, FLAGS)
        E = E_elec+E_nuc

        #Step #10: Test for Convergence
        ΔE::Float64 = E - E_old

        ΔD::Array{Float64,2} = D - D_old
        D_rms::Float64 = √(∑(ΔD,ΔD))

        if (MPI.Comm_rank(comm) == 0)
            println(iter,"     ", E,"     ", ΔE,"     ", D_rms)
        end

        converged = (ΔE <= FLAGS.HF.DELE) && (D_rms <= FLAGS.HF.RMSD)
        iter += 1
        if (iter > FLAGS.HF.NITER) break end
    end

    if (MPI.Comm_rank(comm) == 0)
        println(" ")
        println("----------------------------------------")
        println("   The SCF calculation has converged!   ")
        println("----------------------------------------")
        println("Total SCF Energy: ",E," h")
        println(" ")
    end

    scf.Fock = F
    scf.Density = D
    scf.Coeff = C
    scf.Energy = E

    return scf
end
export rhf_energy

#=
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
=#
function iteration(F::Array{Float64,2}, D::Array{Float64,2}, H::Array{Float64,2},
    ortho::Array{Float64,2}, FLAGS::Flags)

    #Step #8: Build the New Density Matrix
    F_eval::Array{Float64,1} = eigvals(LinearAlgebra.Hermitian(F))

    F_evec::Array{Float64,2} = eigvecs(LinearAlgebra.Hermitian(F))
    F_evec = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

    C::Array{Float64,2} = ortho*F_evec

    for i::Int64 in 1:FLAGS.BASIS.NORB, j::Int64 in 1:i
        D[i,j] = ∑(C[i,1:FLAGS.BASIS.NOCC],C[j,1:FLAGS.BASIS.NOCC])
        D[j,i] = D[i,j]
    end

    #Step #9: Compute the New SCF Energy
    E_elec::Float64 = ∑(D,H + F)

    return (F, D, C, E_elec)
end
#=
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
=#
@inline function index(a::Int64,b::Int64,ioff::Array{Int64,1})
    #@assert (a >= b)
    index::Int64 = ioff[a] + b
    return index
end

#=
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
=#
function twoei(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags, μν_idx::Int64)

    debug_array = []

    norb::Int64 = FLAGS.BASIS.NORB
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))
    F = zeros(norb,norb)
    mutex = Base.Threads.Mutex()

    μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
    ν::Int64 = μν_idx%μ + 1
    μ, ν = (μ >= ν) ? (μ, ν) : (ν, μ)

    Threads.@threads for λσ_idx::Int64 in 1:μν_idx
        λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
        σ::Int64 = λσ_idx%λ + 1
        λ, σ = (λ >= σ) ? (λ, σ) : (σ, λ)

        if (μ < λ)
            μ, λ = λ, μ
            ν, σ = σ, ν
        elseif ((μ == λ) && (ν < σ))
            ν, σ = σ, ν
        end

        μν::Int64 = index(μ,ν,ioff)
        λσ::Int64 = index(λ,σ,ioff)

        μνλσ::Int64 = index(μν,λσ,ioff)

        val::Float64 = (μ == ν) ? 0.5 : 1.0
        val *= (λ == σ) ? 0.5 : 1.0
        val *= ((μ == λ) && (ν == σ)) ? 0.5 : 1.0
        eri::Float64 = val * tei[μνλσ]

        #if (eri <= 1E-10) continue end

        F_priv::Array{Float64,2} = zeros(norb,norb)

        F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
        F_priv[μ,ν] += 4.0 * D[σ,λ] * eri
        F_priv[μ,λ] -= D[ν,σ] * eri
        F_priv[μ,σ] -= D[ν,λ] * eri
        F_priv[ν,λ] -= D[μ,σ] * eri
        F_priv[ν,σ] -= D[μ,λ] * eri

        F_priv[σ,λ] = F_priv[λ,σ]
        F_priv[ν,μ] = F_priv[μ,ν]
        F_priv[λ,μ] = F_priv[μ,λ]
        F_priv[σ,μ] = F_priv[μ,σ]
        F_priv[λ,ν] = F_priv[ν,λ]
        F_priv[σ,ν] = F_priv[ν,σ]

        lock(mutex)
        push!(debug_array,"$μ, $ν, $λ, $σ, $μν_idx, $λσ_idx")
        F += F_priv
        unlock(mutex)
    end

    display(sort(debug_array))

    return F
end

end
