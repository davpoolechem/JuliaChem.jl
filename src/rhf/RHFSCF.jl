Base.include(@__MODULE__,"../math/math.jl")

include("../input/InputFunctions.jl")

using JCStructs

using MPI
using Base.Threads
using Distributed
using LinearAlgebra

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
function rhf_energy(FLAGS::Flags, basis::Basis)
    norb::Int64 = FLAGS.BASIS.NORB
    scf = Data(zeros(norb,norb), zeros(norb,norb), zeros(norb,norb), 0)
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
        F_temp = twoei(F, D, tei, H, FLAGS, basis)

        F = MPI.Allreduce(F_temp,MPI.SUM,comm)
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

    if (iter > FLAGS.HF.NITER)
        if (MPI.Comm_rank(comm) == 0)
            println(" ")
            println("----------------------------------------")
            println("   The SCF calculation not converged.   ")
            println("      Restart data is being output.     ")
            println("----------------------------------------")
            println(" ")
        end

        #restart = RHFRestartData(H, ortho, iter, F, D, C, E)

        return RHFRestartData(H, ortho, iter, F, D, C, E)
    else
        if (MPI.Comm_rank(comm) == 0)
            println(" ")
            println("----------------------------------------")
            println("   The SCF calculation has converged!   ")
            println("----------------------------------------")
            println("Total SCF Energy: ",E," h")
            println(" ")
        end

        #scf = Data(F, D, C, E)

        return Data(F, D, C, E)
    end
end

function rhf_energy(FLAGS::Flags, restart::RHFRestartData)
    norb::Int64 = FLAGS.BASIS.NORB
    comm = MPI.COMM_WORLD

    H::Array{Float64,2} = T+V
    tei::Array{Float64,1} = read_in_tei()

    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("      Continuing RHF iterations...                ")
        println("----------------------------------------          ")
        println(" ")
        println("Iter      Energy                   ΔE                   Drms")
    end

    #start scf cycles: #7-10
    converged::Bool = false
    iter::Int64 = restart.iter
    while(!converged)

        #multilevel MPI+threads parallel algorithm
        F_temp = twoei(F, D, tei, H, FLAGS)

        F = MPI.Allreduce(F_temp,MPI.SUM,comm)
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

    if (iter > FLAGS.HF.NITER)
        if (MPI.Comm_rank(comm) == 0)
            println(" ")
            println("----------------------------------------")
            println("   The SCF calculation not converged.   ")
            println("      Restart data is being output.     ")
            println("----------------------------------------")
            println(" ")
        end

        #restart = RHFRestartData(H, ortho, iter, F, D, C, E)

        return RHFRestartData(H, ortho, iter, F, D, C, E)
    else
        if (MPI.Comm_rank(comm) == 0)
            println(" ")
            println("----------------------------------------")
            println("   The SCF calculation has converged!   ")
            println("----------------------------------------")
            println("Total SCF Energy: ",E," h")
            println(" ")
        end

        #scf = Data(F, D, C, E)

        return Data(F, D, C, E)
    end
end

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
#@inline function index(a::Int64,b::Int64,ioff::Array{Int64,1})
#    index::Int64 = (a > b) ? ioff[a] + b : ioff[b] + a
#    return index
#end

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
    H::Array{Float64,2}, FLAGS::Flags, basis::Basis)

    comm=MPI.COMM_WORLD
    norb::Int64 = FLAGS.BASIS.NORB
    nsh::Int64 = length(basis.shells)
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = zeros(norb,norb)
    mutex = Base.Threads.Mutex()

    for bra_pairs::Int64 in 1:ioff[nsh]
        if(MPI.Comm_rank(comm) == bra_pairs%MPI.Comm_size(comm))
            bra_sh_a::Int64 = ceil(((-1+sqrt(1+8*bra_pairs))/2))
            bra_sh_b::Int64 = bra_pairs%bra_sh_a + 1
            bra::ShPair = ShPair(basis.shells[bra_sh_a], basis.shells[bra_sh_b])

            Threads.@threads for ket_pairs::Int64 in 1:ioff[nsh]
                ket_sh_a::Int64 = ceil(((-1+sqrt(1+8*ket_pairs))/2))
                ket_sh_b::Int64 = ket_pairs%ket_sh_a + 1

                ket::ShPair = ShPair(basis.shells[ket_sh_a], basis.shells[ket_sh_b])
                quartet::ShQuartet = ShQuartet(bra,ket)

                F_priv::Array{Float64,2} = dirfck(D, tei, quartet)

                lock(mutex)
                F += F_priv
                unlock(mutex)
            end
        end
        MPI.Barrier(comm)
    end
    return F
end

function dirfck(D::Array{Float64,2}, tei::Array{Float64,1},quartet::ShQuartet)

    norb = size(D)[1]
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F_priv::Array{Float64,2} = fill(0.0,(norb,norb))

    nμ = quartet.bra.sh_a.nbas
    nν = quartet.bra.sh_b.nbas
    nλ = quartet.ket.sh_a.nbas
    nσ = quartet.ket.sh_b.nbas

    for μμ::Int64 in 0:nμ-1
        μ::Int64 = quartet.bra.sh_a.pos + μμ
        μ_idx::Int64 = nν*nλ*nσ*μμ

        for νν::Int64 in 0:nν-1
            ν::Int64 = quartet.bra.sh_b.pos + νν
            μν_idx::Int64 = μ_idx + nλ*nσ*νν

            if (μ < ν) continue end

            for λλ::Int64 in 0:nλ-1
                λ::Int64 = quartet.ket.sh_a.pos + λλ
                μνλ_idx::Int64 = μν_idx + nσ*λλ

                for σσ::Int64 in 0:nσ-1
                    σ::Int64 = quartet.ket.sh_b.pos + σσ
                    μνλσ::Int64 = μνλ_idx + σσ

                    #println("\"$μ, $ν, $λ, $σ\"")

                    if (λ < σ) continue end

                    val::Float64 = (μ == ν) ? 0.5 : 1.0
                    val::Float64 *= (λ == σ) ? 0.5 : 1.0
                    eri::Float64 = val * tei[μνλσ]

                    if (eri <= 1E-10) continue end

                    F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
                    F_priv[σ,λ] += 4.0 * D[μ,ν] * eri

                    F_priv[μ,λ] -= D[ν,σ] * eri
                    F_priv[μ,σ] -= D[ν,λ] * eri
                    F_priv[ν,λ] -= D[μ,σ] * eri
                    F_priv[ν,σ] -= D[μ,λ] * eri
                end
            end
        end
    end
    return F_priv
end

#=
function dirfck(D::Array{Float64,2}, tei::Array{Float64,1},quartet::ShQuartet)

    norb = size(D)[1]
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F_priv::Array{Float64,2} = fill(0.0,(norb,norb))

    for μν_idx::Int64 in 1:quartet.bra.nbas2
        μ::Int64 = quartet.bra.sh_a.pos + ceil(μν_idx/quartet.bra.nbas2) - 1
        ν::Int64 = quartet.bra.sh_b.pos + μν_idx%quartet.bra.nbas2
        μν::Int64 = index(μ,ν,ioff)

        if (μ < ν) continue end

        for λσ_idx::Int64 in 1:quartet.ket.nbas2
            λ::Int64 = quartet.ket.sh_a.pos + ceil(λσ_idx/quartet.ket.nbas2) - 1
            σ::Int64 = quartet.ket.sh_b.pos + λσ_idx%quartet.ket.nbas2

            if (λ < σ) continue end

            #println("\"$μ, $ν, $λ, $σ\"")

            #μνλσ::Int64 = μν + (ket.sh_b.nbas*λ+σ)
            λσ::Int64 = index(λ,σ,ioff)
            μνλσ::Int64 = index(μν,λσ,ioff)

            val::Float64 = (μ == ν) ? 0.5 : 1.0
            val::Float64 *= (λ == σ) ? 0.5 : 1.0
            eri::Float64 = val * tei[μνλσ]

            if (eri <= 1E-10) continue end

            F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
            F_priv[σ,λ] += 4.0 * D[μ,ν] * eri

            F_priv[μ,λ] -= D[ν,σ] * eri
            F_priv[μ,σ] -= D[ν,λ] * eri
            F_priv[ν,λ] -= D[μ,σ] * eri
            F_priv[ν,σ] -= D[μ,λ] * eri
        end
    end
    return F_priv
end
=#

end
