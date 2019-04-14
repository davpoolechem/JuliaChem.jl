Base.include(@__MODULE__,"../math/math.jl")

include("../input/InputFunctions.jl")

using JCStructs

using MPI
using Base.Threads
using Distributed
using LinearAlgebra

function rhf_energy(FLAGS::Flags)
    scf::Data = rhf_kernel(FLAGS,oneunit(FLAGS.CTRL.PREC))
    return scf
end

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
function rhf_kernel(FLAGS::Flags, type::T) where {T<:Number}
    norb::Int32 = FLAGS.BASIS.NORB
    scf::Data{T} = Data(Matrix{T}(undef,norb,norb),
                        Matrix{T}(undef,norb,norb),
                        Matrix{T}(undef,norb,norb),
                        zero(T))
    comm=MPI.COMM_WORLD

    #Step #1: Nuclear Repulsion Energy
    E_nuc::T = read_in_enuc(type)
    #println(E_nuc)

    #Step #2: One-Electron Integrals
    S::Array{T,2} = read_in_ovr(type)
    T_oei::Array{T,2} = read_in_kei(type)
    V::Array{T,2} = read_in_nai(type)
    H::Array{T,2} = T_oei+V

    #Step #3: Two-Electron Integrals
    tei::Array{T,1} = read_in_tei(type)

    #Step #4: Build the Orthogonalization Matrix
    #println("Initial S matrix:")
    #display(S)
    #println("")
    S_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(S))

    S_eval_diag::Array{T,1} = eigvals(LinearAlgebra.Hermitian(S))
    #println("Initial S_evec matrix:")
    #display(S_evec)
    #println("")
    S_eval::Array{T,2} = zeros(norb,norb)
    for i::Int32 in 1:norb
        S_eval[i,i] = S_eval_diag[i]
    end

    ortho::Array{T,2} = S_evec*(LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)
    #ortho::Array{T,2} = S_evec*(LinearAlgebra.sqrt(LinearAlgebra.inv(S_eval)))*transpose(S_evec)

    #Step #5: Build the Initial (Guess) Density
    F::Array{T,2} = transpose(ortho)*H*ortho
    D = Matrix{T}(undef,norb,norb)
    C = Matrix{T}(undef,norb,norb)

    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("       Starting RHF iterations...                 ")
        println("----------------------------------------          ")
        println(" ")
        println("Iter      Energy                   ΔE                   Drms")
    end

    F, D, C, E_elec = iteration(F, D, H, ortho, FLAGS)
    E::T = E_elec + E_nuc

    if (MPI.Comm_rank(comm) == 0)
        println(0,"     ", E)
    end

    #start scf cycles: #7-10
    converged::Bool = false
    iter::Int32 = 1
    while(!converged)

        #multilevel MPI+threads parallel algorithm
        F_temp = twoei(F, D, tei, H, FLAGS)

        F = MPI.Allreduce(F_temp,MPI.SUM,comm)
        MPI.Barrier(comm)

        F += deepcopy(H)

        #task-based algorithm
        #F = twoei_tasked(F, D, tei, H, FLAGS)

        #println("Initial Fock matrix:")
        #display(F)
        #println("")

        #Step #8: Build the New Density Matrix
        D_old::Array{T,2} = deepcopy(D)
        E_old::T = E

        F = transpose(ortho)*F*ortho
        F, D, C, E_elec = iteration(F, D, H, ortho, FLAGS)
        E = E_elec+E_nuc

        #Step #10: Test for Convergence
        ΔE::T = E - E_old

        ΔD::Array{T,2} = D - D_old
        D_rms::T = √(∑(ΔD,ΔD))

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
     iteration(F::Array{T,2}, D::Array{T,2}, H::Array{T,2}, ortho::Array{T,2})
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
function iteration(F::Array{T,2}, D::Array{T,2}, H::Array{T,2},
    ortho::Array{T,2}, FLAGS::Flags) where {T<:Number}

    #Step #8: Build the New Density Matrix
    F_eval::Array{T,1} = eigvals(LinearAlgebra.Hermitian(F))

    F_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(F))
    F_evec = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

    C::Array{T,2} = ortho*F_evec

    for i::Int32 in 1:FLAGS.BASIS.NORB, j::Int32 in 1:i
        D[i,j] = ∑(C[i,1:FLAGS.BASIS.NOCC],C[j,1:FLAGS.BASIS.NOCC])
        D[j,i] = D[i,j]
    end

    #Step #9: Compute the New SCF Energy
    E_elec::T = ∑(D,H + F)

    return (F, D, C, E_elec)
end
#=
"""
     index(a::Int32,b::Int32)
Summary
======
Triangular indexing determination.

Arguments
======
a = row index

b = column index
"""
=#
@inline function index(a::Int32,b::Int32,ioff::Array{Int32,1})
    index::Int32 = (a > b) ? ioff[a] + b : ioff[b] + a
    return index
end

#=
"""
     twoei(F::Array{T}, D::Array{T}, tei::Array{T}, H::Array{T})
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
function twoei(F::Array{T,2}, D::Array{T,2}, tei::Array{T,1},
    H::Array{T,2}, FLAGS::Flags) where {T<:Number}

    comm=MPI.COMM_WORLD
    norb::Int32 = FLAGS.BASIS.NORB

    ioff::Array{Int32,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = zeros(norb,norb)
    mutex = Base.Threads.Mutex()

    for μν_idx::Int32 in 1:ioff[norb]
        if(MPI.Comm_rank(comm) == μν_idx%MPI.Comm_size(comm))
            μ::Int32 = ceil(((-1+sqrt(1+8*μν_idx))/2))
            ν::Int32 = μν_idx%μ + 1
            μν::Int32 = index(μ,ν,ioff)

            Threads.@threads for λσ_idx::Int32 in 1:ioff[norb]
                λ::Int32 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
                σ::Int32 = λσ_idx%λ + 1

                λσ::Int32 = index(λ,σ,ioff)
                μνλσ::Int32 = index(μν,λσ,ioff)

                val::T = (μ == ν) ? 0.5 : 1.0
                val::T *= (λ == σ) ? 0.5 : 1.0
                eri::T = val * tei[μνλσ]

                if (eri <= 1E-10) continue end

                F_priv::Array{T,2} = zeros(norb,norb)

                F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
                F_priv[σ,λ] += 4.0 * D[μ,ν] * eri

                F_priv[μ,λ] -= D[ν,σ] * eri
                F_priv[μ,σ] -= D[ν,λ] * eri
                F_priv[ν,λ] -= D[μ,σ] * eri
                F_priv[ν,σ] -= D[μ,λ] * eri

                lock(mutex)
                F += F_priv
                unlock(mutex)
            end
        end
        MPI.Barrier(comm)
    end
    return F
end
