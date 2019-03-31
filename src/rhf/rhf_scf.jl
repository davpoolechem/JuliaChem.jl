#includes and Cdule imports
Base.include(@__MODULE__,"../input/input_functions.jl")
Base.include(@__MODULE__,"../math.jl")

import LinearAlgebra
import SparseArrays
import Base.Threads
import LinearAlgebra.eigvecs
import LinearAlgebra.eigvals
import Distributed

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

"""
     rhf_energy(dat::Array{String,1})
Summary
======
Perform the core RHF SCF algorithm.

Arguments
======
dat = Input data file object
"""
function rhf_energy(FLAGS::Flags)
    norb::Int64 = FLAGS.BASIS.NORB
    scf::Data = Data(zeros(norb,norb), zeros(norb,norb), zeros(norb,norb), 0)

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

    #Step #5: Build the Initial (Guess) Density
    F::Array{Float64,2} = transpose(ortho)*H*ortho
    D::Array{Float64,2} = zeros(norb,norb)
    C::Array{Float64,2} = zeros(norb,norb)

    println("----------------------------------------          ")
    println("       Starting RHF iterations...                 ")
    println("----------------------------------------          ")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")

    F, D, C, E_elec = iteration(F, D, H, ortho, FLAGS)
    E::Float64 = E_elec + E_nuc

    println(0,"     ", E)

    #start scf cycles: #7-10
    converged::Bool = false
    iter::Int64 = 1
    while(!converged)

        #Step #7: Compute the New Fock Matrix
        F = twoei_kl_opt(F, D, tei, H, FLAGS)

        #F = deepcopy(H)
        #@sync Distributed.@distributed reducer for μ::Int64 in 1:norb
        #    F_task = twoei_distributed(F, D, tei, H, FLAGS, μ)
        #end
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

        println(iter,"     ", E,"     ", ΔE,"     ", D_rms)

        converged = (ΔE <= FLAGS.HF.DELE) && (D_rms <= FLAGS.HF.RMSD)
        iter += 1
        if (iter > FLAGS.HF.NITER) break end
    end

    println(" ")
    println("----------------------------------------")
    println("   The SCF calculation has converged!   ")
    println("----------------------------------------")
    println("Total SCF Energy: ",E," h")
    println(" ")

    scf.Fock = F
    scf.Density = D
    scf.Coeff = C
    scf.Energy = E

    return scf
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
@inline function index(a::Int64,b::Int64,ioff::Array{Int64,1})
    index::Int64 = (a > b) ? ioff[a] + b : ioff[b] + a
    return index
end

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
function twoei_distributed(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags, μ::Int64)

    norb::Int64 = FLAGS.BASIS.NORB

    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = zeros(norb,norb)
    J::Array{Float64,2} = zeros(norb,norb)
    K::Array{Float64,2} = zeros(norb,norb)

    for ν::Int64 in 1:μ
        μν::Int64 = index(μ,ν,ioff)

        for λσ_idx::Int64 in 1:ioff[norb]
            λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
            σ::Int64 = λσ_idx%λ + 1

            λσ::Int64 = index(λ,σ,ioff)
            μνλσ::Int64 = index(μν,λσ,ioff)

            val::Float64 = (μ == ν) ? 0.5 : 1.0
            val::Float64 *= (λ == σ) ? 0.5 : 1.0

            J[λ,σ] += 2.0 * val * D[μ,ν] * tei[μνλσ]
            J[σ,λ] += 2.0 * val * D[μ,ν] * tei[μνλσ]
            eri::Float64 = val * tei[μνλσ]

            K[μ,λ] += val * D[ν,σ] * tei[μνλσ]
            K[μ,σ] += val * D[ν,λ] * tei[μνλσ]
            K[ν,λ] += val * D[μ,σ] * tei[μνλσ]
            K[ν,σ] += val * D[μ,λ] * tei[μνλσ]
        end
    end

    F = 2*J - K
    return F
end

"""
     twoei_threaded(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step, using multithreading.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""
function twoei_threaded(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags)

    norb::Int64 = FLAGS.BASIS.NORB

    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = deepcopy(H)
    J::Array{Float64,2} = zeros(norb,norb)
    K::Array{Float64,2} = zeros(norb,norb)
    #K::Array{Threads.Atomic{Float64},2} = fill(Threads.Atomic{Float64}(0.0),(norb,norb))

    #Threads.@threads for μν_idx::Int64 in 1:ioff[norb]
    for μν_idx::Int64 in 1:ioff[norb]
        μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
        ν::Int64 = μν_idx%μ + 1

        μν::Int64 = index(μ,ν,ioff)

        for λσ_idx::Int64 in 1:ioff[norb]
            λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
            σ::Int64 = λσ_idx%λ + 1

            λσ::Int64 = index(λ,σ,ioff)
            μνλσ::Int64 = index(μν,λσ,ioff)

            val::Float64 = (μ == ν) ? 0.5 : 1.0
            val::Float64 *= (λ == σ) ? 0.5 : 1.0
            eri::Float64 = val * tei[μνλσ]

            J[λ,σ] += 2.0 * D[μ,ν] * eri
            J[σ,λ] += 2.0 * D[μ,ν] * eri

            @async begin
                K[μ,λ] += D[ν,σ] * eri
                K[μ,σ] += D[ν,λ] * eri
                K[ν,λ] += D[μ,σ] * eri
                K[ν,σ] += D[μ,λ] * eri
            end
            #Threads.atomic_add!(K[μ,λ], D[ν,σ] * eri)
            #Threads.atomic_add!(K[μ,σ], D[ν,λ] * eri)
            #Threads.atomic_add!(K[ν,λ], D[μ,σ] * eri)
            #Threads.atomic_add!(K[ν,σ], D[μ,λ] * eri)
        end
    end

    F += 2*J - K
    #for i in 1:norb, j in 1:norb
    #    F[i,j] = 2*J[i,j] - K[i,j][]
    #end
    return F
end

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
function twoei_multilevel(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags, μ::Int64)

    norb::Int64 = FLAGS.BASIS.NORB
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    #F = zeros(norb,norb)
    μμ::Int64 = ioff[μ]

    coulomb::Array{Threads.Atomic{Float64},2} = fill(Threads.Atomic{Float64}(0.0),(norb,norb))
    exchange::Array{Threads.Atomic{Float64},2} = fill(Threads.Atomic{Float64}(0.0),(norb,norb))

    for ν::Int64 in 1:μ
        νν::Int64 = ioff[ν]
        μν::Int64 = μμ + ν

        Threads.@threads for λ::Int64 in 1:norb
            λλ::Int64 = ioff[λ]
            μλ::Int64 = (μ > λ) ? μμ + λ : μ + λλ
            νλ::Int64 = (ν > λ) ? νν + λ : ν + λλ

            for σ::Int64 in 1:λ
                σσ::Int64 = ioff[σ]
                μσ::Int64 = (μ > σ) ? μμ + σ : μ + σσ
                νσ::Int64 = (ν > σ) ? νν + σ : ν + σσ
                λσ::Int64 = (λ > σ) ? λλ + σ : λ + σσ
                μνλσ::Int64 = index(μν,λσ,ioff)
                μλνσ::Int64 = index(μλ,νσ,ioff)
                μσνλ::Int64 = index(μσ,νλ,ioff)
                #eri = tei[ijkl]
                #eri1 = 2*tei[ijkl] - tei[ikjl]
                #eri2 = 2*tei[ijkl] - tei[iljk]
                #eri = eri1 + eri2
                val::Float64 = (λ == σ) ? 0.5 : 1.0
                #if (k == l)
            #        eri1 *= 0.5
        #            eri2 *= 0.5
        #        end
                #F[i,j] += D[k,l] * eri1
                #F[i,j] += D[l,k] * eri2
                Threads.atomic_add!(coulomb[μ,ν], val*D[λ,σ]*tei[μνλσ])
                Threads.atomic_add!(exchange[μ,ν], val*D[λ,σ]*tei[μλνσ])
                Threads.atomic_add!(coulomb[μ,ν], val*D[σ,λ]*tei[μνλσ])
                Threads.atomic_add!(exchange[μ,ν], val*D[σ,λ]*tei[μσνλ])
            end
        end
        coulomb[ν,μ] = coulomb[μ,ν]
        exchange[ν,μ] = exchange[μ,ν]
        #F[j,i] = F[i,j]
    end
    for i in 1:norb, j in 1:norb
        F[i,j] += 2*coulomb[i,j] - exchange[i,j]
    end
    return F
end

"""
     twoei_threaded(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step, using multithreading.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""
function twoei_tasks(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags)

    norb::Int64 = FLAGS.BASIS.NORB

    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = deepcopy(H)
    J::Array{Float64,2} = zeros(norb,norb)
    K::Array{Float64,2} = zeros(norb,norb)
    #K::Array{Threads.Atomic{Float64},2} = fill(Threads.Atomic{Float64}(0.0),(norb,norb))

    #Threads.@threads for μν_idx::Int64 in 1:ioff[norb]
    @sync @async for μν_idx::Int64 in 1:ioff[norb]
        μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
        ν::Int64 = μν_idx%μ + 1

        μν::Int64 = index(μ,ν,ioff)

        for λσ_idx::Int64 in 1:ioff[norb]
            λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
            σ::Int64 = λσ_idx%λ + 1

            λσ::Int64 = index(λ,σ,ioff)
            μνλσ::Int64 = index(μν,λσ,ioff)

            val::Float64 = (μ == ν) ? 0.5 : 1.0
            val::Float64 *= (λ == σ) ? 0.5 : 1.0
            eri::Float64 = val * tei[μνλσ]

            J[λ,σ] += 2.0 * D[μ,ν] * eri
            J[σ,λ] += 2.0 * D[μ,ν] * eri

            @async begin
                K[μ,λ] += D[ν,σ] * eri
                K[μ,σ] += D[ν,λ] * eri
                K[ν,λ] += D[μ,σ] * eri
                K[ν,σ] += D[μ,λ] * eri
            end
            #Threads.atomic_add!(K[μ,λ], D[ν,σ] * eri)
            #Threads.atomic_add!(K[μ,σ], D[ν,λ] * eri)
            #Threads.atomic_add!(K[ν,λ], D[μ,σ] * eri)
            #Threads.atomic_add!(K[ν,σ], D[μ,λ] * eri)
        end
    end

    F += 2*J - K
    #for i in 1:norb, j in 1:norb
    #    F[i,j] = 2*J[i,j] - K[i,j][]
    #end
    return F
end

function twoei_kl_opt(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64},
    H::Array{Float64,2}, FLAGS::Flags)

    norb::Int64 = FLAGS.BASIS.NORB
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = zeros(norb,norb)

    J::Array{Float64,2} = zeros(norb,norb)
    K::Array{Float64,2} = zeros(norb,norb)

    for i::Int64 in 1:norb
        for j::Int64 in 1:i
            ij::Int64 = index(i,j,ioff)

            for k::Int64 in 1:norb
                for l::Int64 in 1:k

                    kl::Int64 = index(k,l,ioff)
                    ijkl::Int64 = index(ij,kl,ioff)

                    val::Float64 = (i == j) ? 0.5 : 1.0
                    val::Float64 *= (k == l) ? 0.5 : 1.0

                    J[k,l] += 2.0 * val * D[i,j] * tei[ijkl]
                    J[l,k] += 2.0 * val * D[i,j] * tei[ijkl]

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
