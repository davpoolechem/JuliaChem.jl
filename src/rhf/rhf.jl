#includes and Cdule imports
Base.include(@__MODULE__,"../input/datafile.jl")
Base.include(@__MODULE__,"../math.jl")

import LinearAlgebra
import SparseArrays
import LinearAlgebra.eigvecs
import LinearAlgebra.eigvals

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
function rhf_energy(dat::Array{String,1}, FLAGS::Flags)
    println("========================================")
    println("  RESTRICTED CLOSED-SHELL HARTREE-FOCK"  )
    println("========================================")
    println("")

    norb::Int64 = FLAGS.BASIS.NORB
    scf::Data = Data(zeros(norb,norb), zeros(norb,norb), zeros(norb,norb), 0)

    #Step #1: Nuclear Repulsion Energy
    E_nuc::Float64 = read_in_enuc(dat)

    #Step #2: One-Electron Integrals
    S::Array{Float64,2} = read_in_oei(dat,"OVR")
    T::Array{Float64,2} = read_in_oei(dat,"KEI")
    V::Array{Float64,2} = read_in_oei(dat,"NAI")
    H::Array{Float64,2} = T+V

    #Step #3: Two-Electron Integrals
    tei::Array{Float64,1} = read_in_tei(dat)

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

    println("----------------------------------------")
    println("        Starting RHF iterations...")
    println("----------------------------------------")
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
        F = twoei(F, D, tei, H, FLAGS)
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
    println("The SCF calculation has converged!")
    println("----------------------------------------")
    println("Total SCF Energy: ",E," h")
    println(" ")
    println("========================================")
    println("      END RESTRICTED CLOSED-SHELL"       )
    println("              HARTREE-FOCK"              )
    println("========================================")
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
    println("========================================")
    println("            SYSTEM PROPERTIES"           )
    println("========================================")
    println("")
    println("========================================")
    println("          END SYSTEM PROPERTIES"         )
    println("========================================")
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
function twoei(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
    H::Array{Float64,2}, FLAGS::Flags)
    ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:7*8))

    F = deepcopy(H)
    #F = zeros(7,7)
    coulomb::Array{Float64,2} = zeros(7,7)
    exchange::Array{Float64,2} = zeros(7,7)
    for μ::Int64 in 1:FLAGS.BASIS.NORB
        μμ::Int64 = ioff[μ]
        for ν::Int64 in 1:μ
            νν::Int64 = ioff[ν]
            μν::Int64 = μμ + ν
            for λ::Int64 in 1:FLAGS.BASIS.NORB
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
