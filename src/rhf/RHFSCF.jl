Base.include(@__MODULE__,"../math/math.jl")

Base.include(@__MODULE__,"ReadIn.jl")

using JCStructs

using MPI
using Base.Threads
using Distributed
using LinearAlgebra

function rhf_energy(FLAGS::RHF_Flags, basis::Basis, read_in::Dict{String,Any})
    if (FLAGS.SCF.PREC == "Float64")
        return rhf_kernel(FLAGS,basis,read_in,oneunit(Float64))
    elseif (FLAGS.SCF.PREC == "Float32")
        return rhf_kernel(FLAGS,basis,read_in,oneunit(Float32))
    end
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
function rhf_kernel(FLAGS::RHF_Flags, basis::Basis, read_in::Dict{String,Any},
                        type::T) where {T<:AbstractFloat}
    norb::UInt32 = FLAGS.BASIS.NORB
    comm=MPI.COMM_WORLD

    json_debug::Any = ""
    if (FLAGS.SCF.DEBUG == true)
        json_debug = open(FLAGS.CTRL.NAME*"-debug.json","w")
    end

    #Step #1: Nuclear Repulsion Energy
    E_nuc::T = read_in["enuc"]

    #Step #2: One-Electron Integrals
    S::Array{T,2} = read_in_oei(read_in["ovr"], FLAGS)
    T_oei::Array{T,2} = read_in_oei(read_in["kei"], FLAGS)
    V::Array{T,2} = read_in_oei(read_in["nai"], FLAGS)
    H::Array{T,2} = T_oei+V

    if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
        output_H = Dict([("Core Hamiltonian",H)])
        write(json_debug,JSON.json(output_H))
    end

    #Step #3: Two-Electron Integrals
    tei::Array{T,1} = read_in_tei(read_in["tei"], FLAGS)

    #Step #4: Build the Orthogonalization Matrix
    S_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(S))

    S_eval_diag::Array{T,1} = eigvals(LinearAlgebra.Hermitian(S))

    S_eval::Array{T,2} = zeros(norb,norb)
    for i::UInt32 in 1:norb
        S_eval[i,i] = S_eval_diag[i]
    end

    ortho::Array{T,2} = S_evec*(LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)
    if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
        output_ortho = Dict([("Orthogonalization Matrix",ortho)])
        write(json_debug,JSON.json(output_ortho))
    end

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

    if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
        output_F_initial = Dict([("Initial Fock Matrix",F)])
        output_D_initial = Dict([("Initial Density Matrix",D)])

        write(json_debug,JSON.json(output_F_initial))
        write(json_debug,JSON.json(output_D_initial))
    end

    if (MPI.Comm_rank(comm) == 0)
        println(0,"     ", E)

    end

    #start scf cycles: #7-10
    converged::Bool = false
    iter::UInt32 = 1
    while(!converged)

        #multilevel MPI+threads parallel algorithm
        F_temp = twoei(F, D, tei, H, FLAGS, basis)

        F = MPI.Allreduce(F_temp,MPI.SUM,comm)
        MPI.Barrier(comm)

        F += deepcopy(H)

        #println("Initial Fock matrix:")
        if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
            output_iter_data = Dict([("SCF Iteration",iter),("Fock Matrix",F),
                                        ("Density Matrix",D)])

            write(json_debug,JSON.json(output_iter_data))
        end

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

        converged = (ΔE <= FLAGS.SCF.DELE) && (D_rms <= FLAGS.SCF.RMSD)
        iter += 1
        if (iter > FLAGS.SCF.NITER) break end
    end

    if (iter > FLAGS.SCF.NITER)
        if (MPI.Comm_rank(comm) == 0)
            println(" ")
            println("----------------------------------------")
            println("   The SCF calculation not converged.   ")
            println("      Restart data is being output.     ")
            println("----------------------------------------")
            println(" ")
        end

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

        if (FLAGS.SCF.DEBUG == true)
            close(json_debug)
        end

        return Data(F, D, C, E)
    end
end

#=
function rhf_energy(FLAGS::RHF_Flags, restart::RHFRestartData)
    norb::UInt32 = FLAGS.BASIS.NORB
    comm = MPI.COMM_WORLD

    H::Array{T,2} = T+V
    tei::Array{T,1} = read_in_tei()

    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("      Continuing RHF iterations...                ")
        println("----------------------------------------          ")
        println(" ")
        println("Iter      Energy                   ΔE                   Drms")
    end

    #start scf cycles: #7-10
    converged::Bool = false
    iter::UInt32 = restart.iter
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

        converged = (ΔE <= FLAGS.SCF.DELE) && (D_rms <= FLAGS.SCF.RMSD)
        iter += 1
        if (iter > FLAGS.SCF.NITER) break end
    end

    if (iter > FLAGS.SCF.NITER)
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
=#
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
    ortho::Array{T,2}, FLAGS::RHF_Flags) where {T<:AbstractFloat}

    #Step #8: Build the New Density Matrix
    F_eval::Array{T,1} = eigvals(LinearAlgebra.Hermitian(F))

    F_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(F))
    F_evec = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

    C::Array{T,2} = ortho*F_evec

    for i::UInt32 in 1:FLAGS.BASIS.NORB, j::UInt32 in 1:i
        D[i,j] = ∑(C[i,1:FLAGS.BASIS.NOCC],C[j,1:FLAGS.BASIS.NOCC])
        D[j,i] = D[i,j]
    end

    #Step #9: Compute the New SCF Energy
    E_elec::T = ∑(D,H + F)

    return (F, D, C, E_elec)
end
#=
"""
     index(a::UInt32,b::UInt32)
Summary
======
Triangular indexing determination.

Arguments
======
a = row index

b = column index
"""
=#
@inline function index(a::UInt32,b::UInt32,ioff::Array{UInt32,1})
    return ioff[a] + b
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
    H::Array{T,2}, FLAGS::RHF_Flags, basis::Basis) where {T<:AbstractFloat}

    comm=MPI.COMM_WORLD
    norb::UInt32 = FLAGS.BASIS.NORB
    nsh::UInt32 = length(basis.shells)
    ioff::Array{UInt32,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F = zeros(norb,norb)
    mutex = Base.Threads.Mutex()

    for bra_pairs::UInt32 in 1:ioff[nsh]
        if(MPI.Comm_rank(comm) == bra_pairs%MPI.Comm_size(comm))
            bra_sh_a::UInt32 = ceil(((-1+sqrt(1+8*bra_pairs))/2))
            bra_sh_b::UInt32 = bra_pairs%bra_sh_a + 1
            bra::ShPair = ShPair(basis.shells[bra_sh_a], basis.shells[bra_sh_b])

            Threads.@threads for ket_pairs::UInt32 in 1:bra_pairs
                ket_sh_a::UInt32 = ceil(((-1+sqrt(1+8*ket_pairs))/2))
                ket_sh_b::UInt32 = ket_pairs%ket_sh_a + 1

                ket::ShPair = ShPair(basis.shells[ket_sh_a], basis.shells[ket_sh_b])
                quartet::ShQuartet = ShQuartet(bra,ket)

                eri_batch::Array{T,1} = shellquart(D, tei, quartet)
                F_priv::Array{T,2} = dirfck(D, eri_batch, quartet)
                #F_priv::Array{T,2} = zeros(norb,norb)

                lock(mutex)
                #push!(debug_array,"$μ, $ν, $λ, $σ, $μν_idx, $λσ_idx")
                F += F_priv
                unlock(mutex)
            end
        end
    end

    #display(sort(debug_array))

    return F
end

function shellquart(D::Array{T,2}, tei::Array{T,1},quartet::ShQuartet) where {T<:AbstractFloat}
    norb = size(D)[1]
    ioff::Array{UInt32,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    nμ = quartet.bra.sh_a.nbas
    nν = quartet.bra.sh_b.nbas
    nλ = quartet.ket.sh_a.nbas
    nσ = quartet.ket.sh_b.nbas

    eri_batch::Array{T,1} = [ ]

    for μμ::UInt32 in 1:nμ
        μ::UInt32 = quartet.bra.sh_a.pos + (μμ-1)

        for νν::UInt32 in 1:nν
            ν::UInt32 = quartet.bra.sh_b.pos + (νν-1)
            μν = index(μ,ν,ioff)

            for λλ::UInt32 in 1:nλ
                λ::UInt32 = quartet.ket.sh_a.pos + (λλ-1)

                for σσ::UInt32 in 1:nσ
                    σ::UInt32 = quartet.ket.sh_b.pos + (σσ-1)

                    λσ = index(λ,σ,ioff)
                    μνλσ::UInt32 = index(μν,λσ,ioff)

                    push!(eri_batch,tei[μνλσ])
                end
            end
        end
    end
    return deepcopy(eri_batch)
end

function dirfck(D::Array{T,2}, eri_batch::Array{T,1},quartet::ShQuartet) where {T<:AbstractFloat}
    norb = size(D)[1]
    ioff::Array{UInt32,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

    F_priv::Array{T,2} = fill(0.0,(norb,norb))

    nμ = quartet.bra.sh_a.nbas
    nν = quartet.bra.sh_b.nbas
    nλ = quartet.ket.sh_a.nbas
    nσ = quartet.ket.sh_b.nbas

    for μμ::UInt32 in 1:nμ
        μ::UInt32 = quartet.bra.sh_a.pos + (μμ-1)
        μ_idx::UInt32 = nν*nλ*nσ*(μμ-1)

        for νν::UInt32 in 1:nν
            ν::UInt32 = quartet.bra.sh_b.pos + (νν-1)
            μν_idx::UInt32 = μ_idx + nλ*nσ*(νν-1)

            if (μ < ν) continue end

            for λλ::UInt32 in 1:nλ
                λ::UInt32 = quartet.ket.sh_a.pos + (λλ-1)
                μνλ_idx::UInt32 = μν_idx + nσ*(λλ-1)

                for σσ::UInt32 in 1:nσ
                    σ::UInt32 = quartet.ket.sh_b.pos + (σσ-1)
                    μνλσ::UInt32 = μνλ_idx + (σσ-1) + 1

                    if (λ < σ) continue end

                    #println("\"$μ, $ν, $λ, $σ\"")

                    val::T = (μ == ν) ? 0.5 : 1.0
                    val *= (λ == σ) ? 0.5 : 1.0
                    val *= ((μ == λ) && (ν == σ)) ? 0.5 : 1.0
                    eri::T = val * eri_batch[μνλσ]

                    if (eri <= 1E-10) continue end

                    F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
                    F_priv[μ,ν] += 4.0 * D[λ,σ] * eri
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
                end
            end
        end
    end
    return F_priv
end
