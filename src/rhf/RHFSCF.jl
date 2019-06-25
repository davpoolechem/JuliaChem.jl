Base.include(@__MODULE__,"../math/math.jl")

Base.include(@__MODULE__,"ReadIn.jl")

using JCStructs

using MPI
using Base.Threads
using Distributed
using LinearAlgebra
using HDF5

function rhf_energy(FLAGS::RHF_Flags, basis::Basis, read_in::Dict{String,Any})
  if (FLAGS.SCF.PREC == "Float64")
    return rhf_kernel(FLAGS,basis,read_in,oneunit(Float64))
  elseif (FLAGS.SCF.PREC == "Float32")
    return rhf_kernel(FLAGS,basis,read_in,oneunit(Float32))
  end
end


"""
	 rhf_kernel(FLAGS::RHF_Flags, basis::Basis, read_in::Dict{String,Any},
       type::T)
Summary
======
Perform the core RHF SCF algorithm.

Arguments
======
FLAGS = Input flags

basis = Generated basis set

read_in = file required to read in from input file

type = Precision of variables in calculation
"""
function rhf_kernel(FLAGS::RHF_Flags, basis::Basis, read_in::Dict{String,Any},
  type::T) where {T<:AbstractFloat}

  norb::UInt32 = FLAGS.BASIS.NORB
  comm=MPI.COMM_WORLD

  json_debug::Any = ""
  if (FLAGS.SCF.DEBUG == true)
    json_debug = open(FLAGS.CTRL.NAME*"-debug.json","w")
  end

  #== read variables from input if needed ==#
  E_nuc::T = read_in["enuc"]

  S::Array{T,2} = read_in_oei(read_in["ovr"], FLAGS)
  H::Array{T,2} = read_in_oei(read_in["hcore"], FLAGS)

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    output_H = Dict([("Core Hamiltonian",H)])
    write(json_debug,JSON.json(output_H))
  end

  #== build the orthogonalization matrix ==#
  S_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(S))

  S_eval_diag::Array{T,1} = eigvals(LinearAlgebra.Hermitian(S))

  S_eval::Array{T,2} = zeros(norb,norb)
  for i::UInt32 in 1:norb
    S_eval[i,i] = S_eval_diag[i]
  end

  ortho::Array{T,2} = S_evec*
    (LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    output_ortho = Dict([("Orthogonalization Matrix",ortho)])
    write(json_debug,JSON.json(output_ortho))
  end

  #== build the initial matrices ==#
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

  D_old::Array{T,2} = deepcopy(D)

  E::T = E_elec + E_nuc
  E_old::T = E

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    output_F_initial = Dict([("Initial Fock Matrix",F)])
    output_D_initial = Dict([("Initial Density Matrix",D)])

    write(json_debug,JSON.json(output_F_initial))
    write(json_debug,JSON.json(output_D_initial))
  end

  if (MPI.Comm_rank(comm) == 0)
    println(0,"     ", E)
  end

  #=============================#
  #== start scf cycles: #7-10 ==#
  #=============================#
  converged::Bool = false
  iter::UInt32 = 1

  c = h5open("sto3g-h2.h5", "r") do tei
    while(!converged)
      #== build fock matrix ==#
	  F_temp = twoei(F, D, tei, H, FLAGS, basis)

	  F = MPI.Allreduce(F_temp,MPI.SUM,comm)
	  MPI.Barrier(comm)

	  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
        println("Skeleton Fock matrix:")
        display(F)
        println("")

	    output_iter_data = Dict([("SCF Iteration",iter),("Fock Matrix",F),
	      ("Density Matrix",D)])

	    write(json_debug,JSON.json(output_iter_data))
	  end

	  F += deepcopy(H)

      #== obtain new F,D,C matrices ==#
	  F, D, C, E_elec = iteration(F, D_old, H, ortho, FLAGS)
      F = transpose(ortho)*F*ortho

      #== check for convergence ==#
      ΔD::Array{T,2} = D - D_old
	  D_rms::T = √(∑(ΔD,ΔD))

	  E = E_elec+E_nuc
	  ΔE::T = E - E_old

	  if (MPI.Comm_rank(comm) == 0)
	    println(iter,"     ", E,"     ", ΔE,"     ", D_rms)
	  end

	  converged = (ΔE <= FLAGS.SCF.DELE) && (D_rms <= FLAGS.SCF.RMSD)
	  iter += 1
      if (iter > FLAGS.SCF.NITER) break end

      #== if not converged, replace old D and E values for next iteration ==#
      D_old = deepcopy(D)
      E_old = E
    end
  end

  if (iter > FLAGS.SCF.NITER)
    if (MPI.Comm_rank(comm) == 0)
	    println(" ")
      println("----------------------------------------")
      println(" The SCF calculation did not converge.  ")
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

  #== obtain new orbital coefficients ==#
  F_eval::Array{T,1} = eigvals(LinearAlgebra.Hermitian(F))

  F_evec::Array{T,2} = eigvecs(LinearAlgebra.Hermitian(F))
  F_evec = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

  C::Array{T,2} = ortho*F_evec

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    println("New orbitals:")
    display(C)
    println("")
  end

  #== build new density matrix ==#
  for i::UInt32 in 1:FLAGS.BASIS.NORB, j::UInt32 in 1:i
    D[i,j] = ∑(C[i,1:FLAGS.BASIS.NOCC],C[j,1:FLAGS.BASIS.NOCC])
    D[i,j] *= 2
    D[j,i] = D[i,j]
  end

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    println("New density matrix:")
    display(D)
    println("")
  end

  #== compute new SCF energy ==#
  EHF1::T = ∑(D,F)
  EHF2::T = ∑(D,H)
  E_elec::T = (EHF1 + EHF2)/2

  if (FLAGS.SCF.DEBUG == true && MPI.Comm_rank(comm) == 0)
    println("New energy:")
    display("$EHF1, $EHF2")
    println("")
  end

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
@inline function index(a::UInt32,b::UInt32)
  index::UInt32 = (a*(a-1)) >> 1 #bitwise divide by 2
  index += b
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

function twoei(F::Array{T,2}, D::Array{T,2}, tei::HDF5File,
  H::Array{T,2}, FLAGS::RHF_Flags, basis::Basis) where {T<:AbstractFloat}

  comm=MPI.COMM_WORLD
  norb::UInt32 = FLAGS.BASIS.NORB
  nsh::UInt32 = length(basis.shells)
  ioff::Array{UInt32,1} = map((x) -> x*(x-1)/2, collect(1:norb*(norb+1)))

  F = zeros(norb,norb)
  mutex = Base.Threads.Mutex()

  for bra_pairs::UInt32 in nsh*(nsh+1)/2:-1:1
    if(MPI.Comm_rank(comm) == bra_pairs%MPI.Comm_size(comm))
      ish::UInt32 = ceil(((-1+sqrt(1+8*bra_pairs))/2))
      jsh::UInt32 = bra_pairs - ioff[ish]

      if (ish < jsh) continue end

      ijsh::UInt32 = index(ish,jsh)

      Threads.@threads for ket_pairs::UInt32 in bra_pairs:-1:1
        ksh::UInt32 = ceil(((-1+sqrt(1+8*ket_pairs))/2))
        lsh::UInt32 = ket_pairs - ioff[ksh]

        if (ksh < lsh) continue end

        klsh::UInt32 = index(ksh,lsh)
        if (klsh > ijsh) continue end

		bra::ShPair = ShPair(basis.shells[ish], basis.shells[jsh])
		ket::ShPair = ShPair(basis.shells[ksh], basis.shells[lsh])
		quartet::ShQuartet = ShQuartet(bra,ket)

		eri_batch::Array{T,1} = shellquart(D, quartet, tei, mutex)

	    F_priv::Array{T,2} = zeros(norb,norb)
		if (max(eri_batch...) >= 1E-10)
          F_priv = dirfck(D, eri_batch, quartet)
        end

		lock(mutex)
		F += F_priv
		unlock(mutex)
      end
    end
  end

  for iorb::UInt32 in 1:norb, jorb::UInt32 in 1:iorb
    if (iorb != jorb)
      F[iorb,jorb] /= 2
      F[jorb,iorb] = F[iorb,jorb]
    end
  end

  return F
end

function shellquart(D::Array{T,2},quartet::ShQuartet,tei_file::HDF5File, mutex) where {T<:AbstractFloat}
  norb = size(D)[1]

  nμ = quartet.bra.sh_a.nbas
  nν = quartet.bra.sh_b.nbas
  nλ = quartet.ket.sh_a.nbas
  nσ = quartet.ket.sh_b.nbas

  pμ = quartet.bra.sh_a.pos
  pν = quartet.bra.sh_b.pos
  pλ = quartet.ket.sh_a.pos
  pσ = quartet.ket.sh_b.pos

  eri_batch::Array{T,1} = [ ]

  lock(mutex)
  tei_list::Array{Float64,1} = read(tei_file, "tei")
  unlock(mutex)

  for μ::UInt32 in pμ:pμ+(nμ-1), ν::UInt32 in pν:pν+(nν-1)
    μν = index(μ,ν)

	  for λ::UInt32 in pλ:pλ+(nλ-1), σ::UInt32 in pσ:pσ+(nσ-1)
	    λσ = index(λ,σ)
        μνλσ::UInt32 = index(μν,λσ)

        push!(eri_batch,tei_list[μνλσ])
    end
  end

  return deepcopy(eri_batch)
end

function dirfck(D::Array{T,2}, eri_batch::Array{T,1},
  quartet::ShQuartet) where {T<:AbstractFloat}

  norb = size(D)[1]

  F_priv::Array{T,2} = fill(0.0,(norb,norb))

  nμ = quartet.bra.sh_a.nbas
  nν = quartet.bra.sh_b.nbas
  nλ = quartet.ket.sh_a.nbas
  nσ = quartet.ket.sh_b.nbas

  pμ = quartet.bra.sh_a.pos
  pν = quartet.bra.sh_b.pos
  pλ = quartet.ket.sh_a.pos
  pσ = quartet.ket.sh_b.pos

  for μ::UInt32 in pμ:pμ+(nμ-1), ν::UInt32 in pν:pν+(nν-1)
    if (μ < ν) continue end

    μν = index(μ,ν)
	μν_idx::UInt32 = nν*nλ*nσ*(μ-pμ) + nλ*nσ*(ν-pν)

    for λ::UInt32 in pλ:pλ+(nλ-1), σ::UInt32 in pσ:pσ+(nσ-1)
	  if (λ < σ) continue end

	  λσ = index(λ,σ)
	  if (μν < λσ) continue end

	  μνλσ::UInt32 = μν_idx + nσ*(λ-pλ) + (σ-pσ) + 1

	  eri::T = eri_batch[μνλσ]
	  eri *= (μ == ν) ? 0.5 : 1.0
	  eri *= (λ == σ) ? 0.5 : 1.0
	  eri *= ((μ == λ) && (ν == σ)) ? 0.5 : 1.0

	  if (eri <= 1E-10) continue end

	  F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
	  F_priv[μ,ν] += 4.0 * D[λ,σ] * eri
      F_priv[μ,λ] -= D[ν,σ] * eri
	  F_priv[μ,σ] -= D[ν,λ] * eri
	  F_priv[max(ν,λ),min(ν,λ)] -= D[max(μ,σ),min(μ,σ)] * eri
	  F_priv[max(ν,σ),min(ν,σ)] -= D[max(μ,λ),min(μ,λ)] * eri

	  F_priv[σ,λ] = F_priv[λ,σ]
	  F_priv[ν,μ] = F_priv[μ,ν]
	  F_priv[λ,μ] = F_priv[μ,λ]
	  F_priv[σ,μ] = F_priv[μ,σ]
	  F_priv[min(λ,ν),max(λ,ν)] = F_priv[max(ν,λ),min(ν,λ)]
	  F_priv[min(σ,ν),max(σ,ν)] = F_priv[max(ν,σ),min(ν,σ)]
    end
  end

  return F_priv
end
