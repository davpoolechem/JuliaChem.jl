using MATH
using JCModules.Globals

using MPI
#using Base.Threads
#using Distributed
using LinearAlgebra
#using JLD
using HDF5
using PrettyTables

const do_continue_print = false

function rhf_energy(basis::BasisStructs.Basis,
  molecule::Union{Dict{String,Any},Dict{Any,Any}},
  scf_flags)

  return rhf_kernel(basis,molecule,scf_flags)
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
function rhf_kernel(basis::BasisStructs.Basis,
  molecule::Union{Dict{String,Any},Dict{Any,Any}}, scf_flags)

  comm=MPI.COMM_WORLD
  calculation_status = Dict([])

  #== read in some variables from scf input ==#
  debug::Bool = scf_flags["debug"]
  debug_output = debug ? h5open("debug.h5","w") : nothing

  niter::Int = scf_flags["niter"]

  #== read variables from input if needed ==#
  E_nuc::Float64 = molecule["enuc"]

  S = read_in_oei(molecule["ovr"], basis.norb)
  H = read_in_oei(molecule["hcore"], basis.norb)

  if debug && MPI.Comm_rank(comm) == 0
    #println("Overlap matrix:")
    #for shell_group in 0:cld(size(S)[1],5)
    #  ending = min((5*shell_group + 5), size(S)[2])
    #  S_debug = S[:,(5*shell_group + 1):ending]
    #  pretty_table(hcat(collect(1:1:size(S)[1]),S_debug),
    #    vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
    #    formatter = ft_printf("%5.6f", collect(2:1:6)),
    #    highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    #end
    h5write("debug.h5","SCF/0/S", S)
    #println("")


    #println("Hamiltonian matrix:")
    #for shell_group in 0:cld(size(H)[1],5)
    #  ending = min((5*shell_group + 5), size(H)[2])
    #  H_debug = H[:,(5*shell_group + 1):ending]
    #  pretty_table(hcat(collect(1:1:size(H)[1]),H_debug),
    #    vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
    #    formatter = ft_printf("%5.6f", collect(2:1:6)),
    #    highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    #end
    h5write("debug.h5","SCF/0/H", H)
    #println("")
  end

  #== build the orthogonalization matrix ==#
  S_evec = eigvecs(LinearAlgebra.Hermitian(S))

  S_eval_diag = eigvals(LinearAlgebra.Hermitian(S))

  S_eval = zeros(basis.norb,basis.norb)
  for i in 1:basis.norb
    S_eval[i,i] = S_eval_diag[i]
  end

  ortho = similar(H)
  ortho .= S_evec*(LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)

  if debug && MPI.Comm_rank(comm) == 0
  #  println("Ortho matrix:")
  #  for shell_group in 0:cld(size(ortho)[1],5)
  #    ending = min((5*shell_group + 5), size(ortho)[2])
  #    ortho_debug = ortho[:,(5*shell_group + 1):ending]
  #    pretty_table(hcat(collect(1:1:size(ortho)[1]),ortho_debug),
  #      vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
  #      formatter = ft_printf("%5.6f", collect(2:1:6)),
  #      highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
  #  end
    h5write("debug.h5","SCF/0/Ortho", ortho)
  #  println("")
  end

  #== build the initial matrices ==#
  F = H
  F_eval = Vector{Float64}(undef,basis.norb)
  F_evec = similar(F)
  F_mo = similar(F)

  D = similar(F)
  C = similar(F)

  if MPI.Comm_rank(comm) == 0
    println("----------------------------------------          ")
    println("       Starting RHF iterations...                 ")
    println("----------------------------------------          ")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")
  end

  E_elec = 0.0
  E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo, ortho, basis, 0,
    debug)
  F = deepcopy(F_mo)
  #indices_tocopy::CartesianIndices = CartesianIndices((1:size(F,1),
  #  1:size(F,2)))
  #copyto!(F, indices_tocopy, F_mo, indices_tocopy)

  E = E_elec + E_nuc
  E_old = E

  if MPI.Comm_rank(comm) == 0
    println(0,"     ", E)
  end

  #=============================#
  #== start scf cycles: #7-10 ==#
  #=============================#
  F, D, C, E, converged = scf_cycles(F, D, C, E, H, ortho, S, F_eval,
    F_evec, F_mo, E_nuc, E_elec, E_old, basis, scf_flags; debug=debug,
    niter=niter)

  if !converged
    if MPI.Comm_rank(comm) == 0
      println(" ")
      println("----------------------------------------")
      println(" The SCF calculation did not converge.  ")
      println("      Restart data is being output.     ")
      println("----------------------------------------")
      println(" ")
    end

    calculation_fail = Dict(
    "success" => false,
    "error" => Dict(
      "error_type" => "convergence_error",
      "error_message" => " SCF calculation did not converge within $niter
        iterations. "
      )
    )

    merge!(calculation_status, calculation_fail)

  else
    if MPI.Comm_rank(comm) == 0
      println(" ")
      println("----------------------------------------")
      println("   The SCF calculation has converged!   ")
      println("----------------------------------------")
      println("Total SCF Energy: ",E," h")
      println(" ")

      calculation_success = Dict(
      "return_result" => E,
      "success" => true,
      "properties" => Dict(
        "return_energy" => E,
        "nuclear_repulsion_energy" => E_nuc,
        #"scf_iterations" => iter,
        "scf_total_energy" => E
        )
      )

      merge!(calculation_status, calculation_success)
    end
  end

  if debug close(debug_output) end

  return F, D, C, E, calculation_status
end

function scf_cycles(F::Matrix{Float64}, D::Matrix{Float64}, C::Matrix{Float64},
  E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64}, S::Matrix{Float64},
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, F_mo::Matrix{Float64},
  E_nuc::Float64, E_elec::Float64, E_old::Float64, basis::BasisStructs.Basis,
  scf_flags; debug, niter)

  #== read in some more variables from scf flags input ==#
  ndiis::Int = scf_flags["ndiis"]
  dele::Float64 = scf_flags["dele"]
  rmsd::Float64 = scf_flags["rmsd"]

  #== build DIIS arrays ==#
  F_array = fill(similar(F), ndiis)

  e = similar(F)
  test_e = [ similar(F) ]
  e_array = fill(similar(F), ndiis)
  e_array_old = fill(similar(F), ndiis)
  
  test_F = [ similar(F) ]
  F_array_old = fill(similar(F), ndiis)

  #== build arrays needed for post-fock build iteration calculations ==#
  F_temp = similar(F)
  D_old = similar(F)
  ΔD = similar(F)

  #== build arrays needed for dynamic damping ==#
  damp_values = [ 0.25, 0.75 ]
  D_damp = [ similar(F) for i in 1:2 ]
  D_damp_rms = [ zero(Float64), zero(Float64) ]

  #== build variables needed for eri batching ==#
  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3

  quartet_batch_num_old = trunc(Integer,nindices/QUARTET_BATCH_SIZE) + 1

  #== build eri batch arrays ==#
  #eri_sizes::Vector{Int64} = load("tei_batch.jld",
  #  "Sizes/$quartet_batch_num_old")
  #length_eri_sizes::Int64 = length(eri_sizes)

  #@views eri_starts::Vector{Int64} = [1, [ sum(eri_sizes[1:i])+1 for i in 1:(length_eri_sizes-1)]... ]

  #eri_batch::Vector{Float64} = load("tei_batch.jld",
  #  "Integrals/$quartet_batch_num_old")

  #eri_sizes = []
  #eri_starts = []
  #eri_batch = []

  #== execute convergence procedure ==#
  scf_converged = true

  E = scf_cycles_kernel(F, D, C, E, H, ortho, S, E_nuc,
    E_elec, E_old, basis, F_array, e, e_array, e_array_old,
    F_array_old, F_temp, F_eval, F_evec, F_mo, D_old, ΔD, damp_values, D_damp,
    D_damp_rms, scf_converged, quartet_batch_num_old, test_e, test_F; 
    debug=debug, niter=niter, ndiis=ndiis, dele=dele, rmsd=rmsd)

  #== we are done! ==#
  if debug
    h5write("debug.h5","SCF/Final/F", F)
    h5write("debug.h5","SCF/Final/D", D)
    h5write("debug.h5","SCF/Final/C", C)
    h5write("debug.h5","SCF/Final/E", E)
    h5write("debug.h5","SCF/Final/converged", scf_converged)
  end

  return F, D, C, E, scf_converged
end

function scf_cycles_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  C::Matrix{Float64}, E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64},
  S::Matrix{Float64}, E_nuc::Float64, E_elec::Float64, E_old::Float64,
  basis::BasisStructs.Basis,
  F_array::Vector{Matrix{Float64}}, e::Matrix{Float64},
  e_array::Vector{Matrix{Float64}}, e_array_old::Vector{Matrix{Float64}},
  F_array_old::Vector{Matrix{Float64}}, F_temp::Matrix{Float64},
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, F_mo::Matrix{Float64},
  D_old::Matrix{Float64}, ΔD::Matrix{Float64}, damp_values::Vector{Float64},
  D_damp::Vector{Matrix{Float64}}, D_damp_rms::Vector{Float64},
  scf_converged::Bool, quartet_batch_num_old::Int64, 
  test_e::Vector{Matrix{Float64}}, test_F::Vector{Matrix{Float64}}; debug, 
  niter, ndiis, dele, rmsd)

  #== initialize a few more variables ==#
  comm=MPI.COMM_WORLD

  B_dim = 1
  #length_eri_sizes = length(eri_sizes)

  #=================================#
  #== now we start scf iterations ==#
  #=================================#
  iter = 1
  iter_converged = false

  while !iter_converged
    #== reset eri arrays ==#
    #if quartet_batch_num_old != 1 && iter != 1
    #  resize!(eri_sizes,length_eri_sizes)
    #  resize!(eri_starts,length_eri_sizes)

    #  eri_sizes[:] = load("tei_batch.jld",
  #      "Sizes/$quartet_batch_num_old")

    #  @views eri_starts[:] = [1, [ sum(eri_sizes[1:i])+1 for i in 1:(length_eri_sizes-1)]... ]
      #eri_starts[:] = load("tei_batch.jld",
      #  "Starts/$quartet_batch_num_old")
      #@views eri_starts[:] = eri_starts[:] .- (eri_starts[1] - 1)

    #  resize!(eri_batch,sum(eri_sizes))
    #  eri_batch[:] = load("tei_batch.jld","Integrals/$quartet_batch_num_old")
    #end

    #== build fock matrix ==#
    F_temp .= twoei(F, D, H, basis; debug=debug)

    F .= MPI.Allreduce(F_temp,MPI.SUM,comm)
    MPI.Barrier(comm)

    if debug && MPI.Comm_rank(comm) == 0
    #  println("Skeleton Fock matrix:")
    #  for shell_group in 0:cld(size(F)[1],5)
    #    ending = min((5*shell_group + 5), size(F)[2])
    #    F_debug = F[:,(5*shell_group + 1):ending]
    #    pretty_table(hcat(collect(1:1:size(F)[1]),F_debug),
    #      vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
    #      formatter = ft_printf("%5.6f", collect(2:1:6)),
    #      highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    #  end
      h5write("debug.h5","SCF/$iter/F/Skeleton", F)
    #  println("")
    end

    F .+= H

    if debug && MPI.Comm_rank(comm) == 0
    #  println("Total Fock matrix:")
    #  for shell_group in 0:cld(size(F)[1],5)
    #    ending = min((5*shell_group + 5), size(F)[2])
    #    F_debug = F[:,(5*shell_group + 1):ending]
    #    pretty_table(hcat(collect(1:1:size(F)[1]),F_debug),
    #      vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
    #      formatter = ft_printf("%5.6f", collect(2:1:6)),
    #      highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    #  end
      h5write("debug.h5","SCF/$iter/F/Total", F)
    #  println("")
    end

    #== do DIIS ==#
    if ndiis > 0
      e .= F*D*S .- S*D*F

      @views e_array_old .= e_array[1:ndiis]
      test_e[1] = deepcopy(e)
      @views e_array .= vcat(test_e, e_array_old[1:ndiis-1])

      @views F_array_old .= F_array[1:ndiis]
      test_F[1] = deepcopy(F)
      @views F_array .= vcat(test_F, F_array_old[1:ndiis-1])

      if iter > 1
        B_dim += 1
        B_dim = min(B_dim,ndiis)
        try
          DIIS(F, e_array, F_array, B_dim)
        catch
          B_dim = 2
          DIIS(F, e_array, F_array, B_dim)
        end
      end
    end

    #== obtain new F,D,C matrices ==#
    #indices_tocopy = CartesianIndices((1:size(D,1),
    #  1:size(D,2)))
    #copyto!(D_old, indices_tocopy, D, indices_tocopy)
    D_old .= D

    E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo,
      ortho, basis, iter, debug)

    F .= F_mo
    #indices_tocopy = CartesianIndices((1:size(F_mo,1),
    #  1:size(F_mo,2)))
    #copyto!(F, indices_tocopy, F_mo, indices_tocopy)
  
    #== dynamic damping of density matrix ==#
    #D_damp[:] = map(x -> x*D[:,:] + (oneunit(typeof(dele))-x)*D_old[:,:],
    #  damp_values)
    #D_damp_rms = map(x->√(@∑ x-D_old x-D_old), D_damp)

    #x = maximum(D_damp_rms) > oneunit(typeof(dele)) ? minimum(damp_values) :
    #  maximum(damp_values)
    #D[:,:] = x*D[:,:] + (oneunit(typeof(dele))-x)*D_old[:,:]

    #== check for convergence ==#
    ΔD .= D .- D_old
    D_rms = √(@∑ ΔD ΔD)

    E = E_elec+E_nuc
    ΔE = E - E_old

    if MPI.Comm_rank(comm) == 0
      println(iter,"     ", E,"     ", ΔE,"     ", D_rms)
    end

    iter_converged = abs(ΔE) <= dele && D_rms <= rmsd
    iter += 1
    if iter > niter
      scf_converged = false
      break
    end

    #== if not converged, replace old D and E values for next iteration ==#
    E_old = E
  end

  return E
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

function twoei(F::Matrix{Float64}, D::Matrix{Float64}, H::Matrix{Float64},
  basis::BasisStructs.Basis; debug)

  fill!(F,zero(Float64))

  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3 #bitwise divide by 8

  ish_old = 0
  jsh_old = 0
  ksh_old = 0
  lsh_old = 0

  quartet_batch_num_old = ceil(Int64,nindices/QUARTET_BATCH_SIZE)

  #mutex = Base.Threads.ReentrantLock()
  thread_index_counter = nindices

  #for thread in 1:Threads.nthreads()
  #  F_priv = zeros(basis.norb,basis.norb)

    max_shell_am = MAX_SHELL_AM
    eri_quartet_batch = Vector{Float64}(undef,1296)

    bra = ShPair(basis.shells[1], basis.shells[1])
    ket = ShPair(basis.shells[1], basis.shells[1])
    quartet = ShQuartet(bra,ket)

    twoei_thread_kernel(F, D,
      H, basis, thread_index_counter, eri_quartet_batch,
      bra, ket, quartet, nindices,
      quartet_batch_num_old, ish_old, jsh_old, ksh_old, lsh_old; debug=debug)

    #lock(mutex)
    #F .+= F_priv
    #unlock(mutex)
  #end

  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      F[iorb,jorb] /= 2.0
      F[jorb,iorb] = F[iorb,jorb]
    end
  end

  return F
end

function twoei_thread_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  H::Matrix{Float64}, basis::BasisStructs.Basis, 
  thread_index_counter::Int64, 
  eri_quartet_batch::Vector{Float64}, bra::ShPair , ket::ShPair,
  quartet::ShQuartet, nindices::Int64, quartet_batch_num_old::Int64,
  ish_old::Int64, jsh_old::Int64, ksh_old::Int64, lsh_old::Int64; debug)

  comm=MPI.COMM_WORLD

  if debug println("START TWO-ELECTRON INTEGRALS") end
  while true
    #ijkl_index = Threads.atomic_sub!(thread_index_counter, 1)
    ijkl_index = thread_index_counter
    thread_index_counter -= 1
    if ijkl_index < 1 break end

    if MPI.Comm_rank(comm) != ijkl_index%MPI.Comm_size(comm) continue end
    bra_pair = decompose(ijkl_index)
    ket_pair = ijkl_index - triangular_index(bra_pair)

    ish = decompose(bra_pair)
    jsh = bra_pair - triangular_index(ish)

    ksh = decompose(ket_pair)
    lsh = ket_pair - triangular_index(ksh)

    bra.sh_a = basis[ish]
    bra.sh_b = basis[jsh]

    ket.sh_a = basis[ksh]
    ket.sh_b = basis[lsh]

    quartet.bra = bra
    quartet.ket = ket

    if debug
      if do_continue_print println("QUARTET: $ish, $jsh, $ksh, $lsh ($quartet_num):") end
    end

   # quartet_batch_num::Int64 = fld(quartet_num,
   #   QUARTET_BATCH_SIZE) + 1

    #if quartet_batch_num != quartet_batch_num_old
    #  if length(eri_starts) != QUARTET_BATCH_SIZE && length(eri_sizes) != QUARTET_BATCH_SIZE
    #    resize!(eri_sizes,QUARTET_BATCH_SIZE)
    #    resize!(eri_starts,QUARTET_BATCH_SIZE)
    #  end

    #  eri_sizes[:] = load("tei_batch.jld",
    #    "Sizes/$quartet_batch_num")

      #@views eri_starts[:] = [1, [ sum(eri_sizes[1:i])+1 for i in 1:(QUARTET_BATCH_SIZE-1)]... ]
      #eri_starts[:] = load("tei_batch.jld","Starts/$quartet_batch_num")
      #@views eri_starts[:] = eri_starts[:] .- (eri_starts[1] - 1)

    #  resize!(eri_batch,sum(eri_sizes))
    #  eri_batch[:] = load("tei_batch.jld","Integrals/$quartet_batch_num")

    #  quartet_batch_num_old = quartet_batch_num
    #end

    #quartet_num_in_batch::Int64 = quartet_num - QUARTET_BATCH_SIZE*
    #  (quartet_batch_num-1) + 1

    #starting::Int64 = eri_starts[quartet_num_in_batch]
    #ending::Int64 = starting + eri_sizes[quartet_num_in_batch] - 1
    #batch_ending_final::Int64 = ending - starting + 1

    #@views eri_quartet_batch[1:batch_ending_final] = eri_batch[starting:ending]
    #eri_quartet_batch = @view eri_batch[starting:ending]

    shellquart(ish, jsh, ksh, lsh, eri_quartet_batch)

    #eqb_size = bra.sh_a.am*bra.sh_b.am*ket.sh_a.am*ket.sh_b.am
    #@views eri_quartet_batch_abs[1:eqb_size] = abs.(eri_quartet_batch[1:eqb_size])

    dirfck(F, D, eri_quartet_batch, quartet,
      ish, jsh, ksh, lsh, debug)
  end
  if debug println("END TWO-ELECTRON INTEGRALS") end
end

@inline function shellquart(ish::Int64, jsh::Int64, ksh::Int64,
  lsh::Int64, eri_quartet_batch::Vector{Float64})

  #= actually compute integrals =#
  SIMINT.compute_eris(ish, jsh, ksh, lsh, eri_quartet_batch)
end


@inline function dirfck(F_priv::Matrix{Float64}, D::Matrix{Float64},
  eri_batch::Vector{Float64}, quartet::ShQuartet, ish::Int64, jsh::Int64,
  ksh::Int64, lsh::Int64, debug::Bool)

  norb = size(D,1)

  two_same = ish == jsh
  two_same = two_same || (ish == ksh)
  two_same = two_same || (ish == lsh)
  two_same = two_same || (jsh == ksh)
  two_same = two_same || (jsh == lsh)
  two_same = two_same || (ksh == lsh)

  three_same = ish == jsh && jsh == ksh
  three_same = three_same || (ish == jsh && jsh == lsh)
  three_same = three_same || (ish == ksh && ksh == lsh)
  three_same = three_same || (jsh == ksh && ksh == lsh)

  pμ = quartet.bra.sh_a.pos
  nμ = quartet.bra.sh_a.nbas

  pν = quartet.bra.sh_b.pos
  nν = quartet.bra.sh_b.nbas
  
  pλ = quartet.ket.sh_a.pos
  nλ = quartet.ket.sh_a.nbas
  
  pσ = quartet.ket.sh_b.pos
  nσ = quartet.ket.sh_b.nbas

  for μsize in 0:(nμ-1), νsize in 0:(nν-1)
    μμ = μsize + pμ
    νν = νsize + pν

    do_continue_bra = sort_bra(μμ, νν,
      ish, jsh, ksh, lsh, nμ, nν, nλ, nσ, two_same, three_same)
    if do_continue_bra continue end

    for λsize in 0:(nλ-1), σsize in 0:(nσ-1)
      λλ = λsize + pλ
      σσ = σsize + pσ

      #if debug
      #  if do_continue_print print("$μμ, $νν, $λλ, $σσ => ") end
      #end

      μνλσ = 1 + σsize + nσ*λsize + nσ*nλ*νsize +
        nσ*nλ*nν*μsize

      do_continue_screen = abs(eri_batch[μνλσ]) < 1E-10
      if do_continue_screen continue end

      μ, ν = (μμ > νν) ? (μμ, νν) : (νν, μμ)
      λ, σ = (λλ > σσ) ? (λλ, σσ) : (σσ, λλ)

      do_continue_ket = sort_ket(μμ, νν, λλ, σσ,
        ish, jsh, ksh, lsh, nμ, nν, nλ, nσ, two_same, three_same)
      if do_continue_ket continue end

      do_continue_braket, μ, ν, λ, σ = sort_braket(μμ, μ, νν, ν, λλ, λ, σσ, σ,
        ish, jsh, ksh, lsh, nμ, nν, nλ, nσ)
      if do_continue_braket continue end

      eri = eri_batch[μνλσ] 

      #if debug println("$μ, $ν, $λ, $σ, $eri") end
      eri *= (μ == ν) ? 0.5 : 1.0
      eri *= (λ == σ) ? 0.5 : 1.0
      eri *= ((μ == λ) && (ν == σ)) ? 0.5 : 1.0

      #λσ = λ + norb*(σ-1)
      #μν = μ + norb*(ν-1)
      #μλ = μ + norb*(λ-1)
      #μσ = μ + norb*(σ-1)
      #νλ = max(ν,λ) + norb*(min(ν,λ)-1)
      #νσ = max(ν,σ) + norb*(min(ν,σ)-1)

      F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
      F_priv[μ,ν] += 4.0 * D[λ,σ] * eri
      F_priv[μ,λ] -= D[ν,σ] * eri
      F_priv[μ,σ] -= D[ν,λ] * eri
      F_priv[max(ν,λ), min(ν,λ)] -= D[μ,σ] * eri
      F_priv[max(ν,σ), min(ν,σ)] -= D[μ,λ] * eri
    end
  end
end

#=
"""
	 iteration(F::Matrix{Float64}, D::Matrix{Float64}, H::Matrix{Float64}, ortho::Matrix{Float64})
Summary
======
Perform single SCF iteration.

Arguments
======
D = Current iteration's Density Matrix

H = One-electron Hamiltonian Matrix

ortho = Symmetric Orthogonalization Matrix
"""
=#
function iteration(F_μν::Matrix{Float64}, D::Matrix{Float64},
  C::Matrix{Float64}, H::Matrix{Float64}, F_eval::Vector{Float64},
  F_evec::Matrix{Float64}, F_mo::Matrix{Float64}, ortho::Matrix{Float64},
  basis::BasisStructs.Basis, iter::Int, debug::Bool)

  comm=MPI.COMM_WORLD

  #== obtain new orbital coefficients ==#
  F_mo .= transpose(ortho)*F_μν*ortho

  F_eval .= eigvals(LinearAlgebra.Hermitian(F_mo))

  F_evec .= eigvecs(LinearAlgebra.Hermitian(F_mo))
  F_evec .= F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

  #copyto!(F_evec, CartesianIndices((1:size(F_evec,1), 1:size(F_evec,2))),
  #  F_evec[:,sortperm(F_eval)], CartesianIndices((1:size(F_evec,1), 1:size(F_evec,2))))

  C .= ortho*F_evec

  if debug && MPI.Comm_rank(comm) == 0
  #  println("New orbitals:")
  #  for shell_group in 0:cld(size(C)[1],5)
  #    ending = min((5*shell_group + 5), size(C)[2])
  #    C_debug = C[:,(5*shell_group + 1):ending]
  #    pretty_table(hcat(collect(1:1:size(C)[1]),C_debug),
  #      vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
  #      formatter = ft_printf("%5.6f", collect(2:1:6)),
  #      highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
  #  end
    h5write("debug.h5","SCF/$iter/C", C)
  #  println("")
  end

  #== build new density matrix ==#
  nocc = basis.nels >> 1
  norb = basis.norb

  for i in 1:basis.norb, j in 1:basis.norb
    @views D[i,j] = @∑ C[i,1:nocc] C[j,1:nocc]
    #D[i,j] = @∑ C[1:nocc,i] C[1:nocc,j]
    D[i,j] *= 2.0
  end

  #if debug && MPI.Comm_rank(comm) == 0
  #  println("New density matrix:")
  #  for shell_group in 0:cld(size(D)[1],5)
  #    ending = min((5*shell_group + 5), size(D)[2])
  #    D_debug = D[:,(5*shell_group + 1):ending]
  #    pretty_table(hcat(collect(1:1:size(D)[1]),D_debug),
  #      vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
  #      formatter = ft_printf("%5.6f", collect(2:1:6)),
  #      highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
  #  end
  #  h5write("debug.h5","SCF/$iter/D", D)
  #  println("")
  #end

  #== compute new SCF energy ==#
  EHF1 = @∑ D F_μν
  EHF2 = @∑ D H
  E_elec = (EHF1 + EHF2)/2.0

  #if debug && MPI.Comm_rank(comm) == 0
    #println("New energy:")
    #println("$EHF1, $EHF2")
    #h5write("debug.h5","SCF/$iter/E/EHF1", EHF1)
    #h5write("debug.h5","SCF/$iter/E/EHF2", EHF2)
    #h5write("debug.h5","SCF/$iter/E/EHF", E_elec)
    #println("")
  #end

  return E_elec
end
