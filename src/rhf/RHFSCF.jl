using MATH
using JCModules.Globals

using MPI
using Base.Threads
#using Distributed
using LinearAlgebra
using JLD
using PrettyTables

function rhf_energy(basis::BasisStructs.Basis,
  molecule::Union{Dict{String,Any},Dict{Any,Any}}, scf_flags::Dict{String,Any})

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
  molecule::Union{Dict{String,Any},Dict{Any,Any}}, scf_flags::Dict{String,Any})

  comm=MPI.COMM_WORLD
  calculation_status = Dict([])

  #== read variables from input if needed ==#
  E_nuc = molecule["enuc"]

  S = read_in_oei(molecule["ovr"], basis.norb)
  H = read_in_oei(molecule["hcore"], basis.norb)

  if scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0
    println("Overlap matrix:")
    for shell_group in 0:cld(size(S)[1],5)
      ending = min((5*shell_group + 5), size(S)[2])
      S_debug = S[:,(5*shell_group + 1):ending]
      pretty_table(hcat(collect(1:1:size(S)[1]),S_debug),
        vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
        formatter = ft_printf("%5.6f", collect(2:1:6)),
        highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    end
    #display(S)
    println("")

    println("Hamiltonian matrix:")
    #display(H)
    for shell_group in 0:cld(size(H)[1],5)
      ending = min((5*shell_group + 5), size(H)[2])
      H_debug = H[:,(5*shell_group + 1):ending]
      pretty_table(hcat(collect(1:1:size(H)[1]),H_debug),
        vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
        formatter = ft_printf("%5.6f", collect(2:1:6)),
        highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    end
    println("")
  end

  #== build the orthogonalization matrix ==#
  S_evec = eigvecs(LinearAlgebra.Hermitian(S))

  S_eval_diag = eigvals(LinearAlgebra.Hermitian(S))

  S_eval = zeros(basis.norb,basis.norb)
  for i in 1:basis.norb
    S_eval[i,i] = S_eval_diag[i]
  end

  ortho = Matrix{Float64}(undef, basis.norb, basis.norb)
  @views ortho[:,:] = S_evec[:,:]*
    (LinearAlgebra.Diagonal(S_eval)^-0.5)[:,:]*transpose(S_evec)[:,:]

  if scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0
    println("Ortho matrix:")
    #display(ortho)
    for shell_group in 0:cld(size(ortho)[1],5)
      ending = min((5*shell_group + 5), size(ortho)[2])
      ortho_debug = ortho[:,(5*shell_group + 1):ending]
      pretty_table(hcat(collect(1:1:size(ortho)[1]),ortho_debug),
        vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))),
        formatter = ft_printf("%5.6f", collect(2:1:6)),
        highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    end
    println("")
  end

  #== build the initial matrices ==#
  F = H
  F_eval = Vector{Float64}(undef,basis.norb)
  F_evec = Matrix{Float64}(undef,basis.norb,basis.norb)
  F_mo = Matrix{Float64}(undef,basis.norb,basis.norb)

  D = Matrix{Float64}(undef,basis.norb,basis.norb)
  C = Matrix{Float64}(undef,basis.norb,basis.norb)

  if MPI.Comm_rank(comm) == 0
    println("----------------------------------------          ")
    println("       Starting RHF iterations...                 ")
    println("----------------------------------------          ")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")
  end

  E_elec = 0.0
  E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo, ortho, basis,
    scf_flags)
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
  @time F, D, C, E, converged = scf_cycles(F, D, C, E, H, ortho, S, F_eval,
  F_evec, F_mo, E_nuc, E_elec, E_old, basis, scf_flags)

  if !converged
    iter_limit = scf_flags["niter"]

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
      "error_message" => " SCF calculation did not converge within $iter_limit
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

    #if (FLAGS.SCF.debug == true)
    #  close(json_debug)
    #end
  end

  return F, D, C, E, calculation_status
end

function scf_cycles(F::Matrix{Float64}, D::Matrix{Float64}, C::Matrix{Float64},
  E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64}, S::Matrix{Float64},
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, F_mo::Matrix{Float64},
  E_nuc::Float64, E_elec::Float64, E_old::Float64, basis::BasisStructs.Basis,
  scf_flags::Dict{String,Any})

  #== build DIIS arrays ==#
  ndiis = scf_flags["ndiis"]
  F_array = fill(Matrix{Float64}(undef,basis.norb,basis.norb),
    ndiis)

  e = Matrix{Float64}(undef,basis.norb,basis.norb)
  e_array = fill(
    Matrix{Float64}(undef,basis.norb,basis.norb), ndiis)
  e_array_old = fill(
    Matrix{Float64}(undef,basis.norb,basis.norb), ndiis)
  F_array_old = fill(
    Matrix{Float64}(undef,basis.norb,basis.norb), ndiis)

  #== build arrays needed for post-fock build iteration calculations ==#
  F_temp = Matrix{Float64}(undef,basis.norb,basis.norb)
  D_old = Matrix{Float64}(undef,basis.norb,basis.norb)
  ΔD = Matrix{Float64}(undef,basis.norb,basis.norb)

  #== build arrays needed for dynamic damping ==#
  damp_values = [ 0.25, 0.75 ]
  D_damp = [ Matrix{Float64}(undef,basis.norb,basis.norb)
    for i in 1:2 ]
  D_damp_rms = [ zero(Float64), zero(Float64) ]

  #== build variables needed for eri batching ==#
  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3

  quartet_batch_num_old = fld(nindices,
    QUARTET_BATCH_SIZE) + 1

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
    E_elec, E_old, basis, scf_flags, ndiis, F_array, e, e_array, e_array_old,
    F_array_old, F_temp, F_eval, F_evec, F_mo, D_old, ΔD, damp_values, D_damp,
    D_damp_rms, scf_converged, quartet_batch_num_old)

  #== we are done! ==#
  return F, D, C, E, scf_converged
end

function scf_cycles_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  C::Matrix{Float64}, E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64},
  S::Matrix{Float64}, E_nuc::Float64, E_elec::Float64, E_old::Float64,
  basis::BasisStructs.Basis, scf_flags::Dict{String,Any}, ndiis::Int64,
  F_array::Vector{Matrix{Float64}}, e::Matrix{Float64},
  e_array::Vector{Matrix{Float64}}, e_array_old::Vector{Matrix{Float64}},
  F_array_old::Vector{Matrix{Float64}}, F_temp::Matrix{Float64},
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, F_mo::Matrix{Float64},
  D_old::Matrix{Float64}, ΔD::Matrix{Float64}, damp_values::Vector{Float64},
  D_damp::Vector{Matrix{Float64}}, D_damp_rms::Vector{Float64},
  scf_converged::Bool, quartet_batch_num_old::Int64)

  #== initialize a few more variables ==#
  comm=MPI.COMM_WORLD

  iter_limit = scf_flags["niter"]
  dele = scf_flags["dele"]
  rmsd = scf_flags["rmsd"]

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
    F_temp[:,:] = twoei(F, D, H, basis)

    F[:,:] = MPI.Allreduce(F_temp[:,:],MPI.SUM,comm)
    MPI.Barrier(comm)

    if scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0
      println("Skeleton Fock matrix:")
      #display(F)
      for shell_group in 0:cld(size(F)[1],5)
        ending = min((5*shell_group + 5), size(F)[2]) 
        F_debug = F[:,(5*shell_group + 1):ending]
        pretty_table(hcat(collect(1:1:size(F)[1]),F_debug), 
          vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))), 
          formatter = ft_printf("%5.6f", collect(2:1:6)), 
          highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
      end
      println("")
    end

    F[:,:] .+= H[:,:]

    if scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0
      println("Total Fock matrix:")
      #display(F)
      for shell_group in 0:cld(size(F)[1],5)
        ending = min((5*shell_group + 5), size(F)[2]) 
        F_debug = F[:,(5*shell_group + 1):ending]
        pretty_table(hcat(collect(1:1:size(F)[1]),F_debug), 
          vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))), 
          formatter = ft_printf("%5.6f", collect(2:1:6)), 
          highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
      end
      println("")
    end

    #== do DIIS ==#
    if ndiis > 0
      e[:,:] = F[:,:]*D[:,:]*S[:,:] .- S[:,:]*D[:,:]*F[:,:]

      e_array_old[:] = e_array[1:ndiis]
      e_array[:] = [deepcopy(e), e_array_old[1:ndiis-1]...]

      F_array_old[:] = F_array[1:ndiis]
      F_array[:] = [deepcopy(F), F_array[1:ndiis-1]...]

      if iter > 1
        B_dim += 1
        B_dim = min(B_dim,ndiis)
        try
          F[:,:] = DIIS(e_array, F_array, B_dim)
        catch
          B_dim = 2
          F[:,:] = DIIS(e_array, F_array, B_dim)
        end
      end
    end

    #== obtain new F,D,C matrices ==#
    indices_tocopy = CartesianIndices((1:size(D,1),
      1:size(D,2)))
    copyto!(D_old, indices_tocopy, D, indices_tocopy)

    E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo,
      ortho, basis, scf_flags)

    #F = deepcopy(F_mo)
    indices_tocopy = CartesianIndices((1:size(F_mo,1),
      1:size(F_mo,2)))
    copyto!(F, indices_tocopy, F_mo, indices_tocopy)

    #== dynamic damping of density matrix ==#
    #D_damp[:] = map(x -> x*D[:,:] + (oneunit(typeof(dele))-x)*D_old[:,:],
    #  damp_values)
    #D_damp_rms = map(x->√(@∑ x-D_old x-D_old), D_damp)

    #x::T = maximum(D_damp_rms) > oneunit(typeof(dele)) ? minimum(damp_values) :
    #  maximum(damp_values)
    #D[:,:] = x*D[:,:] + (oneunit(typeof(dele))-x)*D_old[:,:]

    #== check for convergence ==#
    @views ΔD[:,:] = D[:,:] .- D_old[:,:]
    D_rms = √(@∑ ΔD ΔD)

    E = E_elec+E_nuc
    ΔE = E - E_old

    if MPI.Comm_rank(comm) == 0
      println(iter,"     ", E,"     ", ΔE,"     ", D_rms)
    end

    iter_converged = (abs(ΔE) <= dele) && (D_rms <= rmsd)
    iter += 1
    if iter > iter_limit
      scf_converged = false
      break
    end

    #== if not converged, replace old D and E values for next iteration ==#
    #D_old = deepcopy(D)
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
  basis::BasisStructs.Basis)

  fill!(F,zero(Float64))

  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3 #bitwise divide by 8

  quartet_batch_num_old = fld(nindices,
    QUARTET_BATCH_SIZE) + 1

  mutex = Base.Threads.Mutex()
  thread_index_counter = Threads.Atomic{Int64}(nindices)

  Threads.@threads for thread in 1:Threads.nthreads()
    F_priv = zeros(basis.norb,basis.norb)

    max_shell_am = MAX_SHELL_AM
    eri_quartet_batch = Vector{Float64}(undef,1296)

    bra = ShPair(basis.shells[1], basis.shells[1])
    ket = ShPair(basis.shells[1], basis.shells[1])
    quartet = ShQuartet(bra,ket)

    twoei_thread_kernel(F, D,
      H, basis, mutex, thread_index_counter, F_priv, eri_quartet_batch,
      bra, ket, quartet, nindices, quartet_batch_num_old)

    lock(mutex)
    F[:,:] .+= F_priv[:,:]
    unlock(mutex)
  end

  for iorb in 1:basis.norb, jorb in 1:basis.norb
    if iorb != jorb F[iorb,jorb] /= 2.0 end
  end

  return F
end

@inline function twoei_thread_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  H::Matrix{Float64}, basis::BasisStructs.Basis, mutex::Base.Threads.Mutex,
  thread_index_counter::Threads.Atomic{Int64}, F_priv::Matrix{Float64},
  eri_quartet_batch::Vector{Float64}, bra::ShPair , ket::ShPair,
  quartet::ShQuartet, nindices::Int64, quartet_batch_num_old::Int64)

  comm=MPI.COMM_WORLD

  println("START TWO-ELECTRON INTEGRALS")
  while true
    ijkl_index = Threads.atomic_sub!(thread_index_counter, 1)
    if ijkl_index < 1 break end

    if MPI.Comm_rank(comm) != ijkl_index%MPI.Comm_size(comm) continue end
    bra_pair = decompose(ijkl_index)
    ket_pair = ijkl_index - triangular_index(bra_pair)

    ish = decompose(bra_pair)
    jsh = bra_pair - triangular_index(ish)

    ksh = decompose(ket_pair)
    lsh = ket_pair - triangular_index(ksh)

    ijsh = triangular_index(ish,jsh)
    klsh = triangular_index(ksh,lsh)

    if klsh > ijsh ish,jsh,ksh,lsh = ksh,lsh,ish,jsh end

    bra.sh_a = basis[ish]
    bra.sh_b = basis[jsh]

    ket.sh_a = basis[ksh]
    ket.sh_b = basis[lsh]

    quartet.bra = bra
    quartet.ket = ket

    qnum_ij = triangular_index(ish, jsh)
    qnum_kl = triangular_index(ksh, lsh)
    quartet_num = triangular_index(qnum_ij, (qnum_kl - 1))

    #println("QUARTET: $ish, $jsh, $ksh, $lsh ($quartet_num):")

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

    shellquart_direct(ish,jsh,ksh,lsh,eri_quartet_batch)

    #if abs(maximum(eri_quartet_batch)) > 1E-10
      dirfck(F_priv, D, eri_quartet_batch, quartet,
        ish, jsh, ksh, lsh)
    #end
  end
  println("END TWO-ELECTRON INTEGRALS")
end

@inline function shellquart_direct(ish::Int64, jsh::Int64, ksh::Int64,
  lsh::Int64, eri_quartet_batch::Vector{Float64})

  SIMINT.retrieve_eris(ish, jsh, ksh, lsh, eri_quartet_batch)
end


@noinline function dirfck(F_priv::Matrix{Float64}, D::Matrix{Float64},
  eri_batch::Vector{Float64}, quartet::ShQuartet, ish::Int64, jsh::Int64,
  ksh::Int64, lsh::Int64)

  norb = size(D)[1]

  spμ = quartet.bra.sh_a.sp
  spν = quartet.bra.sh_b.sp
  spλ = quartet.ket.sh_a.sp
  spσ = quartet.ket.sh_b.sp

  μνλσ = 0

  for spi in 0:spμ, spj in 0:spν
    nμ = 0
    pμ = quartet.bra.sh_a.pos
    if spμ == 1
      nμ = spi == 1 ? 3 : 1
      pμ += spi == 1 ? 1 : 0
    else
      nμ = quartet.bra.sh_a.nbas
    end

    nν = 0
    pν = quartet.bra.sh_b.pos
    if spν == 1
      nν = spj == 1 ? 3 : 1
      pν += spj == 1 ? 1 : 0
    else
      nν = quartet.bra.sh_b.nbas
    end

    for spk in 0:spλ, spl in 0:spσ
      nλ = 0
      pλ = quartet.ket.sh_a.pos
      if spλ == 1
        nλ = spk == 1 ? 3 : 1
        pλ += spk == 1 ? 1 : 0
      else
        nλ = quartet.ket.sh_a.nbas
      end

      nσ = 0
      pσ = quartet.ket.sh_b.pos
      if spσ == 1
        nσ = spl == 1 ? 3 : 1
        pσ += spl == 1 ? 1 : 0
      else
        nσ = quartet.ket.sh_b.nbas
      end

      for μμ in pμ:pμ+(nμ-1), νν in pν:pν+(nν-1)
        μ, ν = μμ,νν
        #if (μμ < νν) μ, ν = ν, μ end

        μν = triangular_index(μμ,νν)

        for λλ in pλ:pλ+(nλ-1), σσ in pσ:pσ+(nσ-1)
          λ, σ = λλ,σσ
          #if (λλ < σσ) λ, σ = σ, λ end

          λσ = triangular_index(λλ,σσ)

          #if (μν < λσ) μ, ν, λ, σ = λ, σ, μ, ν end
          #if (μν < λσ)
          #  μνλσ += 1
          #  continue
          #end
          #print("$μμ, $νν, $λλ, $σσ => ")
          if (μμ < νν)
            μνλσ += 1
            #println("DO CONTINUE")
            continue
          end

          if (λλ < σσ)
            μνλσ += 1
            #println("DO CONTINUE")
            continue
          end

          if (μν < λσ)
            do_continue = false

            do_continue, μ, ν, λ, σ = sort_braket(μμ, νν, λλ, σσ, ish, jsh,
              ksh, lsh, nμ, nν, nλ, nσ)

            if do_continue
              μνλσ += 1
              #println("DO CONTINUE")
              continue
            end
          end

          μνλσ += 1

	        eri = eri_batch[μνλσ]
          #eri::T = 0
          if (abs(eri) <= 1E-10) continue end

          
          printit = μ==23 && ν==8
          printit = printit || (μ==23 && λ==8)
          printit = printit || (μ==23 && σ==8)
          printit = printit || (ν==23 && λ==8)
          printit = printit || (ν==23 && σ==8)
          printit = printit || (λ==23 && σ==8)
          if true 
            println("$μ, $ν, $λ, $σ, $eri")
          end
	        eri *= (μ == ν) ? 0.5 : 1.0
	        eri *= (λ == σ) ? 0.5 : 1.0
	        eri *= ((μ == λ) && (ν == σ)) ? 0.5 : 1.0

	        F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
	        F_priv[μ,ν] += 4.0 * D[λ,σ] * eri
          F_priv[μ,λ] -= D[ν,σ] * eri
	        F_priv[μ,σ] -= D[ν,λ] * eri
          F_priv[ν,λ] -= D[max(μ,σ),min(μ,σ)] * eri
          F_priv[ν,σ] -= D[max(μ,λ),min(μ,λ)] * eri

          if λ != σ F_priv[σ,λ] += 4.0 * D[μ,ν] * eri end
	        if μ != ν F_priv[ν,μ] += 4.0 * D[λ,σ] * eri end
          if μ != λ F_priv[λ,μ] -= D[ν,σ] * eri end
	        if μ != σ F_priv[σ,μ] -= D[ν,λ] * eri end
          if ν != λ F_priv[λ,ν] -= D[max(μ,σ),min(μ,σ)] * eri end
	        if ν != σ F_priv[σ,ν] -= D[max(μ,λ),min(μ,λ)] * eri end
        end
      end
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
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

H = One-electron Hamiltonian Matrix

ortho = Symmetric Orthogonalization Matrix
"""
=#
function iteration(F_μν::Matrix{Float64}, D::Matrix{Float64},
  C::Matrix{Float64}, H::Matrix{Float64}, F_eval::Vector{Float64},
  F_evec::Matrix{Float64}, F_mo::Matrix{Float64}, ortho::Matrix{Float64},
  basis::BasisStructs.Basis, scf_flags::Dict{String,Any})

  comm=MPI.COMM_WORLD

  #== obtain new orbital coefficients ==#
  @views F_mo[:,:] = transpose(ortho)[:,:]*F_μν[:,:]*ortho[:,:]

  F_eval[:] = eigvals(LinearAlgebra.Hermitian(F_mo))

  @views F_evec[:,:] = eigvecs(LinearAlgebra.Hermitian(F_mo))[:,:]
  @views F_evec[:,:] = F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

  #copyto!(F_evec, CartesianIndices((1:size(F_evec,1), 1:size(F_evec,2))),
  #  F_evec[:,sortperm(F_eval)], CartesianIndices((1:size(F_evec,1), 1:size(F_evec,2))))

  @views C[:,:] = ortho[:,:]*F_evec[:,:]

  if (scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0)
    println("New orbitals:")
    #display(C)
    for shell_group in 0:cld(size(C)[1],5)
      ending = min((5*shell_group + 5), size(C)[2]) 
      C_debug = C[:,(5*shell_group + 1):ending]
      pretty_table(hcat(collect(1:1:size(C)[1]),C_debug), 
        vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))), 
        formatter = ft_printf("%5.6f", collect(2:1:6)), 
        highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    end
    println("")
  end

  #== build new density matrix ==#
  nocc = div(basis.nels,2)
  norb = basis.norb

  for i in 1:basis.norb, j in 1:basis.norb
    @views D[i,j] = @∑ C[i,1:nocc] C[j,1:nocc]
    #D[i,j] = @∑ C[1:nocc,i] C[1:nocc,j]
    D[i,j] *= 2
  end

  if (scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0)
    println("New density matrix:")
    #display(D)
    for shell_group in 0:cld(size(D)[1],5)
      ending = min((5*shell_group + 5), size(D)[2]) 
      D_debug = D[:,(5*shell_group + 1):ending]
      pretty_table(hcat(collect(1:1:size(D)[1]),D_debug), 
        vcat( [ "Shell" ], map( x -> "$x", collect((5*shell_group+1):1:ending))), 
        formatter = ft_printf("%5.6f", collect(2:1:6)), 
        highlighters = Highlighter((data,i,j)-> (i < j), crayon"black"))
    end
    println("")
  end

  #== compute new SCF energy ==#
  EHF1 = @∑ D F_μν
  EHF2 = @∑ D H
  E_elec = (EHF1 + EHF2)/2

  if (scf_flags["debug"] == true && MPI.Comm_rank(comm) == 0)
    println("New energy:")
    println("$EHF1, $EHF2")
    println("")
  end

  return E_elec
end
