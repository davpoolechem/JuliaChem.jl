using MATH
using JCModules.Globals

using Base.Threads
using LinearAlgebra
using HDF5
using PrettyTables

const do_continue_print = false 

function rhf_energy(mol::MolStructs.Molecule, basis::BasisStructs.Basis,
  scf_flags::Union{Dict{String,Any},Dict{Any,Any}}; output)

  return rhf_kernel(mol,basis,scf_flags; output=output)
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
@inline function rhf_kernel(mol::MolStructs.Molecule, basis::BasisStructs.Basis,
  scf_flags::Union{Dict{String,Any},Dict{Any,Any}}; output)

  comm=MPI.COMM_WORLD
  calculation_status = Dict([])

  #== read in some variables from scf input ==#
  debug::Bool = scf_flags["debug"]
  debug_output = debug ? h5open("debug.h5","w") : nothing

  niter::Int = scf_flags["niter"]

  #== compute nuclear repulsion energy ==# 
  #E_nuc::Float64 = molecule["enuc"]
  E_nuc = compute_enuc(mol)
  
  #S = read_in_oei(molecule["ovr"], basis.norb)
  
  #== compute one-electron integrals and Hamiltonian ==#
  S = zeros(Float64, (basis.norb, basis.norb))
  compute_overlap(S, basis)
  
  T = zeros(Float64, (basis.norb, basis.norb))
  compute_ke(T, basis)
 
  V = zeros(Float64, (basis.norb, basis.norb))
  compute_nah(V, mol, basis)

  #H = read_in_oei(molecule["hcore"], basis.norb)
  H = T .+ V
 
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-None/E_nuc", E_nuc)
    h5write("debug.h5","SCF/Iteration-None/S", S)
    h5write("debug.h5","SCF/Iteration-None/T", T)
    h5write("debug.h5","SCF/Iteration-None/V", V)
    h5write("debug.h5","SCF/Iteration-None/H", H)
  end

  #== build the orthogonalization matrix ==#
  S_evec = eigvecs(LinearAlgebra.Hermitian(S))

  S_eval_diag = eigvals(LinearAlgebra.Hermitian(S))

  S_eval = zeros(basis.norb,basis.norb)
  for i in 1:basis.norb
    S_eval[i,i] = S_eval_diag[i]
  end
  
  ortho = S_evec*(LinearAlgebra.Diagonal(S_eval)^-0.5)*transpose(S_evec)
  
  ortho_trans = transpose(LinearAlgebra.Hermitian(ortho))

  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-None/Ortho", ortho)
  end

  #== build the initial matrices ==#
  F = H
  F_eval = Vector{Float64}(undef,basis.norb)
  F_evec = similar(F)
  F_mo = similar(F)
  F_part = similar(F)

  D = similar(F)
  C = similar(F)

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("----------------------------------------          ")
    println("       Starting RHF iterations...                 ")
    println("----------------------------------------          ")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")
  end

  E_elec = 0.0
  #@code_warntype iteration(F, D, C, H, F_eval, F_evec, F_mo, F_part, ortho, 
  E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo, F_part, ortho, 
    ortho_trans.data, basis, 0, debug)
  
  F = deepcopy(F) #apparently this is needed for some reason
  F_old = deepcopy(F)

  E = E_elec + E_nuc
  E_old = E

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println(0,"     ", E)
  end

  #=============================#
  #== start scf cycles: #7-10 ==#
  #=============================#
  #@code_warntype scf_cycles(F, D, C, E, H, ortho, ortho_trans.data, S, 
  F, D, C, E, converged = scf_cycles(F, D, C, E, H, ortho, ortho_trans.data, S, 
    F_eval, F_evec, F_mo, F_part, F_old, E_nuc, E_elec, E_old, basis, 
    scf_flags; output=output, debug=debug, niter=niter)

  if !converged
    if MPI.Comm_rank(comm) == 0 && output != "none"
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
    if MPI.Comm_rank(comm) == 0 && output != "none" 
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

  scf = Dict("Fock" => F,                                                       
             "Density" => D,                                                    
             "MO Coeff" => C,                                                   
             "Energy" => E,                                                     
             "Converged?" => converged                                      
            )                                                                   
                                                                                
  return scf 
end

function scf_cycles(F::Matrix{Float64}, D::Matrix{Float64}, C::Matrix{Float64},
  E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64}, 
  ortho_trans::Matrix{Float64}, S::Matrix{Float64}, F_eval::Vector{Float64}, 
  F_evec::Matrix{Float64}, F_mo::Matrix{Float64}, F_part::Matrix{Float64},
  F_old::Matrix{Float64}, E_nuc::Float64, E_elec::Float64, E_old::Float64, 
  basis::BasisStructs.Basis, scf_flags::Union{Dict{String,Any},Dict{Any,Any}}; 
  output, debug, niter)

  #== read in some more variables from scf flags input ==#
  ndiis::Int = scf_flags["ndiis"]
  dele::Float64 = scf_flags["dele"]
  rmsd::Float64 = scf_flags["rmsd"]
  load::String = scf_flags["load"]

  #== build variables needed for eri batching ==#
  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3

  #== build DIIS arrays ==#
  F_array = fill(similar(F), ndiis)

  e = similar(F)
  test_e = [ similar(F) ]
  e_array = fill(similar(F), ndiis)
  e_array_old = fill(similar(F), ndiis)
  
  test_F = [ similar(F) ]
  F_array_old = fill(similar(F), ndiis)

  FD = similar(F)
  FDS = similar(F)
  SDF = similar(F)
  
  #== build arrays needed for post-fock build iteration calculations ==#
  F_temp = similar(F)
  D_old = similar(F)
  ΔD = similar(F)

  #== build matrix of Cauchy-Schwarz upper bounds ==# 
  schwarz_bounds = zeros(Float64,(nsh,nsh)) 
  compute_schwarz_bounds(schwarz_bounds, nsh)

  Dsh = similar(schwarz_bounds)

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

  #@code_warntype scf_cycles_kernel(F, D, C, E, H, ortho, ortho_trans, S, E_nuc,
  E = scf_cycles_kernel(F, D, C, E, H, ortho, ortho_trans, S, E_nuc,
    E_elec, E_old, basis, F_array, e, e_array, e_array_old,
    F_array_old, F_temp, F_eval, F_evec, F_mo, F_part, F_old, D_old, ΔD, 
    scf_converged, test_e, test_F, FD, FDS, SDF, schwarz_bounds, Dsh; 
    output=output, debug=debug, niter=niter, ndiis=ndiis, dele=dele, 
    rmsd=rmsd, load=load)

  #== we are done! ==#
  if debug
    h5write("debug.h5","SCF/Iteration-Final/F", F)
    h5write("debug.h5","SCF/Iteration-Final/D", D)
    h5write("debug.h5","SCF/Iteration-Final/C", C)
    h5write("debug.h5","SCF/Iteration-Final/E", E)
    h5write("debug.h5","SCF/Iteration-Final/converged", scf_converged)
  end

  return F, D, C, E, scf_converged
end

function scf_cycles_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  C::Matrix{Float64}, E::Float64, H::Matrix{Float64}, ortho::Matrix{Float64},
  ortho_trans::Matrix{Float64}, S::Matrix{Float64}, E_nuc::Float64, 
  E_elec::Float64, E_old::Float64, basis::BasisStructs.Basis,
  F_array::Vector{Matrix{Float64}}, e::Matrix{Float64},
  e_array::Vector{Matrix{Float64}}, e_array_old::Vector{Matrix{Float64}},
  F_array_old::Vector{Matrix{Float64}}, F_temp::Matrix{Float64},
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, F_mo::Matrix{Float64}, 
  F_part::Matrix{Float64}, F_old::Matrix{Float64},
  D_old::Matrix{Float64}, ΔD::Matrix{Float64}, scf_converged::Bool,  
  test_e::Vector{Matrix{Float64}}, test_F::Vector{Matrix{Float64}},
  FD::Matrix{Float64}, FDS::Matrix{Float64}, SDF::Matrix{Float64}, 
  schwarz_bounds::Matrix{Float64}, Dsh::Matrix{Float64};
  output, debug, niter, ndiis, dele, rmsd, load)

  #== initialize a few more variables ==#
  comm=MPI.COMM_WORLD

  B_dim = 1
  D_rms = 1.0
  ΔE = 1.0 

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

    #== compress D into shells in Dsh ==#
    for ish in 1:length(basis.shells), jsh in 1:ish
      ipos = basis[ish].pos
      ibas = basis[ish].nbas

      jpos = basis[jsh].pos
      jbas = basis[jsh].nbas

      @views Dsh[ish, jsh] = maximum(D[ipos:(ipos+ibas-1),jpos:(jpos+jbas-1)])
      Dsh[jsh, ish] = Dsh[ish, jsh] 
    end
  
    #== build new Fock matrix ==#
    F_temp .= fock_build(F, D, H, basis, schwarz_bounds, Dsh; 
      debug=debug, load=load)

    F .= MPI.Allreduce(F_temp,MPI.SUM,comm)
    MPI.Barrier(comm)

    if debug && MPI.Comm_rank(comm) == 0
      h5write("debug.h5","SCF/Iteration-$iter/F/Skeleton", F)
    end

    F .+= H

    if debug && MPI.Comm_rank(comm) == 0
      h5write("debug.h5","SCF/Iteration-$iter/F/Total", F)
    end

    #== do DIIS ==#
    if ndiis > 0
      BLAS.symm!('L', 'U', 1.0, F, D, 0.0, FD)
      BLAS.gemm!('N', 'N', 1.0, FD, S, 0.0, FDS)
      
      SDF .= transpose(FDS)
      
      e .= FDS .- SDF

      @views e_array_old .= e_array[1:ndiis]
      test_e[1] .= e
      @views e_array .= vcat(deepcopy(test_e), e_array_old[1:ndiis-1])

      @views F_array_old .= F_array[1:ndiis]
      test_F[1] .= F
      @views F_array .= vcat(deepcopy(test_F), F_array_old[1:ndiis-1])

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

    #== dynamic damping of Fock matrix ==#
    x = ΔE >= 1.0 ? 0.9/log(50,50*ΔE) : 0.9
    F .= x.*F .+ (1.0-x).*F_old 

    F_old .= F

    #== obtain new F,D,C matrices ==#
    D_old .= D

    E_elec = iteration(F, D, C, H, F_eval, F_evec, F_mo, F_part,
      ortho, ortho_trans, basis, iter, debug)

    #== check for convergence ==#
    ΔD .= D .- D_old
    D_rms = √(@∑ ΔD ΔD)

    E = E_elec+E_nuc
    ΔE = E - E_old

    if MPI.Comm_rank(comm) == 0 && output == "verbose"
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
	 fock_build(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
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

@inline function fock_build(F::Matrix{Float64}, D::Matrix{Float64}, 
  H::Matrix{Float64}, basis::BasisStructs.Basis, 
  schwarz_bounds::Matrix{Float64}, Dsh::Matrix{Float64}; debug, load)

  comm = MPI.COMM_WORLD
  
  fill!(F,zero(Float64))

  nsh = length(basis.shells)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3 #bitwise divide by 8
 
  mutex = Base.Threads.ReentrantLock()
  thread_index_counter = Threads.Atomic{Int64}(nindices)
  
  #== simply do calculation for single-rank runs ==#
  if load == "static" 
    ish_old = 0
    jsh_old = 0
    ksh_old = 0
    lsh_old = 0

    max_shell_am = MAX_SHELL_AM
    eri_quartet_batch = Vector{Float64}(undef,1296)

    quartet = ShQuartet(ShPair(basis.shells[1], basis.shells[1]),
      ShPair(basis.shells[1], basis.shells[1]))
 
    simint_workspace = Vector{Float64}(undef,100000)

    while true 
      ijkl_index = Threads.atomic_sub!(thread_index_counter, 1)
 
      if ijkl_index <= 0 break
      elseif MPI.Comm_rank(comm) != ijkl_index%MPI.Comm_size(comm) continue 
      end

      fock_build_thread_kernel(F, D,
        H, basis, eri_quartet_batch, mutex,
        quartet, ijkl_index, simint_workspace, schwarz_bounds, Dsh,
        ish_old, jsh_old, ksh_old, lsh_old; debug=debug)
    end
      
    #lock(mutex)
    #  F .+= F_priv
    #unlock(mutex)
  #== otherwise use dynamic task distribution with a master/slave model ==#
  elseif load == "dynamic"
    batch_size = ceil(Int,nindices/(MPI.Comm_size(comm)*10000)) 

    #== master rank ==#
    if MPI.Comm_rank(comm) == 0 
      #== send out initial tasks to slaves ==#
      task = [ nindices ]
      initial_task = 1
  
      recv_mesg = [ 0 ]
     
      #println("Start sending out initial tasks") 
      while initial_task < MPI.Comm_size(comm)
        #println("Sending task $task to rank $initial_task")
        sreq = MPI.Send(task, initial_task, 1, comm)
        #println("Task $task sent to rank $initial_task") 
        
        task[1] -= batch_size 
        initial_task += 1 
      end
      #println("Done sending out intiial tasks") 

      #== hand out quartets to slaves dynamically ==#
      #println("Start sending out rest of tasks") 
      while task[1] > 0 
        status = MPI.Probe(MPI.MPI_ANY_SOURCE, 1, comm) 
        rreq = MPI.Recv!(recv_mesg, status.source, 1, comm)  
        #println("Sending task $task to rank ", status.source)
        sreq = MPI.Send(task, status.source, 1, comm)  
        #println("Task $task sent to rank ", status.source)
        task[1] -= batch_size 
      end
      #println("Done sending out rest of tasks") 
     
      #== hand out ending signals once done ==#
      #println("Start sending out enders") 
      for rank in 1:(MPI.Comm_size(comm)-1)
        #println("Sending ender to rank $rank")
        sreq = MPI.Send([ -1 ], rank, 0, comm)                           
        #println("Ender sent to rank $rank")
      end      
      #println("Done sending out enders") 
    #== slave ranks perform actual computations on quartets ==#
    elseif MPI.Comm_rank(comm) > 0
      #== intial setup ==#
      recv_mesg = [ 0 ]
      send_mesg = [ 0 ]

      ish_old = 0
      jsh_old = 0
      ksh_old = 0
      lsh_old = 0

      #mutex = Base.Threads.ReentrantLock()
      thread_index_counter = nindices
 
      #for thread in 1:Threads.nthreads()
      #  F_priv = zeros(basis.norb,basis.norb)

      max_shell_am = MAX_SHELL_AM
      eri_quartet_batch = Vector{Float64}(undef,81)

      quartet = ShQuartet(ShPair(basis.shells[1], basis.shells[1]),
        ShPair(basis.shells[1], basis.shells[1]))
   
      simint_workspace = Vector{Float64}(undef,10000)
      
      #== do computations ==# 
      while true 
        #== get shell quartet ==#
        status = MPI.Probe(0, MPI.MPI_ANY_TAG, comm)
        #println("About to recieve task from master")
        rreq = MPI.Recv!(recv_mesg, status.source, status.tag, comm)

        ijkl_index = recv_mesg[1]
        #println(ijkl_index)
        if ijkl_index < 0 break end
        #println("Recieved task $ijkl_index from master")
 
        #for rank in 1:MPI.Comm_size(comm)
        #  if MPI.Comm_rank(comm) == rank
        #    println("IJKL_INDEX: ", ijkl_index)
        #  end
        #end
        #println("NEW BATCH")
        for ijkl in ijkl_index:-1:(max(1,ijkl_index-batch_size+1))
          #println("IJKL: $ijkl")

         fock_build_thread_kernel(F, D,
            H, basis, eri_quartet_batch, mutex,
            quartet, ijkl, simint_workspace, schwarz_bounds, Dsh,
            ish_old, jsh_old, ksh_old, lsh_old; 
            debug=debug)
        end

        send_mesg[1] = MPI.Comm_rank(comm)
        MPI.Send(send_mesg, 0, 1, comm)
      #lock(mutex)
      #F .+= F_priv
      #unlock(mutex)
      end
    end
    MPI.Barrier(comm)
  end

  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      F[iorb,jorb] /= 2.0
      F[jorb,iorb] = F[iorb,jorb]
    end
  end

  return F
end

@inline function fock_build_thread_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  H::Matrix{Float64}, basis::BasisStructs.Basis, 
  eri_quartet_batch::Vector{Float64}, mutex, 
  quartet::ShQuartet, ijkl_index::Int64,
  simint_workspace::Vector{Float64}, schwarz_bounds::Matrix{Float64}, 
  Dsh::Matrix{Float64}, ish_old::Int64, jsh_old::Int64, ksh_old::Int64, 
  lsh_old::Int64; debug)

  comm=MPI.COMM_WORLD
  
  #== determine shells==# 
  bra_pair = decompose(ijkl_index)
  ket_pair = ijkl_index - triangular_index(bra_pair)

  ish = decompose(bra_pair)
  jsh = bra_pair - triangular_index(ish)

  ksh = decompose(ket_pair)
  lsh = ket_pair - triangular_index(ksh)
  
  #== create shell quartet ==#
  quartet.bra.sh_a = basis[ish]
  quartet.bra.sh_b = basis[jsh]
  quartet.ket.sh_a = basis[ksh]
  quartet.ket.sh_b = basis[lsh]

  #== Cauchy-Schwarz screening ==#
  bound = schwarz_bounds[ish, jsh]*schwarz_bounds[ksh, lsh] 

  dijmax = 4.0*Dsh[triangular_index(ish, jsh)]
  dklmax = 4.0*Dsh[triangular_index(ksh, lsh)]
  
  dikmax = Dsh[triangular_index(ish, ksh)]
  dilmax = Dsh[triangular_index(ish, lsh)]
  djkmax = Dsh[triangular_index(jsh, ksh)]
  djlmax = Dsh[triangular_index(jsh, lsh)]
 
  maxden = max(dijmax, dklmax, dikmax, dilmax, djkmax, djlmax)
  #bound *= maxden

  #== fock build for significant shell quartets ==# 
  if abs(bound) >= 1.0E-10 
    #== compute electron repulsion integrals ==#
    compute_eris(ish, jsh, ksh, lsh, eri_quartet_batch, simint_workspace)

    #== contract ERIs into Fock matrix ==#
    contract_eris(F, D, eri_quartet_batch, quartet,
      ish, jsh, ksh, lsh, debug)
  end
    #if debug println("END TWO-ELECTRON INTEGRALS") end
end

@inline function compute_eris(ish::Int64, jsh::Int64, ksh::Int64,
  lsh::Int64, eri_quartet_batch::Vector{Float64},
  simint_workspace::Vector{Float64})

  #= actually compute integrals =#
  SIMINT.compute_eris(ish, jsh, ksh, lsh, eri_quartet_batch, 
    simint_workspace)
end


@inline function contract_eris(F_priv::Matrix{Float64}, D::Matrix{Float64},
  eri_batch::Vector{Float64}, quartet::ShQuartet, ish::Int64, jsh::Int64,
  ksh::Int64, lsh::Int64, debug::Bool)

  norb = size(D,1)

  two_same = (ish == jsh) || (ish == ksh) || (ish == lsh) || (jsh == ksh) || 
    (jsh == lsh) || (ksh == lsh)

  three_same = (ish == jsh && jsh == ksh) || (ish == jsh && jsh == lsh) || 
    (ish == ksh && ksh == lsh) || (jsh == ksh && ksh == lsh)

  four_same = ish == jsh && jsh == ksh && ksh == lsh

  pμ = quartet.bra.sh_a.pos
  nμ = quartet.bra.sh_a.nbas

  pν = quartet.bra.sh_b.pos
  nν = quartet.bra.sh_b.nbas
  
  pλ = quartet.ket.sh_a.pos
  nλ = quartet.ket.sh_a.nbas
  
  pσ = quartet.ket.sh_b.pos
  nσ = quartet.ket.sh_b.nbas

  μνλσ = 0
  for μsize::Int64 in 0:(nμ-1), νsize::Int64 in 0:(nν-1)
    μμ = μsize + pμ
    νν = νsize + pν

    do_continue_bra = sort_bra(μμ, νν,
      ish, jsh, ksh, lsh, nμ, nν, nλ, nσ, two_same, three_same, four_same)
    if do_continue_bra 
      continue 
    end

    μνλσ = nσ*nλ*νsize + nσ*nλ*nν*μsize
    for λsize::Int64 in 0:(nλ-1), σsize::Int64 in 0:(nσ-1)
      λλ = λsize + pλ
      σσ = σsize + pσ

      #if debug
        #if do_continue_print print("$μμ, $νν, $λλ, $σσ => ") end
      #end

      #μνλσ = 1 + σsize + nσ*λsize + nσ*nλ*νsize + nσ*nλ*nν*μsize
      μνλσ += 1 
  
      eri = eri_batch[μνλσ] 
      
      do_continue_screen = abs(eri) < 1.0E-10
      if do_continue_screen 
        #if do_continue_print println("CONTINUE") end
        continue 
      end

      do_continue_ket = sort_ket(μμ, νν, λλ, σσ,
        ish, jsh, ksh, lsh, nμ, nν, nλ, nσ, two_same, three_same, four_same)
      if do_continue_ket 
        #if do_continue_print println("CONTINUE") end
        continue 
      end
  
      μ, ν = (μμ > νν) ? (μμ, νν) : (νν, μμ)
      λ, σ = (λλ > σσ) ? (λλ, σσ) : (σσ, λλ)

      do_continue_braket, μ, ν, λ, σ = sort_braket(μ, ν, λ, σ, ish, jsh, ksh, 
        lsh, nμ, nν, nλ, nσ)
      if do_continue_braket 
        #if do_continue_print println("CONTINUE") end
        continue 
      end

      #if debug println("ERI($μ, $ν, $λ, $σ) = $eri") end
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
  F_evec::Matrix{Float64}, F_mo::Matrix{Float64}, F_part::Matrix{Float64}, 
  ortho::Matrix{Float64}, ortho_trans::Matrix{Float64},
  basis::BasisStructs.Basis, iter::Int, debug::Bool)

  comm=MPI.COMM_WORLD
 
  #== obtain new orbital coefficients ==#
  BLAS.symm!('L', 'U', 1.0, ortho_trans, F_μν, 0.0, F_part)
  BLAS.gemm!('N', 'N', 1.0, F_part, ortho, 0.0, F_mo)
  
  F_eval .= eigvals(LinearAlgebra.Hermitian(F_mo))

  F_evec .= eigvecs(LinearAlgebra.Hermitian(F_mo))
  #@views F_evec .= F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

  #C .= ortho*F_evec
  BLAS.symm!('L', 'U', 1.0, ortho, F_evec, 0.0, C)
  
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-$iter/C", C)
  end

  #== build new density matrix ==#
  nocc = basis.nels >> 1
  norb = basis.norb

  fill!(D, 0.0)
  for i in 1:basis.norb, j in 1:basis.norb
    #@views D[i,j] = @∑ C[i,1:nocc] C[j,1:nocc]
    for iocc in 1:nocc
      D[i,j] += C[i, iocc] * C[j, iocc]
    end
    #D[i,j] = @∑ C[1:nocc,i] C[1:nocc,j]
  end
  D .*= 2.0
 
  #== compute new SCF energy ==#
  EHF1 = @∑ D F_μν
  EHF2 = @∑ D H
  E_elec = (EHF1 + EHF2)/2.0
  
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-$iter/D", D)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF1", EHF1)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF2", EHF2)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF", E_elec)
  end

  return E_elec
end
