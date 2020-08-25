using Base.Threads
using LinearAlgebra
using HDF5
using PrettyTables

const do_continue_print = false 
const print_eri = false 

function rhf_energy(mol::Molecule, basis::Basis,
  scf_flags::Union{Dict{String,Any},Dict{Any,Any}}; output)
  
  debug::Bool = scf_flags["debug"]
  niter::Int = scf_flags["niter"]

  ndiis::Int = scf_flags["ndiis"]
  dele::Float64 = scf_flags["dele"]
  rmsd::Float64 = scf_flags["rmsd"]
  load::String = scf_flags["load"]
  fdiff::Bool = scf_flags["fdiff"]

  return rhf_kernel(mol,basis; output=output, debug=debug, 
    niter=niter, ndiis=ndiis, dele=dele, rmsd=rmsd, load=load, fdiff=fdiff)
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
function rhf_kernel(mol::Molecule, 
  basis::Basis; 
  output::String, debug::Bool, niter::Int, ndiis::Int, 
  dele::Float64, rmsd::Float64, load::String, fdiff::Bool)

  comm=MPI.COMM_WORLD
  calculation_status = Dict([])

  #== read in some variables from scf input ==#
  debug_output = debug ? h5open("debug.h5","w") : nothing

  #== compute nuclear repulsion energy ==# 
  E_nuc = compute_enuc(mol)
  
  jeri_oei_engine = JERI.OEIEngine(mol.mol_cxx, basis.basis_cxx) 
  
  #== compute one-electron integrals and Hamiltonian ==#
  S = zeros(Float64, (basis.norb, basis.norb))
  compute_overlap(S, basis, jeri_oei_engine)
 
  #for i in 1:basis.norb, j in 1:i
  #  println("OVR($i,$j): ", S[i,j])
  #end
  
  T = zeros(Float64, (basis.norb, basis.norb))
  compute_ke(T, basis, jeri_oei_engine)
 
  V = zeros(Float64, (basis.norb, basis.norb))
  compute_nah(V, mol, basis, jeri_oei_engine)

  H = T .+ V
  
  #for i in 1:basis.norb, j in 1:i
  #  println("HAMIL($i,$j): ", H[i,j])
  #end
 
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-None/E_nuc", E_nuc)
    h5write("debug.h5","SCF/Iteration-None/S", S)
    h5write("debug.h5","SCF/Iteration-None/T", T)
    h5write("debug.h5","SCF/Iteration-None/V", V)
    h5write("debug.h5","SCF/Iteration-None/H", H)
  end

  #== build the initial matrices ==#
  F = deepcopy(H)
  F_eval = Vector{Float64}(undef,basis.norb)
  F_evec = similar(F)

  D = similar(F)
  C = similar(F)

  #== allocate workspace matrices ==#
  workspace_a = similar(F)
  workspace_b = similar(F)
  workspace_c = [ similar(F) ]

  #== build the orthogonalization matrix ==#
  workspace_b .= S
  S_eval_diag, workspace_a = eigen!(LinearAlgebra.Hermitian(workspace_b))

  fill!(workspace_b, 0.0)
  for i in 1:basis.norb
    workspace_b[i,i] = S_eval_diag[i]
  end
  
  ortho = workspace_a*(LinearAlgebra.Diagonal(workspace_b)^-0.5)*transpose(workspace_a)
  
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-None/Ortho", ortho)
  end

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("----------------------------------------          ")
    println("       Starting RHF iterations...                 ")
    println("----------------------------------------          ")
    println(" ")
    println("Iter      Energy                   ΔE                   Drms")
  end

  E_elec = 0.0
  E_elec = iteration(F, D, C, H, F_eval, F_evec, workspace_a, 
    workspace_b, ortho, basis, 0, debug)
  
  F_old = deepcopy(F)
  
  E = E_elec + E_nuc
  E_old = E

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println(0,"     ", E)
  end

  #=============================#
  #== start scf cycles: #7-10 ==#
  #=============================#
  F, D, C, E, converged = scf_cycles(F, D, C, E, H, ortho, S, 
    F_eval, F_evec, F_old, workspace_a, workspace_b, workspace_c,
    E_nuc, E_elec, E_old, basis; 
    output=output, debug=debug, niter=niter, ndiis=ndiis, dele=dele,
    rmsd=rmsd, load=load, fdiff=fdiff)

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
  S::Matrix{Float64}, F_eval::Vector{Float64}, 
  F_evec::Matrix{Float64},  F_old::Matrix{Float64},
  workspace_a::Matrix{Float64}, workspace_b::Matrix{Float64}, 
  workspace_c::Vector{Matrix{Float64}}, E_nuc::Float64, E_elec::Float64, 
  E_old::Float64, basis::Basis;
  output::String, debug::Bool, niter::Int, ndiis::Int, 
  dele::Float64, rmsd::Float64, load::String, fdiff::Bool)

  #== read in some more variables from scf flags input ==#
  nsh = length(basis)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3

  #== build DIIS arrays ==#
  F_array = fill(similar(F), ndiis)

  e_array = fill(similar(F), ndiis)
  e_array_old = fill(similar(F), ndiis-1)
  
  F_array_old = fill(similar(F), ndiis-1)

  #FD = similar(F)
  FDS = similar(F)
  #SDF = similar(F)
  
  #== build arrays needed for post-fock build iteration calculations ==#
  #F_temp = similar(F)
  ΔF = similar(F) 
  F_cumul = zeros(size(F)) 
 
  D_old = similar(F)
  ΔD = deepcopy(D) 
  D_input = similar(F)

  #== build matrix of Cauchy-Schwarz upper bounds ==# 
  schwarz_bounds = zeros(Float64,(nsh,nsh)) 
  compute_schwarz_bounds(schwarz_bounds, basis, nsh)

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

  E = scf_cycles_kernel(F, D, C, E, H, ortho, S, E_nuc,
    E_elec, E_old, basis, F_array, e_array, e_array_old,
    F_array_old, F_eval, F_evec, F_old, workspace_a, 
    workspace_b, workspace_c, ΔF, F_cumul, 
    D_old, ΔD, D_input, scf_converged, FDS, 
    schwarz_bounds, Dsh; 
    output=output, debug=debug, niter=niter, ndiis=ndiis, dele=dele, 
    rmsd=rmsd, load=load, fdiff=fdiff)

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
  C::Matrix{Float64}, E::Float64, H::Matrix{Float64}, 
  ortho::Matrix{Float64}, S::Matrix{Float64}, E_nuc::Float64, 
  E_elec::Float64, E_old::Float64, basis::Basis,
  F_array::Vector{Matrix{Float64}}, 
  e_array::Vector{Matrix{Float64}}, e_array_old::Vector{Matrix{Float64}},
  F_array_old::Vector{Matrix{Float64}}, 
  F_eval::Vector{Float64}, F_evec::Matrix{Float64}, 
  F_old::Matrix{Float64}, workspace_a::Matrix{Float64}, 
  workspace_b::Matrix{Float64}, workspace_c::Vector{Matrix{Float64}}, 
  ΔF::Matrix{Float64},
  F_cumul::Matrix{Float64}, D_old::Matrix{Float64}, 
  ΔD::Matrix{Float64}, D_input::Matrix{Float64}, scf_converged::Bool,  
  FDS::Matrix{Float64}, 
  schwarz_bounds::Matrix{Float64}, Dsh::Matrix{Float64}; 
  output, debug, niter, ndiis, dele, rmsd, load, fdiff)

  #== initialize a few more variables ==#
  comm=MPI.COMM_WORLD

  B_dim = 1
  D_rms = 1.0
  ΔE = 1.0 
  cutoff = fdiff ? 5E-11 : 1E-10

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

    #== determine input D and F ==#
    D_input .= fdiff ? ΔD : D
    workspace_b .= fdiff ? ΔF : F

    #== compress D into shells in Dsh ==#
    for ish in 1:length(basis), jsh in 1:ish
      ipos = basis[ish].pos
      ibas = basis[ish].nbas

      jpos = basis[jsh].pos
      jbas = basis[jsh].nbas
      
      max_value = 0.0
      for i in ipos:(ipos+ibas-1), j in jpos:(jpos+jbas-1) 
        max_value = max(max_value, abs(D_input[i,j]))
      end
      Dsh[ish, jsh] = max_value
      Dsh[jsh, ish] = Dsh[ish, jsh] 
    end
  
    #== build new Fock matrix ==#
    workspace_a .= fock_build(workspace_b, D_input, H, basis, schwarz_bounds, Dsh, 
      cutoff, debug, load)

    workspace_b .= MPI.Allreduce(workspace_a,MPI.SUM,comm)
    MPI.Barrier(comm)

    if debug && MPI.Comm_rank(comm) == 0
      h5write("debug.h5","SCF/Iteration-$iter/F/Skeleton", F_input)
    end
 
    if fdiff 
      ΔF .= workspace_b
      F_cumul .+= ΔF
      F .= F_cumul .+ H
    else
      F .= workspace_b .+ H
    end

    if debug && MPI.Comm_rank(comm) == 0
      h5write("debug.h5","SCF/Iteration-$iter/F/Total", F)
    end

    #== do DIIS ==#
    if ndiis > 0
      BLAS.symm!('L', 'U', 1.0, F, D, 0.0, workspace_a)
      BLAS.gemm!('N', 'N', 1.0, workspace_a, S, 0.0, FDS)
      
      transpose!(workspace_b, FDS)
      
      workspace_a .= FDS .- workspace_b #error matrix 

      e_array_old = view(e_array,1:(ndiis-1))                                   
      workspace_c[1] = deepcopy(workspace_a)
      e_array = vcat(workspace_c, e_array_old)                                                                          
      F_array_old = view(F_array,1:(ndiis-1))                                   
      workspace_c[1] = deepcopy(F)
      F_array = vcat(workspace_c, F_array_old)              
      
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

    E_elec = iteration(F, D, C, H, F_eval, F_evec, workspace_a,
      workspace_b, ortho, basis, iter, debug)

    #== check for convergence ==#
    ΔD .= D .- D_old
    D_rms = √(LinearAlgebra.dot(ΔD,ΔD))

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
  H::Matrix{Float64}, basis::Basis, 
  schwarz_bounds::Matrix{Float64}, Dsh::Matrix{Float64},
  cutoff::Float64, debug::Bool, load::String)

  comm = MPI.COMM_WORLD
  
  fill!(F,zero(Float64))

  nsh = length(basis)
  nindices = (nsh*(nsh+1)*(nsh^2 + nsh + 2)) >> 3 #bitwise divide by 8
 
  #mutex = Base.Threads.ReentrantLock()
  #thread_index_counter = Threads.Atomic{Int64}(nindices)
  
  #== simply do calculation for serial runs ==#
  if MPI.Comm_size(comm) == 1  || load == "static"
    top_index = nindices - (MPI.Comm_rank(comm))
    stride = MPI.Comm_size(comm) 
    
    mutex = Base.Threads.ReentrantLock()
    thread_index_counter = Threads.Atomic{Int64}(top_index)
    Threads.@threads for thread in 1:Threads.nthreads() 
      max_am = max_ang_mom(basis) 
      eri_quartet_batch_priv = Vector{Float64}(undef,eri_quartet_batch_size(max_am))
      simint_workspace_priv = Vector{Float64}(undef,get_workmem(0,max_am-1))
    
      F_priv = zeros(size(F))
      while true 
        ijkl = Threads.atomic_sub!(thread_index_counter, stride) 

        if ijkl < 1 break end
        
        fock_build_thread_kernel(F_priv, D,
          H, basis, eri_quartet_batch_priv, #mutex,
          ijkl, simint_workspace_priv, schwarz_bounds, Dsh,
          cutoff, debug)
      end

      lock(mutex)
        F .+= F_priv
      unlock(mutex)
    end
    #MPI.Barrier()
  #== use static task distribution for multirank runs if selected ==#
  elseif MPI.Comm_size(comm) > 1 && load == "dynamic"
    batch_size = ceil(Int,nindices/(MPI.Comm_size(comm)*
      Threads.nthreads()*1000)) 

    #== master rank ==#
    if MPI.Comm_rank(comm) == 0 
      #== send out initial tasks to slaves ==#
      task = [ nindices ]
      initial_task = 1
  
      recv_mesg = [ 0 ]
     
      #println("Start sending out initial tasks") 
      while initial_task < MPI.Comm_size(comm)
        for thread in 1:Threads.nthreads()
          #println("Sending task $task to rank $initial_task")
          sreq = MPI.Send(task, initial_task, thread, comm)
          #println("Task $task sent to rank $initial_task") 
        
          task[1] -= batch_size 
        end
        initial_task += 1
      end
      #println("Done sending out intiial tasks") 

      #== hand out quartets to slaves dynamically ==#
      #println("Start sending out rest of tasks") 
      while task[1] > 0 
        status = MPI.Probe(MPI.MPI_ANY_SOURCE, MPI.MPI_ANY_TAG, 
          comm) 
        rreq = MPI.Recv!(recv_mesg, status.source, status.tag, 
          comm)  
        #println("Sending task $task to rank ", status.source)
        sreq = MPI.Send(task, status.source, status.tag, comm)  
        #println("Task $task sent to rank ", status.source)
        task[1] -= batch_size 
      end
      #println("Done sending out rest of tasks") 
     
      #== hand out ending signals once done ==#
      #println("Start sending out enders") 
      for rank in 1:(MPI.Comm_size(comm)-1)
        for thread in 1:Threads.nthreads()
          sreq = MPI.Send([ -1 ], rank, thread, comm)                           
        end
      end      
      #println("Done sending out enders") 
    #== slave ranks perform actual computations on quartets ==#
    elseif MPI.Comm_rank(comm) > 0
      mutex_mpi = Base.Threads.ReentrantLock()
      mutex_reduce = Base.Threads.ReentrantLock()
      
      Threads.@threads for thread in 1:Threads.nthreads() 
        recv_mesg = [ 0 ]
        send_mesg = [ 0 ]

        max_am = max_ang_mom(basis) 
        eri_quartet_batch_priv = Vector{Float64}(undef,eri_quartet_batch_size(max_am))
        simint_workspace_priv = Vector{Float64}(undef,get_workmem(0,max_am-1))
    
        F_priv = zeros(size(F))
        while true 
          #== get shell quartet ==#
          #status = MPI.Probe(0, MPI.MPI_ANY_TAG, comm)
          #println("About to recieve task from master")
      
          lock(mutex_mpi)
            status = MPI.Probe(0, thread, comm)
            rreq = MPI.Recv!(recv_mesg, status.source, status.tag, comm)
            ijkl_index = recv_mesg[1]
          unlock(mutex_mpi)

          #println(ijkl_index)
          if ijkl_index < 0 break end
          #println("Thread $thread ecieved task $ijkl_index from master")
 
          #for rank in 1:MPI.Comm_size(comm)
          #  if MPI.Comm_rank(comm) == rank
          #    println("IJKL_INDEX: ", ijkl_index)
          #  end
          #end
          #println("NEW BATCH")
          for ijkl in ijkl_index:-1:(max(1,ijkl_index-batch_size+1))
            #println("IJKL: $ijkl")

            fock_build_thread_kernel(F_priv, D,
              H, basis, eri_quartet_batch_priv, #mutex,
              ijkl, simint_workspace_priv, schwarz_bounds, Dsh,
              cutoff, debug)
          end

          lock(mutex_mpi)
            send_mesg[1] = MPI.Comm_rank(comm)
            MPI.Send(send_mesg, 0, thread, comm)
          unlock(mutex_mpi)
        end
      
        lock(mutex_reduce)
          F .+= F_priv
        unlock(mutex_reduce)
      end
    end
    #MPI.Barrier(comm)
  end
  MPI.Barrier(comm)

  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      F[iorb,jorb] /= 2.0
      F[jorb,iorb] = F[iorb,jorb]
    end
  end

  return F
end

@inline function fock_build_thread_kernel(F::Matrix{Float64}, D::Matrix{Float64},
  H::Matrix{Float64}, basis::Basis, 
  eri_quartet_batch::Vector{Float64}, 
  ijkl_index::Int64,
  simint_workspace::Vector{Float64}, schwarz_bounds::Matrix{Float64}, 
  Dsh::Matrix{Float64}, cutoff::Float64, debug::Bool)

  comm=MPI.COMM_WORLD
  
  #== determine shells==# 
  bra_pair = decompose(ijkl_index)
  ket_pair = ijkl_index - triangular_index(bra_pair)

  #quartet.bra = basis.shpair_ordering[bra_pair]
  #quartet.ket = basis.shpair_ordering[ket_pair]
 
  #ish = μsh.shell_id 
  #jsh = νsh.shell_id 
  #ksh = λsh.shell_id 
  #lsh = σsh.shell_id 
  
  ish = decompose(bra_pair)
  jsh = bra_pair - triangular_index(ish)

  ksh = decompose(ket_pair)
  lsh = ket_pair - triangular_index(ksh)

  μsh = basis[ish] 
  νsh = basis[jsh] 
  λsh = basis[ksh] 
  σsh = basis[lsh] 
  
  #icls = unsafe_string(μsh.class)
  #jcls = unsafe_string(νsh.class) 
  #kcls = unsafe_string(λsh.class) 
  #lcls = unsafe_string(σsh.class)

  #println("QUARTET($ish, $jsh, $ksh, $lsh) -> ($icls $jcls | $kcls $lcls)")

  #== Cauchy-Schwarz screening ==#
  bound = schwarz_bounds[ish, jsh]*schwarz_bounds[ksh, lsh] 

  dijmax = 4.0*Dsh[ish, jsh]
  dklmax = 4.0*Dsh[ksh, lsh]
  
  dikmax = Dsh[ish, ksh]
  dilmax = Dsh[ish, lsh]
  djkmax = Dsh[jsh, ksh]
  djlmax = Dsh[jsh, lsh]
 
  maxden = max(dijmax, dklmax, dikmax, dilmax, djkmax, djlmax)
  bound *= maxden

  #== fock build for significant shell quartets ==# 
  if abs(bound) >= cutoff 
    #== compute electron repulsion integrals ==#
    compute_eris(ish, jsh, ksh, lsh, μsh, νsh, λsh, σsh,
      eri_quartet_batch, simint_workspace)

    #== contract ERIs into Fock matrix ==#
    contract_eris(F, D, eri_quartet_batch, ish, jsh, ksh, lsh,
      μsh, νsh, λsh, σsh, debug)
  end
    #if debug println("END TWO-ELECTRON INTEGRALS") end
end

@inline function compute_eris(ish::Int64, jsh::Int64, ksh::Int64, lsh::Int64, 
  μsh::JCModules.Shell, νsh::JCModules.Shell, 
  λsh::JCModules.Shell, σsh::JCModules.Shell,
  eri_quartet_batch::Vector{Float64},
  simint_workspace::Vector{Float64})

  amμ = μsh.am
  amν = νsh.am
  amλ = λsh.am
  amσ = σsh.am

  #fill!(eri_quartet_batch, 0.0)
  #ish = μsh.shell_id
  #jsh = νsh.shell_id
  #ksh = λsh.shell_id
  #lsh = σsh.shell_id

  #= actually compute integrals =#
  SIMINT.compute_eris(ish, jsh, ksh, lsh, eri_quartet_batch, 
    simint_workspace)

  amμ = μsh.am
  amν = νsh.am
  amλ = λsh.am
  amσ = σsh.am

  nμ = μsh.nbas
  nν = νsh.nbas
  nλ = λsh.nbas
  nσ = σsh.nbas

  μνλσ = 0 
  for μsize::Int64 in 0:(nμ-1), νsize::Int64 in 0:(nν-1)
    μνλσ = nσ*nλ*νsize + nσ*nλ*nν*μsize
      
    μnorm = axial_norm_fact[μsize+1,amμ]
    νnorm = axial_norm_fact[νsize+1,amν]

    μνnorm = μnorm*νnorm

    for λsize::Int64 in 0:(nλ-1), σsize::Int64 in 0:(nσ-1)
      μνλσ += 1 
   
      λnorm = axial_norm_fact[λsize+1,amλ]
      σnorm = axial_norm_fact[σsize+1,amσ]
    
      λσnorm = λnorm*σnorm 
      eri_quartet_batch[μνλσ] *= μνnorm*λσnorm
    end 
  end

  #=
  if am[1] == 3 || am[2] == 3 || am[3] == 3 || am[4] == 3
    for idx in 1:nμ*nν*nλ*nσ 
    #for idx in 1:1296
      eri = eri_quartet_batch[idx]
      println("QUARTET($ish, $jsh, $ksh, $lsh): $eri")
    end
  end
  =#
end


@inline function contract_eris(F_priv::Matrix{Float64}, D::Matrix{Float64},
  eri_batch::Vector{Float64}, ish::Int64, jsh::Int64,
  ksh::Int64, lsh::Int64, 
  μsh::JCModules.Shell, νsh::JCModules.Shell, 
  λsh::JCModules.Shell, σsh::JCModules.Shell,
  debug::Bool)

  norb = size(D,1)
  
  #ish = μsh.shell_id
  #jsh = νsh.shell_id
  #ksh = λsh.shell_id
  #lsh = σsh.shell_id

  pμ = μsh.pos
  nμ = μsh.nbas

  pν = νsh.pos
  nν = νsh.nbas
  
  pλ = λsh.pos
  nλ = λsh.nbas
  
  pσ = σsh.pos
  nσ = σsh.nbas

  #amμ = μsh.am
  #amν = νsh.am
  #amλ = λsh.am
  #amσ = σsh.am
  #am = [ amμ, amν, amλ, amσ ]

  μνλσ = 0
  for μsize::Int64 in 0:(nμ-1), νsize::Int64 in 0:(nν-1)
    μμ = μsize + pμ
    νν = νsize + pν

    if μμ < νν && ish == jsh 
      #if do_continue_print println("CONTINUE BRA: $μμ, $νν") end
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
   
      if abs(eri) < 1.0E-10
        #if do_continue_print println("CONTINUE SCREEN") end
        continue 
      end

      if λλ < σσ && ksh == lsh 
        #if do_continue_print println("CONTINUE KET") end
        continue 
      end
  
      μ, ν = (μμ > νν) ? (μμ, νν) : (νν, μμ)
      λ, σ = (λλ > σσ) ? (λλ, σσ) : (σσ, λλ)

      μν = triangular_index(μ,ν)                                                    
      λσ = triangular_index(λ,σ)                                                    
       
      if μν < λσ 
        if ish == ksh && jsh == lsh 
          #if do_continue_print println("CONTINUE BRAKET") end
          continue 
        else
          μ, ν, λ, σ = λ, σ, μ, ν
        end
      end

      #println("QUARTET($ish, $jsh, $ksh, $lsh): $eri")
      #println("ERI($μ, $ν, $λ, $σ) = $eri") 
      
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
  F_evec::Matrix{Float64}, workspace_a::Matrix{Float64}, 
  workspace_b::Matrix{Float64}, ortho::Matrix{Float64}, 
  basis::Basis, iter::Int, debug::Bool)

  comm=MPI.COMM_WORLD
 
  transpose!(workspace_b, LinearAlgebra.Hermitian(ortho)) 

  #== obtain new orbital coefficients ==#
  BLAS.symm!('L', 'U', 1.0, workspace_b, F_μν, 0.0, workspace_a)
  BLAS.gemm!('N', 'N', 1.0, workspace_a, ortho, 0.0, workspace_b)
 
  F_eval, F_evec = eigen!(LinearAlgebra.Hermitian(workspace_b)) 
  
  #@views F_evec .= F_evec[:,sortperm(F_eval)] #sort evecs according to sorted evals

  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-$iter/F_evec/Sorted", F_mo)
  end

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
    for iocc in 1:nocc
      D[i,j] += 2 * C[i, iocc] * C[j, iocc]
    end
  end
 
  #== compute new SCF energy ==#
  EHF1 = LinearAlgebra.dot(D, F_μν)
  EHF2 = LinearAlgebra.dot(D, H)
  E_elec = (EHF1 + EHF2)/2.0
  
  if debug && MPI.Comm_rank(comm) == 0
    h5write("debug.h5","SCF/Iteration-$iter/D", D)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF1", EHF1)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF2", EHF2)
    h5write("debug.h5","SCF/Iteration-$iter/E/EHF", E_elec)
  end

  return E_elec
end
