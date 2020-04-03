using JCModules.BasisStructs
using JCModules.MolStructs
using JCModules.Globals

using Base.Threads
using MATH
using JLD

@inline function sort_bra(μμ::Int, νν::Int, ish::Int, jsh::Int, ksh::Int, 
  lsh::Int, nμ::Int, nν::Int, nλ::Int, nσ::Int, two_same::Bool, 
  three_same::Bool, four_same::Bool)

  do_continue = false

  condition3 = false

  #condition6 = ish == jsh && ((nμ > nλ && nν > nλ) || (nμ > nσ && nν > nσ)) 
  condition6 = ish == jsh 

  if μμ < νν && condition6 
	  do_continue = true
  end
  return do_continue
end

@inline function sort_ket(μμ::Int, νν::Int, λλ::Int, σσ::Int, ish::Int, 
  jsh::Int, ksh::Int, lsh::Int, nμ::Int, nν::Int, nλ::Int, nσ::Int, 
  two_same::Bool, three_same::Bool, four_same::Bool)

  do_continue = false

  #condition3 = two_same && !(ish == ksh && jsh == lsh) && nμ > 1 && nν > 1 && 
  #  nλ > 1 && nσ > 1
  condition3 = false

  condition5 = ish == ksh && jsh == lsh && nμ > nν && nλ > nσ 

  condition7 = ksh == lsh && ((nλ > nμ  && nσ > nμ) || (nλ > nν && nσ > nν)) 

  condition8 = four_same

  if μμ < λλ && condition5
	  do_continue = true
  elseif λλ < σσ && (condition3 || condition7 || condition8)
	  do_continue = true
  end

  return do_continue
end

@inline function sort_braket(μ::Int, ν::Int, λ::Int, σ::Int, ish::Int, 
  jsh::Int, ksh::Int, lsh::Int, nμ::Int, nν::Int, nλ::Int, nσ::Int)

  do_continue = false

  μν = triangular_index(μ,ν)
  λσ = triangular_index(λ,σ)

  if μν < λσ
    three_shell = (nμ == nν && nν == nλ) || (nμ == nν && nν == nσ) || 
      (nμ == nλ && nλ == nσ) || (nν == nλ && nλ == nσ)

    four_shell = nμ == nν && nν == nλ && nλ == nσ

    #condition1 = nμ > 1 && nν > 1 && nλ > 1 && nσ > 1
    condition1 = false
      
    if four_shell && ((ish == ksh && jsh == lsh)) 
      do_continue = true
    #elseif three_shell && (μ < ν || λ < σ)
    #  do_continue = true
    else
	    λ, σ, μ, ν = μ, ν, λ, σ 
    end
    #if !do_continue λ, σ, μ, ν = μ, ν, λ, σ end
  end

  return do_continue, μ, ν, λ, σ
end

#=
function set_up_eri_database(basis::BasisStructs.Basis)
  jldopen("tei_batch.jld", "w") do file
    #== write quartet eri lists to database ==#
    eri_array = load("tei_all.jld")["Integrals"]["All"]

    nsh::Int64 = length(basis.shells)

    eri_array_batch::Vector{Float64} = [ ]
    #eri_array_starts::Vector{Int64} = [ ]
	  eri_array_sizes::Vector{Int64} = [ ]

    eri_start::Int64 = 1
    quartet_batch_num_old::Int64 = 1

    for ish::Int64 in 1:nsh, jsh::Int64 in 1:ish
      ijsh::Int64 = index(ish,jsh)
      qnum_ij = ish*(ish-1)/2 + jsh

      ibas::Int64 = basis.shells[ish].nbas
      jbas::Int64 = basis.shells[jsh].nbas

      ipos::Int64 = basis.shells[ish].pos
      jpos::Int64 = basis.shells[jsh].pos

      for ksh::Int64 in 1:nsh, lsh::Int64 in 1:ksh
        klsh::Int64 = index(ksh,lsh)
        if (klsh > ijsh) continue end

        kbas::Int64 = basis.shells[ksh].nbas
        lbas::Int64 = basis.shells[lsh].nbas

        kpos::Int64 = basis.shells[ksh].pos
        lpos::Int64 = basis.shells[lsh].pos

        qnum_kl::Int64 = ksh*(ksh-1)/2 + lsh
        quartet_num::Int64 = qnum_ij*(qnum_ij-1)/2 + qnum_kl - 1

        qint_ij::Int64 = ibas*(ibas-1)/2 + jbas
        qint_kl::Int64 = kbas*(kbas-1)/2 + lbas

        eri_size::Int64 = 0
        #println("QUARTET: $ish, $jsh, $ksh, $lsh")
        for μμ::Int64 in ipos:ipos+(ibas-1), νν::Int64 in jpos:jpos+(jbas-1)
          μ::Int64, ν::Int64 = μμ,νν
          if (μμ < νν) continue end

          μν::Int64 = index(μμ,νν)

          for λλ::Int64 in kpos:kpos+(kbas-1), σσ::Int64 in lpos:lpos+(lbas-1)
            λ::Int64, σ::Int64 = λλ,σσ
            if (λλ < σσ) continue end

            λσ::Int64 = index(λλ,σσ)

            if (μν < λσ)
              do_continue::Bool = false

              #print("$μμ, $νν, $λλ, $σσ => ")
              do_continue, μ, ν, λ, σ = sort_braket(μμ, νν, λλ, σσ, ish, jsh,
                ksh, lsh, ibas, jbas, kbas, lbas)

              if (do_continue)
                continue
              end
            end
            #println(" $μ, $ν, $λ, $σ")
            eri_size += 1
          end
        end

        quartet_batch_num::Int64 = Int64(floor(quartet_num/
          QUARTET_BATCH_SIZE)) + 1

        if quartet_batch_num != quartet_batch_num_old

          #== write arrays to disk ==#
          write(file, "Integrals/$quartet_batch_num_old",
            eri_array_batch)
        #  write(file, "Starts/$quartet_batch_num_old",
          #  eri_array_starts)
		      write(file, "Sizes/$quartet_batch_num_old",
            eri_array_sizes)

          #== reset variables as needed ==#
          eri_array_batch = [ ]
      #    eri_array_starts = [ ]
		       eri_array_sizes = [ ]

          quartet_batch_num_old = quartet_batch_num
        end

        append!(eri_array_batch,
          @view eri_array[eri_start:eri_start+(eri_size-1)])

      #  eri_start_readin::Int64 = eri_start - QUARTET_BATCH_SIZE*
    #      (quartet_batch_num-1)
      #  push!(eri_array_starts,eri_start_readin)
        push!(eri_array_sizes,eri_size)

        eri_start += eri_size

        #println("$ish, $jsh, $ksh, $lsh, $quartet_num, $eri_size")
      end
    end

    write(file, "Integrals/$quartet_batch_num_old",
      eri_array_batch)
  #  write(file, "Starts/$quartet_batch_num_old",
    #  eri_array_starts)
  	write(file, "Sizes/$quartet_batch_num_old",
	    eri_array_sizes)
  end
end
=#

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
@inline function triangular_index(a::Int,b::Int)
  index = (a*(a-1)) >> 1 #bitwise divide by 2
  index += b
  return index
end

@inline function triangular_index(a::Int)
  return (a*(a-1)) >> 1
end

@inline function decompose(input::Int)
  return ceil(Int,(-1.0+√(1+8*input))/2.0)
  #return ccall((:decompose, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libjeri.so"),
  #  Int64, (Int64,), input)
end

function read_in_enuc()
	return input_enuc()
end

function compute_enuc(mol::MolStructs.Molecule)
  E_nuc = 0.0
  for iatom in 1:length(mol.atoms), jatom in 1:(iatom-1)
    ix = mol.atoms[iatom].atom_center[1] 
    jx = mol.atoms[jatom].atom_center[1] 

    iy = mol.atoms[iatom].atom_center[2] 
    jy = mol.atoms[jatom].atom_center[2] 

    iz = mol.atoms[iatom].atom_center[3]
    jz = mol.atoms[jatom].atom_center[3]
  
    distance = √((jx-ix)^2 + (jy-iy)^2 + (jz-iz)^2) 
    
    E_nuc += mol.atoms[iatom].atom_id*mol.atoms[jatom].atom_id/distance
  end 
  
  return E_nuc
end
 
function compute_overlap(S::Matrix{Float64}, basis::BasisStructs.Basis)
  for ash in 1:length(basis.shells), bsh in 1:ash
    abas = basis.shells[ash].nbas
    bbas = basis.shells[bsh].nbas
    
    apos = basis.shells[ash].pos
    bpos = basis.shells[bsh].pos
       
    S_block = zeros(abas*bbas)
    SIMINT.compute_overlap(ash, bsh, S_block)
    
    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      S[max(iorb,jorb),min(iorb,jorb)] = S_block[idx]
      
      idx += 1 
    end
  end
  
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      S[min(iorb,jorb),max(iorb,jorb)] = S[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

function compute_ke(T::Matrix{Float64}, basis::BasisStructs.Basis)
  for ash in 1:length(basis.shells), bsh in 1:ash
    abas = basis.shells[ash].nbas
    bbas = basis.shells[bsh].nbas
    
    apos = basis.shells[ash].pos
    bpos = basis.shells[bsh].pos
       
    T_block = zeros(Float64, (abas*bbas,))
    SIMINT.compute_ke(ash, bsh, T_block)
    
    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      T[max(iorb,jorb),min(iorb,jorb)] = T_block[idx]
      
      idx += 1 
    end
  end
  
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      T[min(iorb,jorb),max(iorb,jorb)] = T[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

function compute_nah(V::Matrix{Float64}, mol::MolStructs.Molecule, 
  basis::BasisStructs.Basis)
  
  #== define ncenter ==#
  ncenter::Int64 = length(mol.atoms)
  
  Z = Vector{Float64}([])
  x = Vector{Float64}([])
  y = Vector{Float64}([])
  z = Vector{Float64}([])

  for atom in mol.atoms 
    push!(Z, convert(Float64,atom.atom_id))  
    push!(x, atom.atom_center[1])  
    push!(y, atom.atom_center[2])  
    push!(z, atom.atom_center[3])  
  end

  for ash in 1:length(basis.shells), bsh in 1:ash
    abas = basis.shells[ash].nbas
    bbas = basis.shells[bsh].nbas
    
    apos = basis.shells[ash].pos
    bpos = basis.shells[bsh].pos
       
    V_block = zeros(Float64, (abas*bbas,))
    SIMINT.compute_nah(ncenter, Z, x, y, z, ash, bsh, V_block)
    
    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      V[max(iorb,jorb),min(iorb,jorb)] = V_block[idx]
      
      idx += 1 
    end
  end
  
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      V[min(iorb,jorb),max(iorb,jorb)] = V[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

#=
"""
		get_oei_matrix(oei::Array{Float64,2})
Summary
======
Extract one-electron integrals from data file object. Kinetic energy integrals,
overlap integrals, and nuclear attraction integrals can all be extracted.

Arguments
======
oei = array of one-electron integrals to extract
"""
=#
function read_in_oei(oei::Vector{T}, nbf::Int) where T
	nbf2 = (nbf*(nbf+1)) >> 1

	oei_matrix = Matrix{Float64}(undef,(nbf,nbf))
	for ibf in 1:nbf2
    i = decompose(ibf)
    j = ibf - triangular_index(i)

		oei_matrix[i,j] = float(oei[ibf])
		oei_matrix[j,i] = oei_matrix[i,j]
	end

	return oei_matrix
end

function DIIS(F::Matrix{Float64}, e_array::Vector{Matrix{Float64}}, 
  F_array::Vector{Matrix{Float64}}, B_dim::Int64)
  
  B = Matrix{Float64}(undef,B_dim+1,B_dim+1)
  for i in 1:B_dim, j in 1:B_dim
    B[i,j] = @∑ e_array[i] e_array[j]

	  B[i,B_dim+1] = -1
	  B[B_dim+1,i] = -1
	  B[B_dim+1,B_dim+1] =  0
  end
  #DIIS_coeff::Vector{Float64} = [ fill(0.0,B_dim)..., -1.0 ]
  DIIS_coeff::Vector{Float64} = vcat(zeros(B_dim), [-1.0])

  #DIIS_coeff[:], B[:,:], ipiv = LinearAlgebra.LAPACK.gesv!(B, DIIS_coeff)
  DIIS_coeff[:], B[:,:], ipiv = LinearAlgebra.LAPACK.sysv!('U', B, DIIS_coeff)
  
  fill!(F, zero(Float64))
  for index in 1:B_dim
    F .+= DIIS_coeff[index] .* F_array[index]
  end
end

macro eri_quartet_batch_size(max_am)
  return quote
    if $(max_am) == "s"
      1
    elseif $(max_am) == "p"
      81
    elseif $(max_am) == "L"
      256
    elseif $(max_am) == "d"
      1296
    elseif $(max_am) == "f"
      10000
    else throw
    end
  end
end
