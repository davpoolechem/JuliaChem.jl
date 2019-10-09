using JCModules.BasisStructs
using JCModules.Globals

using Base.Threads
using MATH
using JLD

function sort_braket(μμ, t_μ, νν, t_ν, λλ, t_λ, σσ, t_σ, ish, jsh, ksh, lsh,
  nμ, nν, nλ, nσ)

  do_continue = false
  μ,ν,λ,σ = t_μ, t_ν, t_λ, t_σ

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

  condition1 = μμ == λλ && νν == σσ
  condition1 = condition1 || (μμ == νν && λλ == σσ)
  condition1 = condition1 || (μμ == σσ && λλ == νν)
  condition1 = condition1 &&
	nμ > 1 && nν > 1 && nλ > 1 && nσ > 1

  condition2 = μμ == νν && λλ == σσ &&
	ish == jsh && jsh == ksh && ksh == lsh

  condition3 = two_same && !(ish == ksh && jsh == lsh)
  condition3 = condition3 &&
	nμ > 1 && nν > 1 && nλ > 1 && nσ > 1

  condition4 =  nμ > 1 && nν > 1 && nλ > 1
  condition4 =  condition4 || (nμ > 1 && nν > 1 && nσ > 1)
  condition4 =  condition4 || (nμ > 1 && nλ > 1 && nσ > 1)
  condition4 =  condition4 || (nν > 1 && nλ > 1 && nσ > 1)

  condition5 = ish == ksh && jsh == lsh
  condition5 = condition5 &&
	nμ > 1 && nν == 1 && nλ > 1 && nσ == 1

  condition6 = ish == jsh
  condition6 = condition6 &&
	nμ > 1 && nν > 1 && (nλ == 1 || nσ == 1)

  condition7 = ksh == lsh
  condition7 = condition7 &&
	(nμ == 1 || nν == 1) && nλ > 1 && nσ > 1

  μ, ν = (μμ > νν) ? (μμ, νν) : (νν, μμ)
  if μμ < νν && condition1
	do_continue = true
  end
  μν = triangular_index(μμ,νν)

  λ,σ = (λλ > σσ) ? (λλ, σσ) : (σσ, λλ)
  if λλ < σσ && condition1
	do_continue = true
  end
  λσ = triangular_index(λλ,σσ)

  if μμ < λλ && νν < σσ && condition2
	do_continue = true
  elseif (μμ < νν || λλ < σσ) && condition3
	do_continue = true
  elseif (μμ < νν && λλ < σσ) && condition4
	do_continue = true
  elseif μμ < λλ && condition5
	do_continue = true
  elseif μμ < νν && condition6
	do_continue = true
  elseif λλ < σσ && condition7
	do_continue = true
  elseif μν < λσ
	two_shell = nμ == nν
	two_shell = two_shell || (nμ == nλ)
	two_shell = two_shell || (nμ == nσ)
    two_shell = two_shell || (nν == nλ)
    two_shell = two_shell || (nν == nσ)
    two_shell = two_shell || (nλ == nσ)

    three_shell = nμ == nν && nν == nλ
    three_shell = three_shell || (nμ == nν && nν == nσ)
    three_shell = three_shell || (nμ == nλ && nλ == nσ)
    three_shell = three_shell || (nν == nλ && nλ == nσ)

    four_shell = nμ == nν
    four_shell = four_shell && (nν == nλ)
    four_shell = four_shell && (nλ == nσ)

    if four_shell
      if ish == ksh && jsh == lsh
        do_continue = true
      end
    elseif three_shell
      if μμ < νν && λλ < σσ
        do_continue = true
      end
    #elseif two_shell
    end
	if !do_continue λ, σ, μ, ν = t_μ, t_ν, t_λ, t_σ end
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
@inline function triangular_index(a::Integer,b::Integer)
	if a < b a, b = b, a end
  index = (a*(a-1)) >> 1 #bitwise divide by 2
  index += b
  return index
end

@inline function triangular_index(a::Integer)
  return (a*(a-1)) >> 1
end

@noinline function decompose(input)
  return trunc(Integer,cld((-1.0+√(1+8*input)),2.0))
end

function read_in_enuc()
	return input_enuc()
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
function read_in_oei(oei, nbf::Integer)
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

function DIIS(e_array, F_array, B_dim)
  B = Matrix{Float64}(undef,B_dim+1,B_dim+1)
  for i in 1:B_dim, j in 1:B_dim
    B[i,j] = @∑ e_array[i] e_array[j]

	  B[i,B_dim+1] = -1
	  B[B_dim+1,i] = -1
	  B[B_dim+1,B_dim+1] =  0
  end
  DIIS_coeff = [ fill(0.0,B_dim)..., -1.0 ]

  DIIS_coeff[:,:], B[:,:], ipiv = LinearAlgebra.LAPACK.gesv!(B, DIIS_coeff)

  F_DIIS = zeros(size(F_array[1],1),size(F_array[1],2))
  for index in 1:B_dim
    F_DIIS[:,:] .+= DIIS_coeff[index]*F_array[index]
  end

  return F_DIIS
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
