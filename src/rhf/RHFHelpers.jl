using JCModules.BasisStructs

using Base.Threads
using MATH
using JLD

function sort_braket(μμ::Int64, νν::Int64, λλ::Int64, σσ::Int64,
  ish::Int64, jsh::Int64, ksh::Int64, lsh::Int64,
  ibas::Int64, jbas::Int64, kbas::Int64, lbas::Int64)

  do_continue::Bool = false
  μ,ν,λ,σ = μμ,νν,λλ,σσ

  #print("$μμ, $νν, $λλ, $σσ => ")

  two_shell::Bool = ibas == jbas
  two_shell = two_shell || (ibas == kbas)
  two_shell = two_shell || (ibas == lbas)
  two_shell = two_shell || (jbas == kbas)
  two_shell = two_shell || (jbas == lbas)
  two_shell = two_shell || (kbas == lbas)

  three_shell::Bool = ibas == jbas && jbas == kbas
  three_shell = three_shell || (ibas == jbas && jbas == lbas)
  three_shell = three_shell || (ibas == kbas && kbas == lbas)
  three_shell = three_shell || (jbas == kbas && kbas == lbas)

  four_shell::Bool = ibas == jbas
  four_shell = four_shell && (jbas == kbas)
  four_shell = four_shell && (kbas == lbas)

  if four_shell
    three_same::Bool = ish == jsh && jsh == ksh
    three_same = three_same || (ish == jsh && jsh == lsh)
    three_same = three_same || (ish == ksh && ksh == lsh)
    three_same = three_same || (jsh == ksh && ksh == lsh)

    four_same::Bool = ish == jsh
    four_same = four_same && jsh == ksh
    four_same = four_same && ksh == lsh

    if four_same
      #print("\n")
      do_continue = true
    elseif three_same
      if (μμ != λλ && νν != σσ)
        μ,ν,λ,σ = λλ,σσ,μμ,νν
      else
        #print("\n")
        do_continue = true
      end
    else
      #print("\n")
      do_continue = true
    end
  elseif three_shell
    if (μμ != λλ && νν != σσ)
      μ,ν,λ,σ = λλ,σσ,μμ,νν
    else
      #print("\n")
      do_continue = true
    end
  elseif two_shell
    if (ish == ksh && jsh == lsh &&
      μμ != νν && μμ != λλ && μμ != σσ &&
      νν != λλ && νν != σσ &&
      λλ != σσ )
      #print("\n")
      do_continue = true
    elseif (μμ != λλ && νν != σσ)
      μ,ν,λ,σ = λλ,σσ,μμ,νν
    else
      #print("\n")
      do_continue = true
    end
  end

  return do_continue, μ, ν, λ, σ
end

function set_up_eri_database(basis::BasisStructs.Basis)
  jldopen("tei_batch.jld", "w") do file
    #== write quartet eri lists to database ==#
    eri_array = load("tei_all.jld")["Integrals"]["All"]

    nsh::Int64 = length(basis.shells)

    eri_array_batch::Vector{Float64} = [ ]
    eri_array_starts::Vector{Int64} = [ ]
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

        quartets_per_batch::Int64 = 10000
        quartet_batch_num::Int64 = Int64(floor(quartet_num/
          quartets_per_batch)) + 1

        if quartet_batch_num != quartet_batch_num_old

          #== write arrays to disk ==#
          write(file, "Integrals/$quartet_batch_num_old",
            eri_array_batch)
          write(file, "Starts/$quartet_batch_num_old",
            eri_array_starts)
          write(file, "Sizes/$quartet_batch_num_old",
            eri_array_sizes)

          #== reset variables as needed ==#
          eri_array_batch = [ ]
          eri_array_starts = [ ]
          eri_array_sizes = [ ]

          quartet_batch_num_old = quartet_batch_num
        end

        append!(eri_array_batch,
          @view eri_array[eri_start:eri_start+(eri_size-1)])

        eri_start_readin::Int64 = eri_start - quartets_per_batch*
          (quartet_batch_num-1)
        push!(eri_array_starts,eri_start_readin)

        push!(eri_array_sizes,eri_size)

        eri_start += eri_size

        #println("$ish, $jsh, $ksh, $lsh, $quartet_num, $eri_size")
      end
    end

    write(file, "Integrals/$quartet_batch_num_old",
      eri_array_batch)
    write(file, "Starts/$quartet_batch_num_old",
      eri_array_starts)
    write(file, "Sizes/$quartet_batch_num_old",
      eri_array_sizes)
  end
end

function read_in_enuc()
	enuc::Float64 = input_enuc()

	return enuc
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
function read_in_oei(oei::Vector{T}, nbf::Int64) where {T}
	nbf2::Int64 = nbf*(nbf+1)/2

    ioff::Vector{Int64} = map((x) -> x*(x-1)/2, collect(1:nbf*(nbf+1)))

	oei_matrix::Matrix{T} = Matrix{T}(undef,(nbf,nbf))
	Threads.@threads for index::Int64 in 1:nbf2
        i::Int64 = ceil(((-1+sqrt(1+8*index))/2))
        j::Int64 = index - ioff[i]

        eri = oei[index]
		oei_matrix[i,j] = oei[index]
		oei_matrix[j,i] = oei_matrix[i,j]
	end

	return oei_matrix
end

function DIIS(e_array::Vector{Matrix{T}},
  F_array::Vector{Matrix{T}}, B_dim::Int64) where {T<:AbstractFloat}

  B::Matrix{T} = Matrix{T}(undef,B_dim+1,B_dim+1)
  for i::Int64 in 1:B_dim, j::Int64 in 1:B_dim
    B[i,j] = @∑ e_array[i] e_array[j]

	B[i,B_dim+1] = -1
	B[B_dim+1,i] = -1
	B[B_dim+1,B_dim+1] =  0
  end
  DIIS_coeff::Vector{T} = [ fill(0.0,B_dim)..., -1.0 ]

  DIIS_coeff, B, ipiv = LinearAlgebra.LAPACK.gesv!(B, DIIS_coeff)

  F_DIIS::Matrix{T} = zeros(size(F_array[1],1),size(F_array[1],2))
  for index::Int64 in 1:B_dim
    F_DIIS += DIIS_coeff[index]*F_array[index]
  end

  return F_DIIS
end

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
@inline function index(a::Int64,b::Int64)
  index::Int64 = (a*(a-1)) >> 1 #bitwise divide by 2
  index += b
  return index
end
