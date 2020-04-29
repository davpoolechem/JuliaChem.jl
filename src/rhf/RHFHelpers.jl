using JCModules.BasisStructs
using JCModules.MolStructs
using JCModules.Globals

using Base.Threads
using MATH
using JLD

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

function compute_schwarz_bounds(schwarz_bounds::Matrix{Float64}, nsh::Int64)
  eri_quartet_batch = Vector{Float64}(undef,81)
  simint_workspace = Vector{Float64}(undef,10000)

  for ash in 1:nsh, bsh in 1:ash
    SIMINT.compute_eris(ash, bsh, ash, bsh, eri_quartet_batch, 
      simint_workspace)
    
    schwarz_bounds[ash, bsh] = sqrt(maximum(eri_quartet_batch) )
  end

  for ash in 1:nsh, bsh in 1:ash
    if ash != bsh
      schwarz_bounds[min(ash,bsh),max(ash,bsh)] = 
        schwarz_bounds[max(ash,bsh),min(ash,bsh)]
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

>>>>>>> development
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
