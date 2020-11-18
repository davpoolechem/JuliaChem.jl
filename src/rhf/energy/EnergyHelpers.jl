using Base.Threads

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
  #return ceil(Int,(-1.0+√(1+8*input))/2.0)
  return Base.fptosi(Int, Base.ceil_llvm((-1.0 + 
    Base.Math.sqrt_llvm(float(1+8*input)))/2.0))
    #return ccall((:decompose, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libjeri.so"),
  #  Int64, (Int64,), input)
end

function compute_enuc(mol::Molecule)
  E_nuc = 0.0
  for iatom in 1:length(mol), jatom in 1:(iatom-1)
    ix = mol[iatom].atom_center[1] 
    jx = mol[jatom].atom_center[1] 

    iy = mol[iatom].atom_center[2] 
    jy = mol[jatom].atom_center[2] 

    iz = mol[iatom].atom_center[3]
    jz = mol[jatom].atom_center[3]
  
    distance = √((jx-ix)^2 + (jy-iy)^2 + (jz-iz)^2) 
    
    E_nuc += mol[iatom].atom_id*mol[jatom].atom_id/distance
  end 
  
  return E_nuc
end
 
function compute_overlap(S::Matrix{Float64}, basis::Basis,
  jeri_oei_engine)

  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
    
    apos = basis[ash].pos
    bpos = basis[bsh].pos
       
    S_block_JERI = zeros(Float64,(abas*bbas,))
    JERI.compute_overlap_block(jeri_oei_engine, S_block_JERI, ash, bsh, 
      length(S_block_JERI))
    axial_normalization_factor(S_block_JERI, basis[ash], basis[bsh])

    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      S[max(iorb,jorb),min(iorb,jorb)] = S_block_JERI[idx]
      
      idx += 1 
    end
  end
 
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      S[min(iorb,jorb),max(iorb,jorb)] = S[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

function compute_ke(T::Matrix{Float64}, basis::Basis, 
  jeri_oei_engine)

  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
    
    apos = basis[ash].pos
    bpos = basis[bsh].pos
       
    T_block_JERI = zeros(Float64,(abas*bbas,))
    JERI.compute_kinetic_block(jeri_oei_engine, T_block_JERI, ash, bsh, 
      length(T_block_JERI))
    axial_normalization_factor(T_block_JERI, basis[ash], 
      basis[bsh])

    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      T[max(iorb,jorb),min(iorb,jorb)] = T_block_JERI[idx]
      
      idx += 1 
    end
  end
  
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      T[min(iorb,jorb),max(iorb,jorb)] = T[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

function compute_nah(V::Matrix{Float64}, mol::Molecule, 
  basis::Basis, jeri_oei_engine)
  
  #== define ncenter ==#
  #=
  ncenter::Int64 = length(mol)
  
  Z = Vector{Float64}([])
  x = Vector{Float64}([])
  y = Vector{Float64}([])
  z = Vector{Float64}([])

  for atom in mol 
    push!(Z, convert(Float64,atom.atom_id))  
    push!(x, atom.atom_center[1])  
    push!(y, atom.atom_center[2])  
    push!(z, atom.atom_center[3])  
  end
  =#
  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
    
    apos = basis[ash].pos
    bpos = basis[bsh].pos
       
    V_block_JERI = zeros(Float64,(abas*bbas,))
    JERI.compute_nuc_attr_block(jeri_oei_engine, V_block_JERI, ash, bsh, 
      length(V_block_JERI))
    axial_normalization_factor(V_block_JERI, basis[ash], 
      basis[bsh])
  
    idx = 1
    for ibas in 0:abas-1, jbas in 0:bbas-1
      iorb = apos + ibas
      jorb = bpos + jbas
      
      V[max(iorb,jorb),min(iorb,jorb)] = V_block_JERI[idx]
      
      idx += 1 
    end
  end
  
  for iorb in 1:basis.norb, jorb in 1:iorb
    if iorb != jorb
      V[min(iorb,jorb),max(iorb,jorb)] = V[max(iorb,jorb),min(iorb,jorb)]
    end
  end
end

function compute_schwarz_bounds(schwarz_bounds::Matrix{Float64}, 
  basis::Basis, nsh::Int64)

  max_am = max_ang_mom(basis) 
  eri_quartet_batch = Vector{Float64}(undef,eri_quartet_batch_size(max_am))
  jeri_schwarz_engine = JERI.TEIEngine(basis.basis_cxx, basis.shpdata_cxx)

  for ash in 1:nsh, bsh in 1:ash
    fill!(eri_quartet_batch, 0.0)
    
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
    abshp = triangular_index(ash, bsh)
 
    JERI.compute_eri_block(jeri_schwarz_engine, eri_quartet_batch, 
      ash, bsh, ash, bsh, abshp, abshp, abas*bbas, abas*bbas)
 
    schwarz_bounds[ash, bsh] = sqrt(maximum(abs.(eri_quartet_batch)) )
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

function DIIS(F::Matrix{Float64}, e_array::Vector{Matrix{Float64}}, 
  F_array::Vector{Matrix{Float64}}, B_dim::Int64)
  
  B = Matrix{Float64}(undef,B_dim+1,B_dim+1)
  for i in 1:B_dim, j in 1:B_dim
    B[i,j] = LinearAlgebra.BLAS.dot(length(e_array[i]), e_array[i], 1, 
      e_array[j], 1)

	  B[i,B_dim+1] = -1
	  B[B_dim+1,i] = -1
	  B[B_dim+1,B_dim+1] =  0
  end
  #DIIS_coeff::Vector{Float64} = [ fill(0.0,B_dim)..., -1.0 ]
  DIIS_coeff::Vector{Float64} = vcat(zeros(B_dim), [-1.0])

  #DIIS_coeff[:], B[:,:], ipiv = LinearAlgebra.LAPACK.gesv!(B, DIIS_coeff)
  DIIS_coeff[:], B[:,:], ipiv = LinearAlgebra.LAPACK.sysv!('U', B, DIIS_coeff)
  
  #fill!(F, zero(Float64))
  LinearAlgebra.BLAS.scal!(length(F), 0.0, F, 1) 
  for index in 1:B_dim
    F .+= DIIS_coeff[index] .* F_array[index]
  end
end

function axial_normalization_factor(oei, ash, bsh)
  ama = ash.am
  amb = bsh.am

  na = ash.nbas
  nb = bsh.nbas

  ab = 0 
  for asize::Int64 in 0:(na-1), bsize::Int64 in 0:(nb-1)
    ab += 1 
   
    anorm = axial_norm_fact[asize+1,ama]
    bnorm = axial_norm_fact[bsize+1,amb]
    
    abnorm = anorm*bnorm 
    oei[ab] *= abnorm
  end
end

function eri_quartet_batch_size(max_am)
  return am_to_nbas_cart(max_am)^4
end
