using JuliaChem.JCModules

using Base.Threads
using LinearAlgebra

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

function compute_nuc_grad(mol::Molecule)
  nuc_grad = zeros(Float64,(length(mol.atoms),3))
  for iatom in 1:length(mol.atoms), jatom in 1:length(mol.atoms)
    if iatom != jatom
      ix = mol.atoms[iatom].atom_center[1]
      jx = mol.atoms[jatom].atom_center[1]

      iy = mol.atoms[iatom].atom_center[2]
      jy = mol.atoms[jatom].atom_center[2]

      iz = mol.atoms[iatom].atom_center[3]
      jz = mol.atoms[jatom].atom_center[3]

      iatm = mol.atoms[iatom].atom_id
      jatm = mol.atoms[jatom].atom_id

      distance = √((jx-ix)^2 + (jy-iy)^2 + (jz-iz)^2)

      nuc_grad[jatom,1] -= iatm*jatm*(jx-ix)/(distance^3)
      nuc_grad[jatom,2] -= iatm*jatm*(jy-iy)/(distance^3)
      nuc_grad[jatom,3] -= iatm*jatm*(jz-iz)/(distance^3)
    end
  end

  return nuc_grad
end
 
function compute_overlap_grad(mol::Molecule, 
  basis::Basis, W::Matrix{Float64}, jeri_oei_grad_engine)

  #== define initial variables ==#
  natoms = length(mol.atoms)
  WS_grad = zeros(Float64,(natoms,3))
  
  ncoord = length(WS_grad) 
  shell_set = 2

  #== generate S derivative matrices ==#
  S_deriv = Vector{Matrix{Float64}}([ zeros(Float64,(basis.norb, basis.norb)) for i in 1:ncoord ])
  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
 
    apos = basis[ash].pos
    bpos = basis[bsh].pos
  
    iatom = basis[ash].atom_id
    jatom = basis[bsh].atom_id
 
    S_block_JERI = zeros(Float64,(abas*bbas*ncoord)) 
    
    JERI.compute_overlap_grad_block(jeri_oei_grad_engine, S_block_JERI, ash, bsh, 
      abas*bbas)

    #axial_normalization_factor(S_block_JERI, basis[ash], basis[bsh])
    
    #println("Julia Shells: $ash, $bsh")
    #println("Julia Atoms: $iatom, $jatom")

    shlset_idx = 1
    for ishlset in 1:shell_set
      atom = ishlset == 1 ? iatom : jatom

      for icoord in 1:3
        idx = 1
        op = 3*(atom-1) + icoord 

        #println("$ishlset, $icoord => $op, $shlset_idx")
        
        for ibas in 0:abas-1, jbas in 0:bbas-1
          iorb = apos + ibas
          jorb = bpos + jbas

          S_deriv[op][max(iorb,jorb),min(iorb,jorb)] += S_block_JERI[abas*bbas*(op-1) + idx]
          idx += 1
        end
        shlset_idx += 1
      end
    end
  end

  for ideriv in S_deriv
    for iorb in 1:basis.norb, jorb in 1:(iorb-1)
      ideriv[min(iorb,jorb),max(iorb,jorb)] = ideriv[max(iorb,jorb),min(iorb,jorb)]
    end
    #display(ideriv); println()
  end
  
  #== contract with energy-weighted density ==#
  for iatom in 1:natoms
    for icoord in 1:3
      deriv_idx = 3*(iatom-1)  + icoord
      WS_grad[iatom,icoord] = -(LinearAlgebra.dot(W, S_deriv[deriv_idx]))
    end
  end
  return WS_grad
end

function compute_kinetic_grad(mol::Molecule, 
  basis::Basis, P::Matrix{Float64}, jeri_oei_grad_engine)

  #== define initial variables ==#
  natoms = length(mol.atoms)
  PT_grad = zeros(Float64,(natoms,3))
  
  ncoord = length(PT_grad) 
  shell_set = 2

  #== generate T derivative matrices ==#
  T_deriv = Vector{Matrix{Float64}}([ zeros(Float64,(basis.norb, basis.norb)) for i in 1:ncoord ])
  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
 
    apos = basis[ash].pos
    bpos = basis[bsh].pos
  
    iatom = basis[ash].atom_id
    jatom = basis[bsh].atom_id
 
    T_block_JERI = zeros(Float64,(abas*bbas*ncoord)) 
    
    JERI.compute_kinetic_grad_block(jeri_oei_grad_engine, T_block_JERI, ash, bsh, 
      abas*bbas)

    #axial_normalization_factor(S_block_JERI, basis[ash], basis[bsh])
    
    #println("Julia Shells: $ash, $bsh")
    #println("Julia Atoms: $iatom, $jatom")

    shlset_idx = 1
    for ishlset in 1:shell_set
      atom = ishlset == 1 ? iatom : jatom

      for icoord in 1:3
        idx = 1
        op = 3*(atom-1) + icoord 

        #println("$ishlset, $icoord => $op, $shlset_idx")
        
        for ibas in 0:abas-1, jbas in 0:bbas-1
          iorb = apos + ibas
          jorb = bpos + jbas

          T_deriv[op][max(iorb,jorb),min(iorb,jorb)] += T_block_JERI[abas*bbas*(op-1) + idx]
          idx += 1
        end
        shlset_idx += 1
      end
    end
  end

  for ideriv in T_deriv
    for iorb in 1:basis.norb, jorb in 1:(iorb-1)
      ideriv[min(iorb,jorb),max(iorb,jorb)] = ideriv[max(iorb,jorb),min(iorb,jorb)]
    end
    #display(ideriv); println()
  end
  
  #== contract with energy-weighted density ==#
  for iatom in 1:natoms
    for icoord in 1:3
      deriv_idx = 3*(iatom-1)  + icoord
      PT_grad[iatom,icoord] = LinearAlgebra.dot(P, T_deriv[deriv_idx])
    end
  end
  return PT_grad
end

function compute_nuc_attr_grad(mol::Molecule, 
  basis::Basis, P::Matrix{Float64}, jeri_oei_grad_engine)

  #== define initial variables ==#
  natoms = length(mol.atoms)
  PV_grad = zeros(Float64,(natoms,3))
  
  ncoord = length(PV_grad) 
  shell_set = 2 + length(mol.atoms)

  #== generate T derivative matrices ==#
  V_deriv = Vector{Matrix{Float64}}([ zeros(Float64,(basis.norb, basis.norb)) for i in 1:ncoord ])
  for ash in 1:length(basis), bsh in 1:ash
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
 
    apos = basis[ash].pos
    bpos = basis[bsh].pos
  
    iatom = basis[ash].atom_id
    jatom = basis[bsh].atom_id
 
    V_block_JERI = zeros(Float64,(abas*bbas*ncoord)) 
    
    JERI.compute_nuc_attr_grad_block(jeri_oei_grad_engine, V_block_JERI, ash, bsh, 
      abas*bbas)

    #axial_normalization_factor(S_block_JERI, basis[ash], basis[bsh])
    
    #println("Julia Shells: $ash, $bsh")
    #println("Julia Atoms: $iatom, $jatom")

    shlset_idx = 1
    for ishlset in 1:shell_set
      atom = ishlset%2 != 0 ? iatom : jatom

      for icoord in 1:3
        idx = 1
        op = 3*(atom-1) + icoord 

        #println("$ishlset, $icoord => $op, $shlset_idx")
        
        for ibas in 0:abas-1, jbas in 0:bbas-1
          iorb = apos + ibas
          jorb = bpos + jbas

          V_deriv[op][max(iorb,jorb),min(iorb,jorb)] += V_block_JERI[abas*bbas*(op-1) + idx]
          idx += 1
        end
        shlset_idx += 1
      end
    end
  end

  for ideriv in V_deriv
    for iorb in 1:basis.norb, jorb in 1:(iorb-1)
      ideriv[min(iorb,jorb),max(iorb,jorb)] = ideriv[max(iorb,jorb),min(iorb,jorb)]
    end
    #display(ideriv); println()
  end
  
  #== contract with energy-weighted density ==#
  for iatom in 1:natoms
    for icoord in 1:3
      deriv_idx = 3*(iatom-1)  + icoord
      PV_grad[iatom,icoord] = LinearAlgebra.dot(P, V_deriv[deriv_idx])
    end
  end
  return PV_grad
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
    B[i,j] = LinearAlgebra.dot(e_array[i], e_array[j])

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
