using JuliaChem.JCModules

using Base.Threads
using LinearAlgebra
using Printf

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

function axial_normalization_factor(oei, ash, bsh, ncoord)
  ama = ash.am
  amb = bsh.am

  na = ash.nbas
  nb = bsh.nbas

  for icoord in 1:ncoord
    ab = (icoord-1)*na*nb 
    for asize::Int64 in 0:(na-1), bsize::Int64 in 0:(nb-1)
      ab += 1 
   
      anorm = axial_norm_fact[asize+1,ama]
      bnorm = axial_norm_fact[bsize+1,amb]
    
      abnorm = anorm*bnorm 
      oei[ab] *= abnorm
    end
  end
end

function compute_dipole(mol::Molecule, 
  basis::Basis, P::Matrix{Float64}, jeri_prop_engine)

  #== define initial variables ==#
  natoms = length(mol.atoms)
  dipole = zeros(Float64,(3,))
  
  ncoord = length(dipole) 

  #== compute nuclear contribution to dipole ==#
  nuc_dipole = similar(dipole)
  for atom in mol  
    nuc_dipole .+= atom.atom_id .* atom.atom_center
  end

  #== compute electronic contribution to dipole ==# 
  elec_dipole = similar(dipole)
  
  elec_dipole_matrix = Vector{Matrix{Float64}}([ zeros(Float64,(basis.norb, basis.norb)) for i in 1:ncoord ])
  for ash in 1:length(basis), bsh in 1:length(basis)
    abas = basis[ash].nbas
    bbas = basis[bsh].nbas
 
    apos = basis[ash].pos
    bpos = basis[bsh].pos
  
    iatom = basis[ash].atom_id
    jatom = basis[bsh].atom_id
 
    elec_dipole_block_JERI = zeros(Float64,(abas*bbas*ncoord)) 
    
    JERI.compute_dipole_block(jeri_prop_engine, elec_dipole_block_JERI, ash, bsh, 
      abas*bbas)

    axial_normalization_factor(elec_dipole_block_JERI, basis[ash], basis[bsh], ncoord)
    
    #println("Julia Shells: $ash, $bsh")
    #println("Julia Atoms: $iatom, $jatom")

    for icoord in 1:ncoord
      idx = 1

      #println("$ishlset, $icoord => $op, $shlset_idx")
        
      for ibas in 0:abas-1, jbas in 0:bbas-1
        iorb = apos + ibas
        jorb = bpos + jbas

        #elec_dipole_matrix[icoord][iorb,jorb] += elec_dipole_block_JERI[abas*bbas*(icoord-1) + idx]
        elec_dipole[icoord] += P[iorb,jorb] * elec_dipole_block_JERI[abas*bbas*(icoord-1) + idx]
        idx += 1
      end
    end
  end
  
  @printf("          %.8f   %.8f    %.8f         \n",
    nuc_dipole[1], nuc_dipole[2], nuc_dipole[3])
  @printf("          %.8f   %.8f    %.8f         \n",
    elec_dipole[1], elec_dipole[2], elec_dipole[3])

  dipole .= 2.54174623 .* (nuc_dipole .- elec_dipole)  
  return dipole
end


function eri_quartet_batch_size(max_am)
  return am_to_nbas_cart(max_am)^4
end
