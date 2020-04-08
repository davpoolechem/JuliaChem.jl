using JCModules.MolStructs
using MPI
using Base.Threads
using LinearAlgebra

function print_xyz(mol::MolStructs.Molecule)
  #== determine some pre-information ==#
  natoms = length(mol.atoms)
  comm = MPI.COMM_WORLD
  comm = MPI.COMM_WORLD

  #== print coordinates in xyz format ==#
  if (MPI.Comm_rank(comm) == 0)
    println("----------------------------------------          ")
    println("      Printing coordinates (bohr)                 ")
    println("            in xyz format...                      ")
    println("----------------------------------------          ")
    println(" ")
  end
  println(natoms) 
  println() 

  for iatom in mol.atoms
    if (MPI.Comm_rank(comm) == 0)
      println(iatom.symbol,"         ",iatom.atom_center[1],"     ",
        iatom.atom_center[2],"     ",iatom.atom_center[3])
    end 
  end
  
  println() 
end

function analyze_bond_lengths(mol::MolStructs.Molecule)
  #== determine some pre-information ==#
  natoms = length(mol.atoms)
  natoms_length = floor(Int64,natoms*(natoms+1)/2)
  comm = MPI.COMM_WORLD

  #== calculate bond lengths ==#
  bond_lengths = zeros(Float64,(natoms_length,))
  
  ijatom = 1
  for iatom in 1:natoms, jatom in 1:(iatom-1)
    diff = (mol.atoms[iatom].atom_center .- mol.atoms[jatom].atom_center).^2  
    bond_lengths[ijatom] = 0.52917724924*sqrt(reduce(+,diff))
    ijatom += 1
  end

  #== print bond lengths ==#
  if (MPI.Comm_rank(comm) == 0)
      println("----------------------------------------          ")
      println("        Printing bond lengths...                  ")
      println("----------------------------------------          ")
      println(" ")
      println("Atom #1   Atom #2     Bond length")
  end
  
  ijatom = 1
  for iatom in 1:natoms, jatom in 1:(iatom-1)
    if (MPI.Comm_rank(comm) == 0)
      println("   ",iatom,"         ",jatom,"     ",
        bond_lengths[ijatom])
    end
    ijatom += 1
  end

  if (MPI.Comm_rank(comm) == 0)
      println(" ")
  end

  return bond_lengths
end

function analyze_bond_angles(mol::MolStructs.Molecule, 
  bond_lengths::Vector{Float64})
  
  #== determine some pre-information ==#
  natoms = length(mol.atoms)
  natoms_length = floor(Int64,natoms*(natoms+1)*(natoms+2)/6)
  comm=MPI.COMM_WORLD

  #== calculate bond angles ==#
  bond_angles = zeros(Float64, (natoms_length,))

  ijatom = 1 
  ijkatom = 1
  for iatom in 1:natoms
    jkatom = 1
    for jatom in 1:(iatom-1), katom in 1:(jatom-1)
      if bond_lengths[ijatom] < 4.0 && bond_lengths[jkatom] < 4.0
        e_ji = Vector{Float64}(-(mol.atoms[jatom].atom_center .- 
          mol.atoms[iatom].atom_center) ./ bond_lengths[ijatom])

        e_jk = Vector{Float64}(-(mol.atoms[jatom].atom_center .- 
          mol.atoms[katom].atom_center) ./ bond_lengths[jkatom])
     
        e_ji ./= LinearAlgebra.norm(e_ji)
        e_jk ./= LinearAlgebra.norm(e_jk)

        bond_angles[ijkatom] = (360/(2π))*acos(LinearAlgebra.dot(e_ji,e_jk))
      
        ijkatom += 1 
      end
      jkatom += 1 
    end
    ijatom += 1
  end
  
  #= 
  for ijkatom in 1:natoms*natoms*natoms
      iatom::Int64 = ceil(ijkatom/(natoms*natoms))
      jatom::Int64 = ceil(ijkatom/natoms)%natoms+1
      katom::Int64 = (ijkatom%(natoms*natoms))%natoms + 1

      if (iatom > jatom && jatom > katom)
          if (bond_lengths[iatom,jatom] < 4.0 && bond_lengths[jatom,katom] < 4.0)
              e_ji = [ -(coord[jatom,1]-coord[iatom,1])/bond_lengths[iatom,jatom];
                       -(coord[jatom,2]-coord[iatom,2])/bond_lengths[iatom,jatom];
                       -(coord[jatom,3]-coord[iatom,3])/bond_lengths[iatom,jatom]; ]

              e_jk = [ -(coord[jatom,1]-coord[katom,1])/bond_lengths[jatom,katom];
                       -(coord[jatom,2]-coord[katom,2])/bond_lengths[jatom,katom];
                       -(coord[jatom,3]-coord[katom,3])/bond_lengths[jatom,katom]; ]

              e_ji = e_ji/LinearAlgebra.norm(e_ji)
              e_jk = e_jk/LinearAlgebra.norm(e_jk)

              bond_angles[iatom,jatom,katom] = (360/(2π))*acos(LinearAlgebra.dot(e_ji,e_jk))
          end
      end
  end
  =#

  #== print bond angles ==#
  if (MPI.Comm_rank(comm) == 0)
      println("----------------------------------------          ")
      println("         Printing bond angles...                  ")
      println("----------------------------------------          ")
      println(" ")
      println("Atom #1   Atom #2   Atom #3     Bond angle")
  end
  
  ijatom = 1 
  ijkatom = 1
  for iatom in 1:natoms 
    jkatom = 1
    for jatom in 1:(iatom-1), katom in 1:(jatom-1)
      if bond_lengths[ijatom] < 4.0 && bond_lengths[jkatom] < 4.0
        if (MPI.Comm_rank(comm) == 0)
          println("   ",iatom,"         ",jatom,"         ",katom,"     ",bond_angles[ijkatom])
        end
        ijkatom += 1
      end
      jkatom += 1 
    end
    ijatom += 1
  end
  if (MPI.Comm_rank(comm) == 0)
      println(" ")
  end

  return bond_lengths
end

#=
function test(ijkatom::Int64)
    natoms = 3

    iatom::Int64 = ceil(ijkatom/(natoms*natoms))
    jatom::Int64 = ceil(ijkatom/natoms)%natoms+1
    katom::Int64 = (ijkatom%(natoms*natoms))%natoms + 1
    #println(iatom,", ",jatom,", ",katom)
    return (iatom, jatom, katom)
end
"""
     coordinate_analysis(coord::Array{Float64,2})
Summary
======
Perform the core molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
"""
function coordinate_analysis(mol::MolStructs.Molecule)
  #== print coordinates ==#
  print_xyz(mol) 

  #== compute bond lengths ==#
  bond_lengths = analyze_bond_lengths(mol)

  #== compute bond angles ==#
  bond_angles = analyze_bond_angles(mol,bond_lengths)
end
=#
