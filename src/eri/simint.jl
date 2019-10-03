module SIMINT

using JCModules
using StaticArrays

function initialize()
  ccall((:initialize_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
#  ccall((:initialize_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
end
export initialize

function finalize()
  ccall((:finalize_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
#  ccall((:finalize_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
end
export finalize

function reset()
  ccall((:reset_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
#  ccall((:reset_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid, ())
end
export reset

function get_julia_shell_info(shell::JCModules.BasisStructs.Shell)
  ccall( (:get_julia_shell_info_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:get_julia_shell_info_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Ref{JCModules.BasisStructs.Shell},), Ref(shell) )
end
export get_julia_shell_info

function get_simint_shell_info(shell::Int64)
  ccall( (:get_simint_shell_info_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:get_simint_shell_info_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Int64,), shell )
end
export get_simint_shell_info

function allocate_shell_array(basis::JCModules.BasisStructs.Basis)
  nshell::Int64 = length(basis.shells)

  nshell_simint::Int64 = 0
  for shell in basis.shells
    nshell_simint += shell.nbas == 4 ? 2 : 1
  end

  ccall( (:allocate_shell_array_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:allocate_shell_array_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Int64,Int64), nshell, nshell_simint )

  return nshell_simint
end
export allocate_shell_array

function add_shell(shell::JCModules.BasisStructs.Shell)
  ccall( (:add_shell_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:add_shell_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Ref{JCModules.BasisStructs.Shell},), Ref(shell) )
end
export add_shell

function normalize_shells()
  ccall( (:normalize_shells_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:normalize_shells_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    () )
end
export normalize_shells

function unnormalize_shell(shell::JCModules.BasisStructs.Shell)
  for iprim::Int64 in 1:shell.nprim
    ee::Float64 = 2*shell.exponents[iprim]
    facs::Float64 = (pi/ee)^1.5
    shell.coefficients[iprim] /= sqrt(facs)
  end
end
export unnormalize_shell

#= compute shell pairs =#
function create_ij_shell_pair(ish::Int64, jsh::Int64)
  ccall( (:create_ij_shell_pair_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"),
    Cvoid, (Int64, Int64), ish, jsh)
end
export create_ij_shell_pair

function allocate_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:allocate_kl_shell_pair_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export allocate_kl_shell_pair

function create_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:create_kl_shell_pair_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export create_kl_shell_pair

function fill_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:fill_kl_shell_pair_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export fill_kl_shell_pair

function compute_eris(eri::Vector{Float64})
#  ccall( (:compute_eris_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
  ccall( (:compute_eris_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Ptr{Float64},), eri)
end
export compute_eris

#=
function retrieve_eris(ish::Int64, jsh::Int64, ksh::Int64, lsh::Int64,
  eri::Vector{T}) where {T<:AbstractFloat}

  ccall( (:retrieve_eris_c, "/home/david/Server/Documents/Programming/julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
#  ccall( (:retrieve_eris_c, "/export/home/david/projects/Julia/JuliaChem.jl/src/eri/libsimint.so"), Cvoid,
    (Int64, Int64, Int64, Int64, Ptr{T}),
    ish, jsh, ksh, lsh, eri)
end
export retrieve_eris
=#

end
