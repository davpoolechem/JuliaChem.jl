module SIMINT

using JuliaChem.JCModules
using StaticArrays

@inline function initialize()
  ccall((:initialize_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid, ())
end
export initialize

@inline function finalize()
  ccall((:finalize_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid, ())
end
export finalize

@inline function reset()
  ccall((:reset_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid, ())
end
export reset

@inline function get_julia_shell_info(shell::Shell)
  ccall( (:get_julia_shell_info_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Ref{Shell},), Ref(shell) )
end
export get_julia_shell_info

@inline function get_simint_shell_info(shell::Int64)
  ccall( (:get_simint_shell_info_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64,), shell )
end
export get_simint_shell_info

@inline function allocate_shell_array(basis::Basis)
  nshell::Int64 = length(basis.shells)

  nshell_simint::Int64 = 0
  for shell in basis.shells
    nshell_simint += shell.nbas == 4 ? 2 : 1
  end

  ccall( (:allocate_shell_array_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64,Int64), nshell, nshell_simint )

  return nshell_simint
end
export allocate_shell_array

@inline function add_shell(shell::Shell)
  ccall( (:add_shell_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Ref{Shell},), Ref(shell) )
end
export add_shell

@inline function get_workmem(derivative_order::Int64,
  max_ang_momentum::Int64)
  workmem = ccall( (:get_workmem_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Int64,
    (Int64,Int64), derivative_order, max_ang_momentum )
  return workmem
end
export get_workmem

@inline function normalize_shells()
  ccall( (:normalize_shells_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    () )
end
export normalize_shells

@inline function precompute_shell_pair_data()
  ccall( (:precompute_shell_pair_data_c,
    "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"),
    Cvoid, () )
end
export normalize_shells

@inline function unnormalize_shell(shell::Shell)
  for iprim::Int64 in 1:shell.nprim
    ee::Float64 = 2*shell.exponents[iprim]
    facs::Float64 = (pi/ee)^1.5
    shell.coefficients[iprim] /= sqrt(facs)
  end
end
export unnormalize_shell

@inline function compute_overlap(ash, bsh, ovr::Vector{Float64})
  ccall( (:compute_overlap_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64, Int64, Ptr{Float64},), ash, bsh, ovr)
end
export compute_overlap

@inline function compute_ke(ash, bsh, ke::Vector{Float64})
  ccall( (:compute_ke_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64, Int64, Ptr{Float64},), ash, bsh, ke)
end
export compute_ke

@inline function compute_nah(ncenter, Z, x, y, z, ash, bsh, 
 nah::Vector{Float64})
 ccall( (:compute_nah_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
   Ptr{Float64}, Int64, Int64, Ptr{Float64},), 
   ncenter, Z, x, y, z, ash, bsh, nah)
end
export compute_nah

@inline function create_ij_shell_pair(ish::Int64, jsh::Int64)
  ccall( (:create_ij_shell_pair_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"),
    Cvoid, (Int64, Int64), ish, jsh)
end
export create_ij_shell_pair

@inline function allocate_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:allocate_kl_shell_pair_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export allocate_kl_shell_pair

@inline function create_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:create_kl_shell_pair_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export create_kl_shell_pair

@inline function fill_kl_shell_pair(ksh::Int64, lsh::Int64)
  ccall( (:fill_kl_shell_pair_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"),
    Cvoid, (Int64, Int64), ksh, lsh)
end
export fill_kl_shell_pair

@inline function compute_eris(ish, jsh, ksh, lsh, eri::Vector{Float64}, work::Vector{Float64})
  ccall( (:compute_eris_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64, Int64, Int64, Int64, Ptr{Float64},Ptr{Float64}), ish, jsh, ksh, lsh, eri, work)
end
export compute_eris

@inline function retrieve_eris(ish::Int64, jsh::Int64, ksh::Int64, lsh::Int64,
  eri::Vector{T}) where {T<:AbstractFloat}
  ccall( (:retrieve_eris_c, "/export/home/david/projects/Julia/JuliaChem.jl/deps/libjeri.so"), Cvoid,
    (Int64, Int64, Int64, Int64, Ptr{T}),
    ish, jsh, ksh, lsh, eri)
end
export retrieve_eris

end
