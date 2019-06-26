#=
function twoei_base(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
	F = deepcopy(H)
	J::Array{Float64,2} = zeros(7,7)
	K::Array{Float64,2} = zeros(7,7)

	for i::Int64 in 1:7
		for j::Int64 in 1:7
			for k::Int64 in 1:7
				for l::Int64 in 1:7
					ij::Int64 = index(i,j)
					kl::Int64 = index(k,l)
					ijkl::Int64 = index(ij,kl)
					#J_idx = 7*(μ-1)+ν
					#F[i,j] += D[k,l] * (2*tei[ijkl] - tei[ikjl])
					#println(ijkl,"        ",tei[ijkl])
					J[i,j] += D[k,l] * tei[ijkl]
					K[i,k] += D[j,l] * tei[ijkl]
				end
			end
		end
	end
	F += 2*J - K
	return F
end

function twoei_ij(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
	F = deepcopy(H)
	J::Array{Float64,2} = zeros(7,7)
	K::Array{Float64,2} = zeros(7,7)
	for i::Int64 in 1:7
		for j::Int64 in 1:i
			for k::Int64 in 1:7
				for l::Int64 in 1:7
					ij::Int64 = index(i,j)
					kl::Int64 = index(k,l)
					ijkl::Int64 = index(ij,kl)
					#J_idx = 7*(μ-1)+ν
					#F[i,j] += D[k,l] * (2*tei[ijkl] - tei[ikjl])
					#F[i,j] += 2 * D[k,l] * tei[ijkl]
					#F[i,k] += -1 * D[j,l] * tei[ijkl]
					val = (i == j) ? 0.5 : 1
					J[i,j] += val * D[k,l] * tei[ijkl]
					J[j,i] += val * D[k,l] * tei[ijkl]
					#J[k,l] += val * D[i,j] * tei[ijkl]
					#J[l,k] += val * D[i,j] * tei[ijkl]
					if (i == j)
						K[i,k] += D[j,l] * tei[ijkl]
					else
						K[i,k] += val * D[j,l] * tei[ijkl]
						K[j,k] += val * D[i,l] * tei[ijkl]
					end
				end
			end
			#F[j,i] = F[i,j]
		end
	end
	F += 2*J - K
	return F
end

function twoei_kl(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
	F = deepcopy(H)
	#F = zeros(7,7)
	J::Array{Float64,2} = zeros(7,7)
	K::Array{Float64,2} = zeros(7,7)
	for i::Int64 in 1:7
		for j::Int64 in 1:i
			for k::Int64 in 1:7
				for l::Int64 in 1:k
					ij::Int64 = index(i,j)
					kl::Int64 = index(k,l)
					ijkl::Int64 = index(ij,kl)
					#val = (k == l) ? 0.5 : 1
					#F[i,j] += val * D[k,l] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
					val = (i == j) ? 0.5 : 1
					val *= (k == l) ? 0.5 : 1
					#J[i,j] += val * D[k,l] * tei[ijkl]
					#J[j,i] += val * D[k,l] * tei[ijkl]
					#J[i,j] += val * D[l,k] * tei[ijkl]
					#J[j,i] += val * D[l,k] * tei[ijkl]
					J[k,l] += val * D[i,j] * tei[ijkl]
					J[k,l] += val * D[j,i] * tei[ijkl]
					J[l,k] += val * D[i,j] * tei[ijkl]
					J[l,k] += val * D[j,i] * tei[ijkl]
					K[i,k] += val * D[j,l] * tei[ijkl]
					K[i,l] += val * D[j,k] * tei[ijkl]
					K[j,k] += val * D[i,l] * tei[ijkl]
					K[j,l] += val * D[i,k] * tei[ijkl]
				end
			end
			#F[j,i] = F[i,j]
		end
	end
	#F += H
	F += 2*J - K
	return F
end

function twoei_kl_opt(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64},
	H::Array{Float64,2}, FLAGS::Flags)

	norb::Int64 = FLAGS.BASIS.NORB
	ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

	F = zeros(norb,norb)

	J::Array{Float64,2} = zeros(norb,norb)
	K::Array{Float64,2} = zeros(norb,norb)

	for i::Int64 in 1:norb
		for j::Int64 in 1:i
			ij::Int64 = index(i,j,ioff)

			for k::Int64 in 1:norb
				for l::Int64 in 1:k

					kl::Int64 = index(k,l,ioff)
					ijkl::Int64 = index(ij,kl,ioff)

					val::Float64 = (i == j) ? 0.5 : 1.0
					val::Float64 *= (k == l) ? 0.5 : 1.0

					J[k,l] += 2.0 * val * D[i,j] * tei[ijkl]
					J[l,k] += 2.0 * val * D[i,j] * tei[ijkl]

					K[i,k] += val * D[j,l] * tei[ijkl]
					K[i,l] += val * D[j,k] * tei[ijkl]
					K[j,k] += val * D[i,l] * tei[ijkl]
					K[j,l] += val * D[i,k] * tei[ijkl]
			#        if (i == j && k == l)
		#                K[i,k] += D[j,l] * tei[ijkl]
	#                elseif (i == j)
#                        K[i,k] += D[j,l] * tei[ijkl]
					 #    K[i,l] += D[j,k] * tei[ijkl]
					 #elseif (k == l)
			#            K[i,k] += D[j,l] * tei[ijkl]
			#            K[j,k] += D[i,l] * tei[ijkl]
			#        else
			#        end
				end
			end
			#F[j,i] = F[i,j]
		end
	end
	#F += H
	F = H + 2*J - K
	return F
end

function twoei_fock2a(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64}, H::Array{Float64,2})
	#F = deepcopy(H)
	F = zeros(7,7)
	J::Array{Float64,2} = zeros(7,7)
	Ki::Array{Float64,2} = zeros(7,7)
	Kj::Array{Float64,2} = zeros(7,7)
	for i::Int64 in 1:7
		for j::Int64 in 1:i
			ij::Int64 = index(i,j)
			for k::Int64 in 1:7
				for l::Int64 in 1:k
					kl::Int64 = index(k,l)
					ijkl::Int64 = index(ij,kl)
					val = (i == j) ? 0.5 : 1
					val *= (k == l) ? 0.5 : 1
					J[k,l] += 2 * val * D[i,j] * tei[ijkl]
					J[l,k] += 2 * val * D[i,j] * tei[ijkl]
					if (k == l)
						Ki[i,k] += 2 * val * D[j,l] * tei[ijkl]
						Kj[j,k] += 2 * val * D[i,l] * tei[ijkl]
					else
						Ki[i,k] += val * D[j,l] * tei[ijkl]
						Ki[i,l] += val * D[j,k] * tei[ijkl]
						Kj[j,k] += val * D[i,l] * tei[ijkl]
						Kj[j,l] += val * D[i,k] * tei[ijkl]
					end
			#        if (i == j && k == l)
		#                K[i,k] += D[j,l] * tei[ijkl]
	#                elseif (i == j)
#                        K[i,k] += D[j,l] * tei[ijkl]
					   #  K[i,l] += D[j,k] * tei[ijkl]
					 #elseif (k == l)
			#            K[i,k] += D[j,l] * tei[ijkl]
			#            K[j,k] += D[i,l] * tei[ijkl]
			#        else
			#        end
				end
			end
			#F[j,i] = F[i,j]
		end
	end
	#F += H
	F = H + 2*J - Ki - Kj
	return F
end

function twoei_edit(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
	F = deepcopy(H)
	#F = zeros(7,7)
	for i::Int64 in 1:7, j::Int64 in 1:i
		for k::Int64 in 1:7, l::Int64 in 1:k
			ij::Int64 = index(i,j)
			kl::Int64 = index(k,l)
			ijkl::Int64 = index(ij,kl)
			ik::Int64 = index(i,k)
			jl::Int64 = index(j,l)
			ikjl::Int64 = index(ik,jl)
			il::Int64 = index(i,l)
			jk::Int64 = index(j,k)
			iljk::Int64 = index(il,jk)
			eri = tei[ijkl]
			val = (k == l) ? 0.5 : 1
			val = (j == l) ? 0.5 : 1
			val3 = (j == k) ? 0.5 : 1
			F[i,j] += val * D[k,l] * (2*tei[ijkl] - tei[ikjl])
			F[i,j] += val * D[l,k] * (2*tei[ijkl] - tei[iljk])
			#F[i,j] += 4 * val * D[k,l] * (tei[ijkl])
			#F[i,k] -= val * D[j,l] * (tei[ijkl])
			#F[i,l] -= val * D[j,k] * (tei[ijkl])
			#F[i,j] += 2*D[k,l] * eri
			#F[i,k] -= D[j,l] * eri
			#F[k,i] = F[i,k]
			#F[l,i] = F[i,l]
		end
		F[j,i] = F[i,j]
		#F[i,j] += H[i,j]
	end
	return F
end

function twoei_six(F::Array{Float64,2}, D::Array{Float64,2},
	tei::Array{Float64,1}, H::Array{Float64,2}, FLAGS::Flags)

	ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:7*8))
	F = deepcopy(H)
	norb::Int64 = FLAGS.BASIS.NORB

	J::Array{Float64,2} = zeros(7,7)
	K::Array{Float64,2} = zeros(7,7)
	for i::Int64 in 1:7
		for j::Int64 in 1:i
			for k::Int64 in 1:j
				for l::Int64 in 1:k
					ij::Int64 = index(i,j,ioff)
					kl::Int64 = index(k,l,ioff)
					if (ij < kl) continue end
					ijkl::Int64 = index(ij,kl,ioff)
					ik::Int64 = index(i,k,ioff)
					jl::Int64 = index(j,l,ioff)
					ikjl::Int64 = index(ik,jl,ioff)
					il::Int64 = index(i,k,ioff)
					jk::Int64 = index(j,l,ioff)
					iljk::Int64 = index(il,jk,ioff)
					val::Float64 = (k == l) ? 0.5 : 1.0
					F[i,j] += val * D[k,l] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
					F[k,l] += val * D[i,j] * (4*tei[ijkl] - tei[ikjl] - tei[iljk])
					#val = (i == j) ? 0.5 : 1
					#val *= (k == l) ? 0.5 : 1
					#val *= (i == k && j == l) ? 0.5 : 1
					#J[i,j] += val * D[k,l] * tei[ijkl]
					#J[j,i] += val * D[k,l] * tei[ijkl]
					#J[i,j] += val * D[l,k] * tei[ijkl]
					#J[j,i] += val * D[l,k] * tei[ijkl]
					#J[k,l] += val * D[i,j] * tei[ijkl]
					#J[l,k] += val * D[i,j] * tei[ijkl]
					#J[k,l] += val * D[j,i] * tei[ijkl]
					#J[l,k] += val * D[j,i] * tei[ijkl]
					#K[i,k] += val * D[j,l] * tei[ijkl]
					#K[i,l] += val * D[j,k] * tei[ijkl]
					#K[j,k] += val * D[i,l] * tei[ijkl]
					#K[j,l] += val * D[i,k] * tei[ijkl]
				end
			end
			F[j,i] = F[i,j]
		end
	end
	F += H
	#F += 2*J - K
	return F
end

"""
	 twoei(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""
function twoei_distributed(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
	H::Array{Float64,2}, FLAGS::Flags, μν_idx::Int64)

	norb::Int64 = FLAGS.BASIS.NORB

	ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

	F = zeros(norb,norb)

	μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
	ν::Int64 = μν_idx%μ + 1
	μν::Int64 = index(μ,ν,ioff)

	for λσ_idx::Int64 in 1:ioff[norb]
		λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
		σ::Int64 = λσ_idx%λ + 1

		λσ::Int64 = index(λ,σ,ioff)
		μνλσ::Int64 = index(μν,λσ,ioff)

		val::Float64 = (μ == ν) ? 0.5 : 1.0
		val::Float64 *= (λ == σ) ? 0.5 : 1.0
		eri::Float64 = val * tei[μνλσ]

		if (eri <= 1E-10) continue end

		F[λ,σ] += 4.0 * D[μ,ν] * eri
		F[σ,λ] += 4.0 * D[μ,ν] * eri
		F[μ,λ] -= D[ν,σ] * eri
		F[μ,σ] -= D[ν,λ] * eri
		F[ν,λ] -= D[μ,σ] * eri
		F[ν,σ] -= D[μ,λ] * eri
	end

	return F
end

"""
	 twoei_threaded(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step, using multithreading.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""

function twoei_threaded(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
	H::Array{Float64,2}, FLAGS::Flags)

	norb::Int64 = FLAGS.BASIS.NORB
	ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))
	F = deepcopy(H)
	mutex = Base.Threads.Mutex()

	Threads.@threads for μν_idx::Int64 in 1:ioff[norb]
		μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
		ν::Int64 = μν_idx%μ + 1

		μν::Int64 = index(μ,ν,ioff)

		F_priv::Array{Float64,2} = zeros(norb,norb)
		for λσ_idx::Int64 in 1:ioff[norb]
			λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
			σ::Int64 = λσ_idx%λ + 1

			λσ::Int64 = index(λ,σ,ioff)
			μνλσ::Int64 = index(μν,λσ,ioff)

			val::Float64 = (μ == ν) ? 0.5 : 1.0
			val::Float64 *= (λ == σ) ? 0.5 : 1.0
			eri::Float64 = val * tei[μνλσ]

			if (eri <= 1E-10) continue end

			F_priv[λ,σ] += 4.0 * D[μ,ν] * eri
			F_priv[σ,λ] += 4.0 * D[μ,ν] * eri

			F_priv[μ,λ] -= D[ν,σ] * eri
			F_priv[μ,σ] -= D[ν,λ] * eri
			F_priv[ν,λ] -= D[μ,σ] * eri
			F_priv[ν,σ] -= D[μ,λ] * eri
		end

		lock(mutex)
		F += F_priv
		unlock(mutex)
	end

	return F
end

"""
	 twoei_threaded(F::Array{Float64}, D::Array{Float64}, tei::Array{Float64}, H::Array{Float64})
Summary
======
Perform Fock build step, using multithreading.

Arguments
======
F = Current iteration's Fock Matrix

D = Current iteration's Density Matrix

tei = Two-electron integral array

H = One-electron Hamiltonian Matrix
"""
function twoei_tasked(F::Array{Float64,2}, D::Array{Float64,2}, tei::Array{Float64,1},
	H::Array{Float64,2}, FLAGS::Flags)

	norb::Int64 = FLAGS.BASIS.NORB

	ioff::Array{Int64,1} = map((x) -> x*(x+1)/2, collect(1:norb*(norb+1)))

	F = deepcopy(H)
	J::Array{Float64,2} = zeros(norb,norb)
	K::Array{Float64,2} = zeros(norb,norb)
	#K::Array{Threads.Atomic{Float64},2} = fill(Threads.Atomic{Float64}(0.0),(norb,norb))

	#Threads.@threads for μν_idx::Int64 in 1:ioff[norb]
	@sync @async for μν_idx::Int64 in 1:ioff[norb]
		μ::Int64 = ceil(((-1+sqrt(1+8*μν_idx))/2))
		ν::Int64 = μν_idx%μ + 1

		μν::Int64 = index(μ,ν,ioff)

		for λσ_idx::Int64 in 1:ioff[norb]
			λ::Int64 = ceil(((-1+sqrt(1+8*λσ_idx))/2))
			σ::Int64 = λσ_idx%λ + 1

			λσ::Int64 = index(λ,σ,ioff)
			μνλσ::Int64 = index(μν,λσ,ioff)

			val::Float64 = (μ == ν) ? 0.5 : 1.0
			val::Float64 *= (λ == σ) ? 0.5 : 1.0
			eri::Float64 = val * tei[μνλσ]

			if (eri <= 1E-10) continue end

			J[λ,σ] += 2.0 * D[μ,ν] * eri
			J[σ,λ] += 2.0 * D[μ,ν] * eri

			K[μ,λ] += D[ν,σ] * eri
			K[μ,σ] += D[ν,λ] * eri
			K[ν,λ] += D[μ,σ] * eri
			K[ν,σ] += D[μ,λ] * eri
			#Threads.atomic_add!(K[μ,λ], D[ν,σ] * eri)
			#Threads.atomic_add!(K[μ,σ], D[ν,λ] * eri)
			#Threads.atomic_add!(K[ν,λ], D[μ,σ] * eri)
			#Threads.atomic_add!(K[ν,σ], D[μ,λ] * eri)
		end
	end

	F += 2*J - K
	#for i in 1:norb, j in 1:norb
	#    F[i,j] = 2*J[i,j] - K[i,j][]
	#end
	return F
end

=#
