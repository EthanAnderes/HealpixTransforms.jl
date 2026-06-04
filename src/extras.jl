

# white noie
# ==========

function white_noise(tmH::Unionℍ)
	sqΩ  = .√(Ωpix(tmH))
	return Xmap(tmH, randn(eltype_in(tmH), size_in(tmH)) ./ sqΩ) 
end

# easy contruction of DiagOp's in ℓ,m space
# =========================================
"""
```
spectra2DiagOp(tmH::Unionℍ, fℓ::Function)
spectra2DiagOp(tmH::ℍ0,  fℓ)
spectra2DiagOp(tmH::ℍ2,  eℓ, bℓ)
spectra2DiagOp(tmH::ℍ02, tℓ, eℓ, bℓ)
```
all return `C, C_ℓ`. `C` is the field operator and `C_ℓ` has columns, effectivly the hcat of 
input spectra, but is isolated to `0:tmH.lmax`. These `C_ℓ` are generally used for plotting.
"""
function spectra2DiagOp end

function spectra2DiagOp(tmH::Unionℍ, fℓ::Function)
	lmax = tmH.lmax
	ℓrng = 0:lmax
	ls,  = lm(lmax)
	fls = fℓ.(ls) .+ Xfourier(tmH)[!]
	DiagOp(Xfourier(tmH, fls)), fℓ.(ℓrng)
end

function spectra2DiagOp(tmH::ℍ0, fℓ::AbstractVector)
	lmax = tmH.lmax
	ℓrng = 0:lmax
	@assert length(fℓ) >= lmax+1
	ls,  = lm(lmax)
	fls = fℓ[ls .+ 1]
	DiagOp(Xfourier(tmH, fls)), fℓ[ℓrng .+ 1]
end

function spectra2DiagOp(tmH::ℍ2, eℓ::AbstractVector, bℓ::AbstractVector)
	lmax = tmH.lmax
	ℓrng = 0:lmax
	@assert length(eℓ) >= lmax+1
	@assert length(bℓ) >= lmax+1
	ls,  = lm(lmax)
	els = eℓ[ls .+ 1] 
	bls = bℓ[ls .+ 1] 
	DiagOp(Xfourier(tmH, hcat(els,bls))), hcat(eℓ[ℓrng .+ 1], bℓ[ℓrng .+ 1])
end

function spectra2DiagOp(tmH::ℍ02, tℓ::AbstractVector, eℓ::AbstractVector, bℓ::AbstractVector)
	lmax = tmH.lmax
	ℓrng = 0:lmax
	@assert length(tℓ) >= lmax+1
	@assert length(eℓ) >= lmax+1
	@assert length(bℓ) >= lmax+1
	ls,  = lm(lmax)
	tls = tℓ[ls .+ 1] 
	els = eℓ[ls .+ 1] 
	bls = bℓ[ls .+ 1] 
	DiagOp(Xfourier(tmH, hcat(tls,els,bls))), hcat(tℓ[ℓrng .+ 1], eℓ[ℓrng .+ 1], bℓ[ℓrng .+ 1])
end

# Anafast wrapper
# ===============
"""
anafast returns NamedTuple{:TT} for input <: Xfield{ℍ0}
anafast returns NamedTuple{:EE, :BB, :EB} for input <: Xfield{ℍ2}
anafast returns NamedTuple{:TT, :EE, :BB, :TE, :EB, :TB} for input <: Xfield{ℍ02}
"""
function anafast end

# spin 0
function anafast(f::Xfield{<:ℍ0}; lmax=fieldtransform(f).lmax)
	(; TT=pyimport("healpy").anafast(f[:]; pol=false, lmax=lmax))
end
function anafast(f::Xfield{<:ℍ0}, g::Xfield{<:ℍ0}; lmax=fieldtransform(f).lmax)
	(; TT=pyimport("healpy").anafast((f[:], g[:]); pol=false, lmax=lmax))
end

# spin 2
function anafast(f::Xfield{<:ℍ2}; lmax=fieldtransform(f).lmax, Umult=1)
	# set Umult to -1 if you want to convert IAU <-> Cosmo pol convention
	qu   = Array.(eachcol(f[:]))
	tqu  = (zero(qu[1]), qu[1], Umult * qu[2]) 
	ana_out = pyimport("healpy").anafast(tqu; pol=true, lmax=lmax)
	(;EE=ana_out[2,:], BB=ana_out[3,:], EB=ana_out[5,:])
end
function anafast(f::Xfield{<:ℍ2}, g::Xfield{<:ℍ2}; lmax=fieldtransform(f).lmax, Umult=1)
	qu_f  = Array.(eachcol(f[:])) 
	qu_g  = Array.(eachcol(g[:])) 
	tqu_f = (zero(qu_f[1]), qu_f[1], Umult * qu_f[2]) 
	tqu_g = (zero(qu_g[1]), qu_g[1], Umult * qu_g[2]) 
	ana_out = pyimport("healpy").anafast(tqu_f, tqu_g; pol=true, lmax=lmax)
	(;EE=ana_out[2,:], BB=ana_out[3,:], EB=ana_out[5,:])
end

# spin 02
function anafast(f::Xfield{<:ℍ02}; lmax=fieldtransform(f).lmax, Umult=1)
	tqu  = Array.(eachcol(f[:]))
	tqu[3] *= Umult
	ana_out = pyimport("healpy").anafast(tqu; pol=true, lmax=lmax)
	return (;
		TT=ana_out[1,:], 
		EE=ana_out[2,:], 
		BB=ana_out[3,:], 
		TE=ana_out[4,:], 
		EB=ana_out[5,:], 
		TB=ana_out[6,:]
	)
end
function anafast(f::Xfield{<:ℍ02}, g::Xfield{<:ℍ02}; lmax=fieldtransform(f).lmax, Umult=1)
	tqu_f = Array.(eachcol(f[:])) 
	tqu_g = Array.(eachcol(g[:]))
	tqu_f[3] *= Umult
	tqu_g[3] *= Umult
	ana_out = pyimport("healpy").anafast(tqu_f, tqu_g; pol=true, lmax=lmax) 
	return (;
		TT=ana_out[1,:], 
		EE=ana_out[2,:], 
		BB=ana_out[3,:], 
		TE=ana_out[4,:], 
		EB=ana_out[5,:], 
		TB=ana_out[6,:]
	)
end

# Gaussian beam construction
# ==========================

"""
gaussian_beam(l; fwhm_arcmin) acts as an a_lm multiplier 
"""
function gaussian_beam(l; fwhm_arcmin)
	fwhm_rad = deg2rad(fwhm_arcmin/60)
	σ² = fwhm_rad^2 / 8 / log(2)
	return exp(-σ²*l*(l+1)/2)
end 

# alm2map_der1 wrapper (needs updating)
# =====================================

function ∇(alm::Array{Complex{T},1}, h::ℍ0{T}) where {T<:Real}
    ot = pyimport("healpy").alm2map_der1(alm, h.nside, lmax=h.lmax)
    ax = ot[1,:]
    ∂θ_ax = ot[2,:]
    inv_sinθ_∂φ_ax = ot[3,:]
    return ∂θ_ax, inv_sinθ_∂φ_ax, ax
end

