module HealpixTransforms

using FFTW
using LinearAlgebra
using PyCall
using XFields
import XFields: plan, size_in, size_out, eltype_in, eltype_out, Xmap, Xfourier

const module_dir  = joinpath(@__DIR__, "..") |> normpath
const UNSEEN = -1.6375e30 

F64 = Float64
C64 = Complex{Float64}

# Transforms ℍ0 and ℍ02 # Can you make spin a type parameter?
# =========================================
export ℍ0, ℍ02

struct ℍ0 <: Transform{F64,1}
    nside::Int
    lmax::Int
    iter::Int 
    function ℍ0(nside; lmax::Int=3*nside-1, iter::Int=0)
    	new(nside, lmax, iter)
    end
end 

struct ℍ02 <: Transform{F64,2}
    nside::Int
    lmax::Int
    iter::Int 
	function ℍ02(nside; lmax=3*nside-1, iter=0)
		new(nside, lmax, iter)
    end
end 

Unionℍ = Union{ℍ0, ℍ02}

@inline size_in(h::ℍ0)  = (n_pix(h),)
@inline size_in(h::ℍ02) = (n_pix(h),3) # TQU

@inline size_out(h::ℍ0)  = (n_lm(h),)
@inline size_out(h::ℍ02) = (n_lm(h),3) #TEB

@inline eltype_in(h::ℍ0)  = F64
@inline eltype_in(h::ℍ02) = F64

@inline eltype_out(h::ℍ0)  = C64
@inline eltype_out(h::ℍ02) = C64

# There is nothing to pre-process for healpix transform.
# So we simply set plan to the idenitity and then specify
# how plan operates on the Array storage
@inline plan(h::Unionℍ) = h

function Base.:*(h::ℍ0, tx::Array{F64,1}) 
    hp  = pyimport("healpy") 
    hp.map2alm(tx, lmax=h.lmax, iter=h.iter, pol=false)::Array{C64,1}
end

function Base.:\(h::ℍ0, tlm::Array{C64,1})
    hp  = pyimport("healpy")
    hp.alm2map(tlm, h.nside, lmax=h.lmax, pol=false, verbose=false)::Array{F64,1}
end

function Base.:*(h::ℍ02, tqux::Array{F64,2})
    hp  = pyimport("healpy") 
	teblm = hp.map2alm(tqux', lmax=h.lmax, iter=h.iter, pol=true)
    Array(transpose(teblm))::Array{C64,2}
end

function Base.:\(h::ℍ02, teblm::Array{C64,2})
    hp  = pyimport("healpy")  
	tqux  = hp.alm2map(transpose(teblm), h.nside, lmax=h.lmax, pol=true, verbose=false)
    Array(transpose(tqux))::Array{F64,2}
end

# Extra 
# =====================================
include("grid.jl")

# export n_pix, n_rings, n_lm, Ωpix, pix, lm

spin0(h::ℍ02) = ℍ0(h.nside, lmax=h.lmax, iter=h.iter)

spin02(h::ℍ0) = ℍ02(h.nside, lmax=h.lmax, iter=h.iter)

function ∇(alm::Array{C64,1}, h::ℍ0)
    hp  = pyimport("healpy")  
    ot = hp.alm2map_der1(alm, h.nside, lmax=h.lmax)
    ax = ot[1,:]
    ∂θ_ax = ot[2,:]
    inv_sinθ_∂φ_ax = ot[3,:]
    return ∂θ_ax, inv_sinθ_∂φ_ax, ax
end



end
