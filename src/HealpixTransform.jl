module HealpixTransform

using FFTW
using LinearAlgebra
using PyCall
#using PyPlot
#using Interpolations
using XFields: Transform
import XFields: plan, size_in, size_out, eltype_in, eltype_out

const module_dir  = joinpath(@__DIR__, "..") |> normpath
const hp  = pyimport("healpy") 
const hps = pyimport("healpy.sphtfunc") 

const UNSEEN = -1.6375e30 

RN  = Union{Float32,Float64}

# Transforms ℍ0{<:RN} and ℍ02{<:RN}
# =========================================
export ℍ0, ℍ02, ℍ0plan, ℍ02plan

struct ℍ0{Tf<:RN} <: Transform{Tf,1}
    nside::Int
    lmax::Int
    iter::Int 
    function ℍ0{Tf}(nside, lmax=3*nside-1, iter=0) where Tf
    	new{Tf}(nside, lmax, iter)
    end
end 

struct ℍ02{Tf<:RN} <: Transform{Tf,2}
    nside::Int
    lmax::Int
	function ℍ02{Tf}(nside, lmax=3*nside-1, iter=0) where Tf
		new{Tf}(nside, lmax, iter)
    end
end 

Unionℍ{Tf} = Union{ℍ0{Tf}, ℍ02{Tf}}

@inline size_in(h::ℍ0)  = (n_pix(nside),)
@inline size_in(h::ℍ02) = (n_pix(nside),3) # TQU

@inline size_out(h::ℍ0)  = (n_lm(h.lmax),)
@inline size_out(h::ℍ02) = (n_lm(h.lmax),3) #TEB

@inline eltype_in(h::ℍ0{Tf})  where {Tf}  = Tf
@inline eltype_in(h::ℍ02{Tf}) where {Tf}  = Tf

@inline eltype_out(h::ℍ0{Tf})  where {Tf} = Complex{Tf}
@inline eltype_out(h::ℍ02{Tf}) where {Tf} = Complex{Tf}

# There is nothing to pre-process for healpix transform.
# So we simply set plan to the idenitity and then specify
# how plan operates on the Array storage
@inline plan(h::Unionℍ{Tf}) where {Tf<:RN} = h

function Base.:*(h::ℍ0{Tf}, tx::Array{Tf,1}) where Tf 
	hps.map2alm(tx, lmax=h.lmax, iter=h.iter, pol=false)::Array{Complex{Tf},1} 
end

function Base.:\(h::ℍ0{Tf}, tlm::Array{Complex{Tf},1}) where Tf 
	hps.alm2map(tlm, h.nside, lmax=h.lmax, pol=false, verbose=false)::Array{Tf,1}
end

function Base.:*(h::ℍ02{Tf}, tqux::Array{Tf,2}) where Tf 
    #ebℓm_vec = hps._sphtools.map2alm_spin_healpy(tqux, 2, lmax=lmax) 
    tup_tqux  = (tqux[:,1], tqux[:,2], tqux[:,3])
	tup_teblm = hps.map2alm(tup_tqux, lmax=h.lmax, iter=h.iter, pol=true)
    return hcat(tup_teblm[1], tup_teblm[2], tup_teblm[3])::Array{Complex{Tf},2}
end

function Base.:\(h::ℍ02{Tf}, teblm::Array{Complex{Tf},2}) where Tf 
    #qux_vec = hps._sphtools.alm2map_spin_healpy(ebℓm, nside, 2, lmax)
    tup_telmb = (telmb[:,1], telmb[:,2], telmb[:,3])
	tup_tqux = hps.alm2map(tup_telmb, h.nside, lmax=h.lmax, pol=true, verbose=false)
    return hcat(tup_tqux[1], tup_tqux[2], tup_tqux[3])::Array{Tf,2}
end


# Extra grid information
# =====================================
export n_pix, n_rings, n_lm, Ωpix, pix, lm

include("grid.jl")


end
