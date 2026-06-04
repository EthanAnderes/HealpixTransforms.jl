module HealpixTransforms

using FFTW
using LinearAlgebra
using PyCall
using PyPlot
using XFields
import XFields: plan, size_in, size_out, eltype_in, eltype_out, Xmap, Xfourier
import LinearAlgebra: \, *

const module_dir  = joinpath(@__DIR__, "..") |> normpath
const UNSEEN = -1.6375e30 

F64 = Float64
C64 = Complex{Float64}

# Import healpy from PyCall
# ========================


# Transforms ℍ0 and ℍ02 # Can you make spin a type parameter?
# =========================================
export ℍ0, ℍ2, ℍ02

struct ℍ0{T<:Real} <: Transform
    nside::Int
    lmax::Int
    iter::Int 
    function ℍ0{T}(nside; lmax::Int=3*nside-1, iter::Int=3) where {T}
    	new{T}(nside, lmax, iter)
    end
    function ℍ0(nside; lmax::Int=3*nside-1, iter::Int=3)
        new{F64}(nside, lmax, iter)
    end

end 

struct ℍ2{T<:Real} <: Transform
    nside::Int
    lmax::Int
    iter::Int 
    function ℍ2{T}(nside; lmax=3*nside-1, iter::Int=3) where {T}
        new{T}(nside, lmax, iter)
    end
    function ℍ2(nside; lmax=3*nside-1, iter::Int=3)
        new{F64}(nside, lmax, iter)
    end
end 

struct ℍ02{T<:Real} <: Transform
    nside::Int
    lmax::Int
    iter::Int 
	function ℍ02{T}(nside; lmax=3*nside-1, iter::Int=3) where {T}
		new{T}(nside, lmax, iter)
    end
    function ℍ02(nside; lmax=3*nside-1, iter::Int=3)
        new{F64}(nside, lmax, iter)
    end
end 

Unionℍ{T} = Union{ℍ0{T}, ℍ2{T}, ℍ02{T}}

@inline size_in(h::ℍ0)   = (n_pix(h),)
@inline size_in(h::ℍ2)   = (n_pix(h),2) # QU
@inline size_in(h::ℍ02)  = (n_pix(h),3) # TQU

@inline size_out(h::ℍ0)   = (n_lm(h),)  
@inline size_out(h::ℍ2)   = (n_lm(h),2) # EB
@inline size_out(h::ℍ02)  = (n_lm(h),3) # TEB

@inline eltype_in(h::ℍ0{T})  where {T} = T
@inline eltype_in(h::ℍ2{T})  where {T} = T
@inline eltype_in(h::ℍ02{T}) where {T} = T

@inline eltype_out(h::ℍ0{T})  where {T} = Complex{T}
@inline eltype_out(h::ℍ2{T})  where {T} = Complex{T}
@inline eltype_out(h::ℍ02{T}) where {T} = Complex{T}

# There is nothing to pre-process for healpix transform.
# So we simply set plan to the idenitity and then specify
# how plan operates on the Array storage
@inline plan(h::Unionℍ) = h

# ℍ0
function *(h::ℍ0{T}, tx::Array{T,1}) where {T<:Real} 
    pyimport("healpy").map2alm(tx, lmax=h.lmax, iter=h.iter, pol=false)::Array{Complex{T},1}
end

function \(h::ℍ0{T}, tlm::Array{Complex{T},1}) where {T<:Real} 
    pyimport("healpy").alm2map(tlm, h.nside, lmax=h.lmax, pol=false)::Array{T,1}
end

# ℍ2

function *(h::ℍ2{T}, qux::Array{T,2})::Array{Complex{T},2} where {T<:Real} 
    elm, blm = pyimport("healpy").map2alm_spin((qux[:,1], qux[:,2]), 2, lmax=h.lmax)
    hcat(elm, blm)
end

function \(h::ℍ2{T}, eblm::Array{Complex{T},2})::Array{T,2} where {T<:Real} 
    mmax = h.lmax
    qx, ux = pyimport("healpy").sphtfunc.alm2map_spin((eblm[:,1], eblm[:,2]), h.nside, 2, h.lmax, mmax)
    hcat(qx, ux)
end

# ℍ02

function *(h::ℍ02{T}, tqux::Array{T,2}) where {T<:Real} 
	teblm = pyimport("healpy").map2alm(tqux', lmax=h.lmax, iter=h.iter, pol=true)
    Array(transpose(teblm))::Array{Complex{T},2}
end

function \(h::ℍ02{T}, teblm::Array{Complex{T},2}) where {T<:Real} 
	tqux  = pyimport("healpy").alm2map(transpose(teblm), h.nside, lmax=h.lmax, pol=true)
    Array(transpose(tqux))::Array{T,2}
end

# conveniet methods

spin0(h::Unionℍ{T})  where {T}  = ℍ0{T}(h.nside, lmax=h.lmax, iter=h.iter)
spin2(h::Unionℍ{T})  where {T}  = ℍ2{T}(h.nside, lmax=h.lmax, iter=h.iter)
spin02(h::Unionℍ{T}) where {T} = ℍ02{T}(h.nside, lmax=h.lmax, iter=h.iter)

# Extras
# =====================================
include("grid.jl")
include("pixels.jl")
include("plot.jl")
include("extras.jl")

end
