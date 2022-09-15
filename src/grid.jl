
npix2nside(n_pix::Int) = Int(sqrt(n_pix/12))

n_pix(h::Unionℍ)  = n_pix(h.nside)
n_pix(nside::Int) = 12nside^2

n_rings(h::Unionℍ)  = n_rings(h.nside)
n_rings(nside::Int) = 4nside - 1

n_lm(h::Unionℍ) = n_lm(h.lmax)
n_lm(lmax::Int) = lmax * (lmax + 1) ÷ 2 + lmax + 1

Ωpix(h::Unionℍ{T}) where {T} = Ωpix(h.nside;T)
Ωpix(nside::Int;T=Float32) = T(4π / n_pix(nside)) # π / (3 nside^2)

function pix(h::Unionℍ{T}) where {T}
	θ, φ = pix(h.nside; T)
	return θ, φ
end
function pix(Nside::Int;T=Float32)
	θ_col, φ_col, idx_col, Δφ_col, nφ_col = θ_φ_idx_4_rings(Nside)
	n_col = length(θ_col)

	npx = n_pix(Nside)
	θ = zeros(T,npx)
	φ = zeros(T,npx)

	for  i_col in 1:length(idx_col)
		i_px      = idx_col[i_col]
		if i_col == length(idx_col)
			i_next_px = npx
		else  
			i_next_px = idx_col[i_col+1]
		end
		n_ring = length(i_px:i_next_px)
		
		θ[i_px:i_next_px] .= θ_col[i_col] 
		φ[i_px:i_next_px] .= φ_col[i_col] .+  Δφ_col[i_col] .* (0:n_ring-1)
	end

	return θ, φ
end
#=
using HealpixTransforms 
using HealpixTransforms: θ_φ_idx_4_rings, n_pix
using PyCall

function pix_test(nside::Int)
	hp  = pyimport("healpy") 
	θ, φ  = hp.pix2ang(nside, 0:(n_pix(nside)-1))
	return θ, φ
end

Nside = 2048

@time θ_orig, φ_orig = pix_test(Nside);
@time θ_new, φ_new   = pix(Nside);

θ_new ≈ θ_orig
φ_new ≈ φ_orig
=#


# l,m <-> index  
# -----------------------------------

index2lm(i, h::Unionℍ) = index2lm(i, h.lmax)

function index2lm(i, lmax)
	hp  = pyimport("healpy") 
	hp.Alm.getlm(lmax, i - 1)
end

lm2index(l, m, h::Unionℍ) = lm2index(l, m, h.lmax)

function lm2index(l, m, lmax)
	hp  = pyimport("healpy")
	hp.Alm.getidx(lmax, l, m) .+ 1
end

lm(h::Unionℍ) = lm(h.lmax)

function lm(lmax::Int) 
	nlm = n_lm(lmax)
	i = 1:nlm 
	m = @. ceil(Int, (2lmax + 1)/2 - √((2lmax + 1)^2 - 8*(i - 1 - lmax))/2 )
	l = @. i - 1 - m*(2lmax + 1 - m) ÷ 2
	return l, m 
end


"""
`alm2triangle(alms) -> alms_mat` arrange the healpix alms into matrix form.
Rows correspond to l (from 0:lmax).
Columns correspond to m (from 0:lmax). 
Example:
```
import HealpixTransforms as HT
Nside = 2048
lmax = 3Nside - 1
ls, ms = HT.lm(lmax) 
ls_mat = HT.alm2triangle(ls)
ms_mat = HT.alm2triangle(ms)
```
"""
function alm2triangle(alms::AbstractVector{T}) where T
	hp  = pyimport("healpy")
    lmax         = hp.sphtfunc.Alm.getlmax(length(alms))
    alms_mat     = fill(promote(T(0),NaN32)[2], lmax+1, lmax+1)
    start_ℓindex = end_ℓindex = 0
    for i = 1:lmax+1
        start_ℓindex = end_ℓindex + 1
        end_ℓindex   = start_ℓindex + (lmax+1) - i
        rng = start_ℓindex:end_ℓindex
        alms_mat[end-length(rng)+1:end,i] = alms[rng]
    end
    return Array(alms_mat)
end

# Note: to reduce mmax, just cut the end columns of the output of 
# alm2triangle.
function triangle2alm(tri_alms::AbstractMatrix{T}) where T
	hp   = pyimport("healpy")
	lmax = size(tri_alms,1)-1
	mmax = size(tri_alms,2)-1
	nlms = hp.sphtfunc.Alm.getsize(lmax, mmax)
    alms = fill(promote(T(0),NaN32)[2], nlms)
    indx_crr = 1
    for col = 1:mmax+1
    	x  = tri_alms[:,col]
    	fx = isfinite.(x)
    	alms[indx_crr:(indx_crr - 1 + sum(fx))] = x[fx]
    	indx_crr += sum(fx)
    end
    return alms
end


#= test ...
import HealpixTransforms as HT

hp   = pyimport("healpy")

Nside = 2048
lmax  = 3*(Nside)+1
l     = 0:lmax

l_mmax2, m_mmax2 = let
	mmax  = 2
	nlms  = HP.sphtfunc.Alm.getsize(lmax, mmax)
    ls = zeros(Int, nlms)
    ms = zeros(Int, nlms)
    for i in 1:nlms
        l4i, m4i = hp.sphtfunc.Alm.getlm(lmax, i-1)
        ls[i] = l4i
        ms[i] = m4i
    end
    ls, ms 
end

l, m = let
	mmax  = lmax
	nlms  = HP.sphtfunc.Alm.getsize(lmax, mmax)
    ls = zeros(Int, nlms)
    ms = zeros(Int, nlms)
    for i in 1:nlms
        l4i, m4i = hp.sphtfunc.Alm.getlm(lmax, i-1)
        ls[i] = l4i
        ms[i] = m4i
    end
    ls, ms 
end

lmat = HT.alm2triangle(l)
l′ = HT.triangle2alm(lmat)
l_mmax2′ = HT.triangle2alm(lmat[:,1:3])

sum(abs2, l′ .- l)
sum(abs2, l_mmax2′ .- l_mmax2)

=# 



# Rings 
# -----------------------------------


ring_region(i, nside) = 
        i < nside               ? "north cap" :
        nside <= i < 2nside     ? "north belt" :
        2nside == i             ? "equator" :
        2nside < i <= 3nside    ? "south belt" :
        3nside < i < 4nside     ? "south cap" : error("not a valid index for input Nside")

function isvalid_nside(nside::Int)
    ex = log(2,nside)
    return ex == Int(ex)
end


function split_rings(healpix_array::Vector{T}, nside::Int) where T<:Number

    @assert isvalid_nside(nside)

    n_rings = 4*nside - 1
    maxn_azimuth = 4*nside

    ring_vector = Vector{T}[]

    start_ring_index = end_ring_index = 0
    for i = 1:n_rings
        start_ring_index = end_ring_index + 1
        if ring_region(i,nside)=="north cap"
            end_ring_index  = start_ring_index + (4i-1)
        elseif ring_region(i,nside) ∈ ("north belt", "equator", "south belt")
            end_ring_index  = start_ring_index + (maxn_azimuth-1)
        else
            end_ring_index  = start_ring_index + (4*(n_rings-i+1)-1)
        end
        idx = start_ring_index:end_ring_index
        push!(ring_vector, healpix_array[idx])
    end

    ring_vector
end



function rings2rows(healpix_array::Vector{T}, nside::Int) where T<:Number

    @assert isvalid_nside(nside)

    n_rings = 4*nside - 1
    maxn_azimuth = 4*nside

    ring_matrix  = zeros(T, n_rings, maxn_azimuth)

    vector_of_rings = split_rings(healpix_array, nside)
    for i = 1:n_rings
        ring_i = vector_of_rings[i]
        cidx = 1:length(ring_i)
        ring_matrix[i,cidx]  = ring_i
    end

    ring_matrix
end


# function rings2rows(healpix_array::Vector{T}, nside::Int) where T<:Number

#     @assert isvalid_nside(nside)

#     n_rings = 4*nside - 1
#     maxn_azimuth = 4*nside

#     ring_matrix  = zeros(T, n_rings, maxn_azimuth)

#     start_ring_index = end_ring_index = 0
#     for i = 1:n_rings
#         start_ring_index = end_ring_index + 1
#         if ring_region(i,nside)=="north cap"
#             end_ring_index  = start_ring_index + (4i-1)
#         elseif ring_region(i,nside) ∈ ("north belt", "equator", "south belt")
#             end_ring_index  = start_ring_index + (maxn_azimuth-1)
#         else
#             end_ring_index  = start_ring_index + (4*(n_rings-i+1)-1)
#         end
#         idx = start_ring_index:end_ring_index
#         cidx = 1:length(idx)
#         ring_matrix[i,cidx]  = healpix_array[idx]
#     end

#     ring_matrix
# end



"""
θ_φ_idx_4_rings(Nside::Int) -> (θ, φ, idx, Δφ, nφ) where 
θ is a vector of the polar coordinate of each healpix ring
φ is a vector of the azimuthal coordinate of the first pixel in each ring
idx is the index of the first pixel in each ring
Δφ is a vector that records azimuthal pixel spacing in each ring
nφ is a vector that records the number of grid elements in each ring
"""
function θ_φ_idx_4_rings(nside::Int; T=Float64)

	@assert isvalid_nside(nside)

    n_rings      = 4*nside - 1
    maxn_azimuth = 4*nside

    vθ    = zeros(T, n_rings)
    vφ    = zeros(T, n_rings)
    vΔφ   = zeros(T, n_rings)
    vidx  = zeros(Int, n_rings)
    vnφ   = zeros(Int, n_rings)

    start_ring_index = end_ring_index = 0
    for i = 1:n_rings
        start_ring_index = end_ring_index + 1
        if ring_region(i,nside)=="north cap"
            end_ring_index  = start_ring_index + (4i-1)
            vθ[i] = acos(1 - i^2/nside/nside/3)
            Δφ  = π / 2 / i # 2π / nring / 2 ... citation https://arxiv.org/pdf/astro-ph/0409513.pdf
            s   = 1
            nφ  = 4i
        elseif ring_region(i,nside) ∈ ("north belt", "equator")
            end_ring_index  = start_ring_index + (maxn_azimuth-1)
            vθ[i] = acos(4/3 - 2i/nside/3)
            ## s        = mod(i-nside+1,2)
            s        = mod(i-nside,2) + 1 # from Erratum
            Δφ       = π / 2 / nside
            nφ       = 4nside
        elseif ring_region(i,nside) == "south belt"
            end_ring_index  = start_ring_index + (maxn_azimuth-1)
            vθ[i]  = acos(-(4/3 - 2*(n_rings-i+1)/nside/3))
            ## s         = mod(n_rings-i+1-nside+1,2)
            s         = mod(n_rings-i+1-nside,2) + 1  # from Erratum
            Δφ        = π / 2 / nside
            nφ        = 4nside
        else
            end_ring_index  = start_ring_index + (4*(n_rings-i+1)-1)
            vθ[i]  = acos(-(1 - (n_rings-i+1)^2/nside/nside/3))
            Δφ        = π / 2 / (n_rings-i+1) # 2π / nring / 2
            s         = 1
            nφ        = 4*(n_rings-i+1)
        end
        vidx[i] = start_ring_index
        vnφ[i]  = nφ
        vφ[i]   = Δφ * (1 - s // 2)
        vΔφ[i]  = Δφ
    end

    vθ, vφ, vidx, vΔφ, vnφ
end

θ_φ_idx_4_rings(h::Unionℍ{T}) where {T} = θ_φ_idx_4_rings(h.nside; T)

# Extract equitorial belt 
# -----------------------------------

# Fixme: this excludes the north cap and south cap boundary ... which technically are part of the equitorial blet
n_eqrings(nside::Int)    = 2nside - 1
eqring_n_pix(nside::Int) = 4nside
eqbelt_n_pix(nside::Int) = n_eqrings(nside) * eqring_n_pix(nside)

function idx_eqbelt(nside::Int)
	# number of pixels:
	#   12 Nside²
	# number of equitorial pixels:
	#	(2 Nside - 1)*(4 Nside) == 8 Nside² - 4 Nside # ... this might be (2 Nside + 1)*(4 Nside) instead
	# number of cap pixels:
	#	 12 Nside² - 8 Nside² + 4 Nside == 4 Nside² + 4 Nside
	# number of north cap pixels:
	#    2 Nside² + 2 Nside
	str_idx = 2 * nside^2 + 2 * nside + 1
	end_idx = 2 * nside^2 + 2 * nside + eqbelt_n_pix(nside)
	return CartesianIndex.(Matrix(reshape(str_idx:end_idx, eqring_n_pix(nside), n_eqrings(nside))'))
end

size_eqbelt(nside) = (n_eqrings(nside), eqring_n_pix(nside)) 


function pixfull_eqbelt(nside)
	φ    = (rw,col) -> (π/2/nside)*(col - mod(rw, 2))
	cosθ = (rw,col) -> 4//3 - 2(rw + nside - 1)//(3nside)   
	θ    = (rw,col) -> acos(cosθ(rw,col))
	nrw, ncol = size_eqbelt(nside) 
	return θ.(1:nrw, (1:ncol)'), φ.(1:nrw, (1:ncol)')
end

function pix_eqbelt(nside)
	φ    = (rw,col) -> (π/2/nside)*(col - mod(rw, 2))
	cosθ = (rw,col) -> 4//3 - 2(rw + nside - 1)//(3nside)   
	θ    = (rw,col) -> acos(cosθ(rw,col))
	nrw, ncol = size_eqbelt(nside) 
	return θ.(1:nrw, 1), φ.(1, (1:ncol)')
end

function eqbelt(hp::Vector)
	n_pix   = length(hp)
	nside   = npix2nside(n_pix)
	idx_eqb = idx_eqbelt(nside)
	ncol    = eqring_n_pix(nside)
	krng    = (0:(ncol÷2))'
	shft    = cis.(.- π .* krng ./ ncol) 
	eq 		= hp[idx_eqb]
	eqk 	= rfft(eq,(2,))
	eqk[2:2:end,:] .*= shft
	return irfft(eqk, ncol, (2,))
end

function eqbelt_2_healpix!(hpx::Vector, eq::Matrix)
	npix   = length(hpx)
	nside   = npix2nside(npix)
	idx_eqb = idx_eqbelt(nside)
	ncol    = eqring_n_pix(nside)
	krng    = (0:(ncol÷2))'
	shft    = cis.(π .* krng ./ ncol) # reverse the shift
	eqk     = rfft(eq,(2,))
	eqk[2:2:end,:] .*= shft
	hpx[idx_eqb]  = irfft(eqk,ncol,(2,)) # fixme ... it would be best to not allocate here
	hpx
end


function eqbelt_2_healpix(eq::Matrix{T}; nside::Int) where T<:Number
	hpx  = zeros(T, n_pix(nside)) 
	eqbelt_2_healpix!(hpx, eq)
end




## Spherical coordinate helpers
## -------------------------------------------------------------- 
## Polar angle $\theta \in [0,\pi]$, Azmuth angle $\varphi \in [0,2\pi]$. 


function n̂(θ,φ)
	sφ, cφ = sincos(φ)
	sθ, cθ = sincos(θ)
	nx = sθ * cφ
	ny = sθ * sφ
	nz = cθ
	return (nx, ny, nz)
end

α_n̂(n1, n2)    = acos(dot(n1,n2))

cosα_n̂(n1, n2) = dot(n1,n2)

function α_θφ(θ1,φ1,θ2,φ2) 
	Δφ = φ1 - φ2
	Δθ = θ1 - θ2
	sθ1, sθ2 = sin(θ1), sin(θ2)
	# return acos(cos(θ1)*cos(θ2) + sin(θ1)*sin(θ2)*cos(Δφ))
	return 2asin(√(sin(Δθ/2)^2 + sθ1 * sθ2 * sin(Δφ/2)^2))
end 

function cosα_θφ(θ1,φ1,θ2,φ2) 
	Δθ = θ1 - θ2
	Δφ = φ1 - φ2
	sθ1, sθ2 = sin(θ1), sin(θ2)
	return 1 - 2 * (sin(Δθ/2)^2 + sθ1 * sθ2 * sin(Δφ/2)^2)
end 



## Viz
## -------------------------------------------------------------- 


for fun ∈ (:mollview, :cartview, :azeqview)
	quote 
		function $fun(
				hpmap; 
				vmin=-maximum(abs.(hpmap)), 
				vmax=maximum(abs.(hpmap)), 
				sub=nothing,
				return_projected_map=false,
				xsize=1000,
				title="title"
			)
			hp  = pyimport("healpy") 
			hp.visufunc.$fun(
				hpmap, sub=sub, 
				min=vmin, max=vmax,
				return_projected_map=return_projected_map,
				xsize=xsize, title=title,
			)
		end
	end |> eval
end



# Rotator
# -------------------------------------------------------------- 
# z == right hand rotation about vector θ==0 
# y == right hand rotation about vector (θ,φ) = (π/2,3π/4)
function rotate_map_ZYZ(hp_map, z1::Real, y2::Real, z3::Real)
	hp  = pyimport("healpy") 
	rot = hp.rotator.Rotator(rot=(z1, y2, z3), eulertype="Y", deg=false)
	rot.rotate_map_pixel(hp_map)
end

function rotate_alm_ZYZ(hp_alm, z1::Real, y2::Real, z3::Real)
	hp  = pyimport("healpy") 
	rot = hp.rotator.Rotator(rot=(z1, y2, z3), eulertype="Y", deg=false)
	rot.rotate_map_alms(hp_map)
end


# fits reader: TODO
# -------------------------------------------------------------- 

# cmbTQU_prerot  = hp.read_map(cmb_file, verbose=false, field= (0,1,2));

# const one_K_in_mK = 1e+6
# nside = 512  |> hp.nside2npix |> hp.pixelfunc.get_min_valid_nside 
# cmb_file  = "HealpixHelper/downloads/cmbs4_06b_llcdm_f095_b23_ellmin30_map_0512_mc_0000.fits"
# eulertype = "Y" # Z rotation, then Y, then Z 
# deg   = false
# Zrot, Yrot, Xrot = deg2rad.((40, -50, -40)) # clockwise rotation
# rot = hp.rotator.Rotator(rot=(Zrot, Yrot, Xrot), eulertype=eulertype, deg=deg)
# cmbTQU_prerot  = hp.read_map(cmb_file, verbose=false, field= (0,1,2));
# obsp_prerot  = .!(UNSEEN .== cmbTQU_prerot[1,:])
# obsp   = rot.rotate_map_pixel(obsp_prerot);
# cmbTQU   = rot.rotate_map_alms(cmbTQU_prerot) .* one_K_in_mK
# cmbT = cmbTQU[1,:] .* obsp
# cmbQ = cmbTQU[2,:] .* obsp
# cmbU = cmbTQU[3,:] .* obsp
