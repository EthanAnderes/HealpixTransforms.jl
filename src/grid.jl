
npix2nside(n_pix::Int) = Int(sqrt(n_pix/12))

n_pix(h::Unionℍ)  = n_pix(h.nside)
n_pix(nside::Int) = 12nside^2

n_rings(h::Unionℍ)  = n_rings(h.nside)
n_rings(nside::Int) = 4nside - 1

n_lm(h::Unionℍ) = n_lm(h.lmax)
n_lm(lmax::Int) = lmax * (lmax + 1) ÷ 2 + lmax + 1

Ωpix(h::Unionℍ) = Ωpix(h.lmax)
Ωpix(nside::Int) = 4π / n_pix(nside)

function pix(h::Unionℍ)
	θ, φ = pix(h.nside)
	return θ, φ
end
function pix(nside::Int)
	hp  = pyimport("healpy") 
	θ, φ  = hp.pix2ang(nside, 0:(n_pix(nside)-1))
	return θ, φ
end

lm(h::Unionℍ) = lm(h.lmax)
function lm(lmax::Int) 
	nlm = n_lm(lmax)
	i = 1:nlm 
	m = @. ceil(Int, (2lmax + 1)/2 - √((2lmax + 1)^2 - 8*(i - 1 - lmax))/2 )
	l = @. i - 1 - m*(2lmax + 1 - m) ÷ 2
	return l, m 
end



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

function rings2rows(healpix_array::Vector{T}, nside::Int) where T<:Number

    @assert isvalid_nside(nside)

    n_rings = 4*nside - 1
    maxn_azimuth = 4*nside

    ring_matrix  = zeros(T, n_rings, maxn_azimuth)

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
        cidx = 1:length(idx)
        ring_matrix[i,cidx]  = healpix_array[idx]
    end

    ring_matrix
end


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

function eqbelt_2_healpix!(hp::Vector, eq::Matrix)
	n_pix   = length(hp)
	nside   = npix2nside(n_pix)
	idx_eqb = idx_eqbelt(nside)
	ncol    = eqring_n_pix(nside)
	krng    = (0:(ncol÷2))'
	shft    = cis.(π .* krng ./ ncol) # reverse the shift
	eqk     = rfft(eq,(2,))
	eqk[2:2:end,:] .*= shft
	hp[idx_eqb]  = irfft(eqk,ncol,(2,)) # fixme ... it would be best to not allocate here
	hp
end


# TODO: removing this in favor of taking a single healpix argument
# function eqbelt(healpix_array::Array{T,d}) where {T<:Real,d}
# 	n_pix = size(healpix_array,1)
# 	nside   = npix2nside(n_pix)

# 	idx_eqb  = idx_eqbelt(nside)

# 	ncol = eqring_n_pix(nside)
# 	krng = (0:(ncol÷2))'
# 	shft = cis.(.- π .* krng ./ ncol) 

# 	array_of_healpix_maps = map(eachcol(healpix_array)) do fi
# 	    fmap = fi[idx_eqb]
# 	    fk = rfft(fmap,(2,))
# 	    fk[2:2:end,:] .*= shft
# 	    irfft(fk,size(fmap,2),(2,))
# 	end
# 	healpix_maps = cat(array_of_healpix_maps..., dims=(3,))

# 	# # this only works with FFTW not MKL
# 	# f = healpix_array[idx_eqb,:]
# 	# fk = rfft(f,2)
# 	# for i=1:d
# 	# 	fk[2:2:end,:,i] .*= cis.(.- π .* krng ./ ncol) 
# 	# end
# 	# healpix_maps = irfft(fk,ncol,2)

# 	#Base.Slice(Base.OneTo(1))
# 	# Will this result be type stable??
# 	if d == 1
# 		return healpix_maps[:,:,1]
# 	else
# 		return healpix_maps
# 	end

# end





#%% Spherical coordinate helpers
#%% -------------------------------------------------------------- 
#%% Polar angle $\theta \in [0,\pi]$, Azmuth angle $\varphi \in [0,2\pi]$. 


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



#%% Viz
#%% -------------------------------------------------------------- 


for fun ∈ (:mollview, :cartview, :orthview, :gnomview)
	quote 
		function $fun(
				hpmap; 
				vmin=-maximum(abs.(hpmap)), 
				vmax=maximum(abs.(hpmap)), 
				sub=nothing,
				xsize=1000, 
				title="title"
			)
			hp  = pyimport("healpy") 
			hp.visufunc.$fun(hpmap,sub=sub,min=vmin,max=vmax,xsize=xsize,title=title)
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
