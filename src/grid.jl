

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

# Extract equitorial belt 
# -----------------------------------

n_eqrings(nside::Int)    = 2nside - 1
eqring_n_pix(nside::Int) = 4nside
eqbelt_n_pix(nside::Int) = n_eqrings(nside) * eqring_n_pix(nside)


function idx_eqbelt(nside::Int)
	# number of pixels:
	#   12 Nside²
	# number of equitorial pixels:
	#	(2 Nside - 1)*(4 Nside) == 8 Nside² - 4 Nside 
	# number of cap pixels:
	#	 12 Nside² - 8 Nside² + 4 Nside == 4 Nside² + 4 Nside
	# number of north cap pixels:
	#    2 Nside² + 2 Nside
	str_idx = 2 * nside^2 + 2 * nside + 1
	end_idx = 2 * nside^2 + 2 * nside + eqbelt_n_pix(nside)
	return CartesianIndex.(Matrix(reshape(str_idx:end_idx, eqring_n_pix(nside), n_eqrings(nside))'))
end

size_eqbelt(nside) = (n_eqrings(nside), eqring_n_pix(nside)) 


function θφ_eqbelt(nside)
	φ    = (rw,col) -> (π/2/nside)*(col - mod(rw, 2))
	cosθ = (rw,col) -> 4//3 - 2(rw + nside - 1)//(3nside)   
	θ    = (rw,col) -> acos(cosθ(rw,col))
	nrw, ncol = size_eqbelt(nside) 
	return θ.(1:nrw, (1:ncol)'), φ.(1:nrw, (1:ncol)')
end

function θφ_eqbelt_align(nside)
	φ    = (rw,col) -> (π/2/nside)*(col - mod(rw, 2))
	cosθ = (rw,col) -> 4//3 - 2(rw + nside - 1)//(3nside)   
	θ    = (rw,col) -> acos(cosθ(rw,col))
	nrw, ncol = size_eqbelt(nside) 
	return θ.(1:nrw, 1), φ.(1, (1:ncol)')
end





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



#%% Extract equitorial belt
#%% -------------------------------------------------------------- 


function get_eq_belt(healpix_array::Array{T,d}) where {T<:Real,d}
	n_pix = size(healpix_array,1)
	nside   = npix2nside(n_pix)

	idx_eqb  = idx_eqbelt(nside)

	ncol = eqring_n_pix(nside)
	krng = (0:(ncol÷2))'
	shft = cis.(.- π .* krng ./ ncol) 

	array_of_healpix_maps = map(eachcol(healpix_array)) do fi
	    fmap = fi[idx_eqb]
	    fk = rfft(fmap,(2,))
	    fk[2:2:end,:] .*= shft
	    irfft(fk,size(fmap,2),(2,))
	end
	healpix_maps = cat(array_of_healpix_maps..., dims=(3,))

	# # this only works with FFTW not MKL
	# f = healpix_array[idx_eqb,:]
	# fk = rfft(f,2)
	# for i=1:d
	# 	fk[2:2:end,:,i] .*= cis.(.- π .* krng ./ ncol) 
	# end
	# healpix_maps = irfft(fk,ncol,2)

	θ, φ = θφ_eqbelt_align(nside)

	#Base.Slice(Base.OneTo(1))
	# Will this result be type stable??
	if d == 1
		return healpix_maps[:,:,1], θ, φ
	else
		return healpix_maps, θ, φ
	end

end
