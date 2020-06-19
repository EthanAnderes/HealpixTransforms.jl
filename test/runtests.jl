using HealpixTransform
using LinearAlgebra
using Test

@testset "Basic test transforming pixel white noise" begin
	
	nside = 1024 # 2048
	sph0  = ℍ0(nside, iter=0)
	sph02 = ℍ02(nside, iter=0)

	tx    = randn(n_pix(sph0)) ./ √Ωpix(sph0)
	@time tlm = sph0 * tx
	@time tx1 = sph0 \ tlm
	[tx tx1]
	@show norm(tx - tx1, 2)

	ot = HealpixTransform.∇(tlm, sph0)

	tqux   = randn(n_pix(sph02), 3) ./ √Ωpix(sph02)
	@time teblm = sph02 * tqux
	@time tqux1 = sph02 \ teblm

end



@testset "CMB simulations and equatorial belt" begin
	
	using HealpixTransform
	using LinearAlgebra
	using CMBspectra
	using Interpolations 
	using PyPlot 
	using XFields

	clTfun, clϕfun = let
		cld = CMBspectra.camb_cls(;
	        lmax    = 8050, 
	        r       = 0.1,
	        ωb      = 0.0224567, 
	        ωc      = 0.118489, 
	        τ       = 0.055, 
	        θs      = 0.0103, 
	        logA    = 3.043, 
	        ns      = 0.968602, 
	        Alens   = 1.0, 
	        k_pivot = 0.002,
	        AccuracyBoost  = 2,
	        lSampleBoost   = 4,
	        lAccuracyBoost = 4,
	    )

		l    = cld[:phi][:ell]
		clT   = cld[:unlen_scalar][:Ctt] ./ cld[:unlen_scalar][:factor_on_cl_cmb] .|> XFields.nan2zero
		clϕ   = cld[:phi][:Cϕϕ] ./ cld[:phi][:factor_on_cl_phi] .|> XFields.nan2zero

	    extrap = (l, Cl) -> CubicSplineInterpolation(l, Cl, extrapolation_bc = Line())
	    l_extrap = function (l, Cl)
	       	iCl = extrap(l, log.(l.^4 .* Cl))
	       	return l -> exp(iCl(l)) / l^4
	    end

	    clTfun = l_extrap(l[3]:l[end], clT[3:end])
	    clϕfun = l_extrap(l[3]:l[end], clϕ[3:end])
	    
	    clTfun, clϕfun
	end


	nside = 1024 # 2048
	sph0  = ℍ0(nside, iter=0)
	sph02 = ℍ02(nside, iter=0)


	CT, Cϕ = let sph0 = sph0
		l, m = lm(sph0)
		CTlm = clTfun.(l)
		Cϕlm = clϕfun.(l)
		CTlm[l .<= 2] .= 0
		Cϕlm[l .<= 2] .= 0

		CT = DiagOp(Xfourier(sph0, CTlm))
		Cϕ = DiagOp(Xfourier(sph0, Cϕlm))

		CT, Cϕ
	end


	T, ϕ = let sph0 = sph0, CT = CT, Cϕ = Cϕ
		zTx = randn(eltype_in(sph0),size_in(sph0)) ./ √Ωpix(sph0)
		zϕx = randn(eltype_in(sph0),size_in(sph0)) ./ √Ωpix(sph0)

		T = √CT * Xmap(sph0, zTx)
		ϕ = √Cϕ * Xmap(sph0, zϕx)

		T,ϕ
	end

	eqT, θ, φ = HealpixTransform.get_eq_belt(T[:])
	eqϕ,      = HealpixTransform.get_eq_belt(ϕ[:])

	eqT |> matshow; colorbar()
	eqϕ |> matshow; colorbar()

end



## use these for testing ...
# θ, φ  = hp.pix2ang(nside, 0:(hp.nside2npix(nside)-1))
# Ωx    = abs2(hp.pixelfunc.nside2resol(nside, arcmin=false))
# npix  = hp.nside2npix(nside)
# l, m  = hp.sphtfunc.Alm.getlm(lmax)

