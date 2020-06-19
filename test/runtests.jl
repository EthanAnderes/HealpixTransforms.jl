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

	clTfun, clEfun, clBfun, clϕfun = let
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
		clTs  = cld[:unlen_scalar][:Ctt] ./ cld[:unlen_scalar][:factor_on_cl_cmb] .|> XFields.nan2zero
		clTt  = cld[:unlen_tensor][:Ctt] ./ cld[:unlen_tensor][:factor_on_cl_cmb] .|> XFields.nan2zero
		clEs   = cld[:unlen_scalar][:Cee] ./ cld[:unlen_scalar][:factor_on_cl_cmb] .|> XFields.nan2zero
		clEt   = cld[:unlen_tensor][:Cee] ./ cld[:unlen_tensor][:factor_on_cl_cmb] .|> XFields.nan2zero
		clBt   = cld[:unlen_tensor][:Cee] ./ cld[:unlen_scalar][:factor_on_cl_cmb] .|> XFields.nan2zero
		clϕ   = cld[:phi][:Cϕϕ] ./ cld[:phi][:factor_on_cl_phi] .|> XFields.nan2zero

	    extrap = (l, Cl) -> CubicSplineInterpolation(l, Cl, extrapolation_bc = Line())
	    l_extrap = function (l, Cl)
	       	iCl = extrap(l, log.(l.^4 .* Cl))
	       	return l -> exp(iCl(l)) / l^4
	    end

	    clTfun = l_extrap(l[3]:l[end], clTs[3:end] + clTt[3:end])
	    clEfun = l_extrap(l[3]:l[end], clEs[3:end] + clEt[3:end])
	    clBfun = l_extrap(l[3]:l[end], clBt[3:end])
	    clϕfun = l_extrap(l[3]:l[end], clϕ[3:end])
	    
	    clTfun, clEfun, clBfun, clϕfun
	end


	nside = 1024 # 2048
	sph0  = ℍ0(nside, iter=0)
	sph02 = ℍ02(nside, iter=0)


	CT, CE, CB, Cϕ = let sph0 = sph0
		l, m = lm(sph0)
		CTlm = clTfun.(l)
		CElm = clEfun.(l)
		CBlm = clBfun.(l)
		Cϕlm = clϕfun.(l)
		CTlm[l .<= 2] .= 0
		CElm[l .<= 2] .= 0
		CBlm[l .<= 2] .= 0
		Cϕlm[l .<= 2] .= 0

		CT = DiagOp(Xfourier(sph0, CTlm))
		CE = DiagOp(Xfourier(sph0, CElm))
		CB = DiagOp(Xfourier(sph0, CBlm))
		Cϕ = DiagOp(Xfourier(sph0, Cϕlm))

		CT, CE, CB, Cϕ
	end


	T, E, B, ϕ, TEB = let 
		zTlm = randn(eltype_out(sph0),size_out(sph0))
		zElm = randn(eltype_out(sph0),size_out(sph0))
		zBlm = randn(eltype_out(sph0),size_out(sph0))
		zϕlm = randn(eltype_out(sph0),size_out(sph0))
		# TODO add the cross correlation ...
		T  = √CT * Xfourier(sph0, zTlm)
		E  = √CE * Xfourier(sph0, zElm)
		B  = √CB * Xfourier(sph0, zBlm)
		ϕ  = √Cϕ * Xfourier(sph0, zϕlm)

		TEB = Xfourier(sph02, hcat(T[!],E[!],B[!]))

		T, E, B, ϕ, TEB
	end

	TQUvec = TEB[:]
	
	eqQ, θ, φ = HealpixTransform.get_eq_belt(TQUvec[:,2])
	eqU,  = HealpixTransform.get_eq_belt(TQUvec[:,3])
	eqT,  = HealpixTransform.get_eq_belt(T[:])
	eqE,  = HealpixTransform.get_eq_belt(E[:])
	eqB,  = HealpixTransform.get_eq_belt(B[:])
	eqϕ,  = HealpixTransform.get_eq_belt(ϕ[:])


	eqT |> matshow; colorbar()
	eqE |> matshow; colorbar()
	eqB |> matshow; colorbar()
	eqQ |> matshow; colorbar()
	eqU |> matshow; colorbar()
	eqϕ |> matshow; colorbar()

end



## use these for testing ...
# θ, φ  = hp.pix2ang(nside, 0:(hp.nside2npix(nside)-1))
# Ωx    = abs2(hp.pixelfunc.nside2resol(nside, arcmin=false))
# npix  = hp.nside2npix(nside)
# l, m  = hp.sphtfunc.Alm.getlm(lmax)

