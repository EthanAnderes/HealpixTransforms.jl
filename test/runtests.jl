using HealpixTransformss
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

	ot = HealpixTransforms.∇(tlm, sph0)

	tqux   = randn(n_pix(sph02), 3) ./ √Ωpix(sph02)
	@time teblm = sph02 * tqux
	@time tqux1 = sph02 \ teblm

end



## use these for testing ...
# θ, φ  = hp.pix2ang(nside, 0:(hp.nside2npix(nside)-1))
# Ωx    = abs2(hp.pixelfunc.nside2resol(nside, arcmin=false))
# npix  = hp.nside2npix(nside)
# l, m  = hp.sphtfunc.Alm.getlm(lmax)

