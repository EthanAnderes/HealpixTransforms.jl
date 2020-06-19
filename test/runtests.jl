using HealpixTransform
using Test

@testset "HealpixTransform.jl" begin
    # Write your tests here.
end


## use these for testing ...
# θ, φ  = hp.pix2ang(nside, 0:(hp.nside2npix(nside)-1))
# Ωx    = abs2(hp.pixelfunc.nside2resol(nside, arcmin=false))
# npix  = hp.nside2npix(nside)
# l, m  = hp.sphtfunc.Alm.getlm(lmax)
