using FFTW 
# for some reason using I need `using FFTW` *before* `using HealpixHelper`
# to get the equatorial belt test set to pass with mkl for FFTW
using HealpixHelper
using Test
using PyCall 
const HH = HealpixHelper
const hp = pyimport("healpy") 

@testset "pixel count and indexing" begin
	let	n4 = Nside(4)
		@test HH.n_pix(n4) == 192

		@test HH.idx_eqbelt(n4)[1,1]     == CartesianIndex(41,)
		@test HH.idx_eqbelt(n4)[end,1]   == CartesianIndex(137,)
		@test HH.idx_eqbelt(n4)[end,end] == CartesianIndex(152,)
		@test HH.idx_eqbelt(n4)[1,end]   == CartesianIndex(56,)
		
		θn4, φn4 = HH.θφ_eqbelt(n4)

		@test φn4[1,1] == 0.0
		@test φn4[2,1] > 0.0
	end
end




@testset "spherical coord helpers" begin
	let θ1=π*rand(), θ2=π*rand(), φ1=2π*rand(), φ2=2π*rand()
		@test sum(abs.(HH.n̂(π/2,0) .- (1,0,0)))≈0.0 atol=1e-10

		n1 = HH.n̂(θ1,φ1)
		n2 = HH.n̂(θ2,φ2)

		@test HH.α_n̂(n1, n2) ≈ HH.α_θφ(θ1,φ1,θ2,φ2)       atol=1e-10
		@test HH.cosα_n̂(n1, n2)       ≈ cos(HH.α_n̂(n1, n2)) atol=1e-10
		@test HH.cosα_θφ(θ1,φ1,θ2,φ2) ≈ cos(HH.α_n̂(n1, n2)) atol=1e-10
	end
end



@testset "extract equatorial belt" begin
	let nsd=Nside(256)
		hp_θ, hp_φ  = hp.pix2ang(nsd.s, 0:(hp.nside2npix(nsd.s)-1))
		@test nsd == HH.npix2nside(length(hp_φ))

		ftst = φ -> sin(φ * 10)
		gtst = φ -> cos(φ * 10)

		hp_map  = ftst.(hp_φ)
		hp_2map = hcat(ftst.(hp_φ), gtst.(hp_φ))
		
		eq_map, eq_θ, eq_φ = get_eq_belt(hp_map)
		eq_2map,           = get_eq_belt(hp_2map)
		
		# eq_map |> matshow
		# eq_2map[:,:,1] |> matshow
		# eq_2map[:,:,2] |> matshow

		@test maximum(abs.(eq_map .- ftst.(eq_φ))) ≈ 0.0 atol=1e-10
		@test maximum(abs.(eq_2map[:,:,1] .- ftst.(eq_φ))) ≈ 0.0 atol=1e-10
		@test maximum(abs.(eq_2map[:,:,2] .- gtst.(eq_φ))) ≈ 0.0 atol=1e-10
	end
end


