# Modules
# =========================================

import HealpixTransforms as HT

using  EAZTransforms
using  EAZTransforms: pix, freq, nyq, Ωpix 
import EAZTransforms as EZ

using XFields 
using CMBrings: map_plot

using LBblocks: @sblock
using PyPlot

# Test pixel
# =========================================

# Set a dense EAZ grid

eaz0, region, θspan, φspan = @sblock let 

    ## --- set φ grid parameters: φspan and nφ
    φspan = deg2rad.((-180, 180)) 
    nφ    = 4096 # 18000

    ## --- set θ grid parameter: θ
    θspan, region  = (0.005, 0.785), :north_cap
    # θspan, region  = (0.785, 1.57) , :north_eq
    # θspan, region  = (1.57, 2.356) , :south_eq
    # θspan, region  = (2.356, 3.14), :south_cap
    ## --- option 1: healpix elevations
    # Nside    = 2048*2
    # θfull    = HT.θ_φ_idx_4_rings(Nside)[1]
    # ri_start = findmin(@. abs2(θfull - θspan[1]))[end]
    # ri_end   = findmin(@. abs2(θfull - θspan[2]))[end]
    # θ        = θfull[ri_start:ri_end]
    ## --- option 2: equiangle elevations
    nθ     = 4000
    θ = range(θspan[1], θspan[2], nθ)

    return EAZ0{Float64}(θ, φspan, nφ), region, θspan, φspan
end

# Choose three consequtive healpix rings to determine pixel shape

θ_center, θ_north, θ_south, Δφ_center = @sblock let region, θspan, φspan

    # ---- set parameters of the pixel
    Nside  = 64 # 128 # 256 # 512
    
    #  ---- specify within region dec
    # i.e. inside = 0.01 * Nside -> close to northern boundary of region
    # i.e. inside = 0.99 * Nside -> close to southern boundary of region
    inside = 0.1 * Nside # must be less than Nside (withing region location)

    # -----
    jrgn =  (region==:north_cap) ? 0 :
            (region==:north_eq) ? 1 :
            (region==:south_eq) ? 2 :
            (region==:south_cap) ? 3 :
            error("`region` not recognized")
    i   = round(Int,inside) + jrgn*Nside
    θ, φ, idx, Δφ, nφ = HT.θ_φ_idx_4_rings(Nside)

    return θ[i], θ[i-1], θ[i+1], Δφ[i]
end


#------

pix_field = HT.pixel.(
    EZ.θ(eaz0), EZ.φ(eaz0)';
    φ_center = 0, 
    θ_center,
    θ_north, 
    θ_south,
) |> x -> Xmap(eaz0, x)


map_plot(pix_field;  title1="Healpix pixel", vmin=-1, vmax=1)


#------
# Here we plot the pixel boundaries

figure()

φ_center = π/4
φ_left    = π/4 - Δφ_center/2
φ_right   = π/4 + Δφ_center/2
φs = range(φ_center - Δφ_center/2, φ_center + Δφ_center/2, 1000)

plot(φs, acos.(HT.zp.(φs, θ_north)), ":")
plot(φs, acos.(HT.zp.(φs, θ_south)), ":")
plot(φs, acos.(HT.zm.(φs, θ_north)), "-")
plot(φs, acos.(HT.zm.(φs, θ_south)), "--")

plot(φ_left,   θ_center, "*")
plot(φ_center, θ_north, "o")
plot(φ_center, θ_south, "x")
plot(φ_right,  θ_center, "*")

#------
# compare with the formulas given in the healpix paper
# in the north cap

# vertices of pixel boundaries
θv(k,j,Nside) = acos(1 - ((k+j) / Nside)^2 / 3)
φv(k,j,Nside) = π * k / (k+j) / 2

Nside = 2*2048
k = 10   # North cap
@assert 1 ≤ k ≤ (Nside-1)÷2 # 1 ≤ 2k ≤ Nside-1 ⟹ 1 ≤ k ≤ (Nside-1)÷2
θ_north  = θv(k,k,Nside)
θ_south  = θv(k+1,k+1,Nside)
θ_center = θv(k,k+1,Nside) # same as θv(k+1,k,Nside)
φ_center = φv(k,k,Nside) # ≈ π / 4
φ_left   = φv(k,k+1,Nside)
φ_right  = φv(k+1,k,Nside)
Δφ_center = φ_right - φ_left
φs = range(φ_left - Δφ_center/2, φ_right + Δφ_center/2, 1000)

figure()
# healpix paper equations for the boundary lines
# ... only for north cap
zp_test(k,φ,Nside) = 1 - (k * (π / 2 / φ) / Nside)^2 / 3
zm_test(k,φ,Nside) = 1 - (k * (π / 2 / (π / 2 - φ)) / Nside)^2 / 3
plot(φs, acos.(zp_test.(k,φs,Nside)))
plot(φs, acos.(zp_test.(k+1,φs,Nside)))
plot(φs, acos.(zm_test.(k,φs,Nside)))
plot(φs, acos.(zm_test.(k+1,φs,Nside)))

# HealpixTransform.jl boundary lines
plot(φs, acos.(HT.zp.(φs, θ_north)), ":")
plot(φs, acos.(HT.zp.(φs, θ_south)), ":")
plot(φs, acos.(HT.zm.(φs, θ_north)), ":")
plot(φs, acos.(HT.zm.(φs, θ_south)), ":")

plot(φ_left,   θ_center, "*")
plot(φ_center, θ_north, "o")
plot(φ_center, θ_south, "x")
plot(φ_right,  θ_center, "*")


