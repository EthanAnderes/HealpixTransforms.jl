# Healpix pixel indicator
# ====================================

"""
Indicator for a healpix pixel, shifted to have center at `φ_center`

```julia
function pixel(
    θ, φ;
    θ_center,
    θ_north, 
    θ_south,
    φ_center 
) -> true or false
```

The pixel is defined by `θ_center, θ_north, θ_south` which correspond 
to the polar elevation of the center of the pixel, the northern tip of 
the pixel and the southern tip of the pixel. Note: `θ_center, θ_north, 
θ_south` need to correspond to polar elevation of three  consecutive 
healpix rings. The spacing of these rings encodes the nside of the pixelation. 

The keyword variable `φ_center` shifts the healpix pixel in the 
azimuthal coordinate so the coordinates of the pixel scenter 
are `(θ_center, φ_center)`.

Warning: In the polar caps the diamond shape of the pixels on each 
ring changes slightly with aziumth. The method `pixel` uses the same 
diamond shape on each  healpix ring which approximates the healpix 
pixles centered near `φ=π/4` in the polar regions.

"""
function pixel(
    	θ, φ;
    	θ_center,
    	θ_north, 
    	θ_south,
    	φ_center, 
        Δφ_center, # strickly not necessary but makes things easier 
	)

    if angle_separation(φ - φ_center) ≥ Δφ_center 
        return false
    elseif (θ < θ_north) | (θ > θ_south)
        return false
    else
        φsft = mod(φ - φ_center + π/4, 2π)
        zp₀  = zp(φsft, θ_north)  
        zp₁  = zp(φsft, θ_south)
        zm₀  = zm(φsft, θ_north)  
        zm₁  = zm(φsft, θ_south)  
        @assert zp₁ ≤ zp₀
        @assert zm₁ ≤ zm₀
        zθ = cos(θ)
        return (zp₁ ≤ zθ ≤ zp₀) && (zm₁ ≤ zθ ≤ zm₀)
    end

end

# Helper functions


# This formula for angle_separation is the geodesic formula simplified 
# for two points (θ₁, φ₁) & (θ₂, φ₂) on the sphere at the equator
# where Δφ = φ₂ - φ₁, θ₁ = θ₂ = π/2
# 
# copied from CirculantCov.jl
#
angle_separation(Δφ)  = 2asin(abs(sin(Δφ/2)))


zp_Cap(φ, θₒ) = 1 - (1-cos(θₒ)) * (π/4/φ)^2 
zm_Cap(φ, θₒ) = 1 - (1-cos(θₒ)) * (π/4/(π/2 - φ))^2 
zp_Equ(φ, θₒ) = cos(θₒ) + (φ*8/3/π - 2/3)
zm_Equ(φ, θₒ) = cos(θₒ) - (φ*8/3/π - 2/3)

function zp(φ, θₒ)
    zₒ = cos(θₒ)
    if zₒ > 2/3
        return zp_Cap(φ, θₒ)
    elseif 0 ≤ zₒ ≤ 2/3
        return zp_Equ(φ, θₒ)
    elseif -2/3 ≤ zₒ < 0
        return - zm_Equ(φ, π-θₒ)
    else
        @assert zₒ < -2/3
        return - zm_Cap(φ, π-θₒ)
    end
end

function zm(φ, θₒ)
    zₒ = cos(θₒ)
    if zₒ > 2/3
        return zm_Cap(φ, θₒ)
    elseif 0 ≤ zₒ ≤ 2/3
        return zm_Equ(φ, θₒ)
    elseif -2/3 ≤ zₒ < 0
        return - zp_Equ(φ, π-θₒ)
    else
        @assert zₒ < -2/3
        return - zp_Cap(φ, π-θₒ)
    end
end
