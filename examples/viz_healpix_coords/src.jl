
# Visualizing the healpix coordinates
# ========================

using HealpixTransforms
using PyPlot

const HT = HealpixTransforms

# Set fourier grid parameters
# --------------------------

trn = let nside = 1024 # 512 # 1024 # 2048 
    ℍ0(nside, iter=1)
end


# Rotations in healpix
# -----------------------------------------

# I want to be able to move the SPT patch to the equitorial belt. Lets take a look at the rotation 
# functionality in healpix


#-
θfull, φfull  = pix(trn)
HT.mollview(θfull, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φfull, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")


# From the top plot we see that in mollview the northpole ($\theta = 0$) is shown at the top.

# The bottom plot appears counter intuitive at first glance. 
# Azmimuth $\varphi = 0$ is right at the center of the view, 
# as expected, but $\varphi \uparrow$ corresponds to moving left rather than 
# right as would be typical. To understand it I think the best way is to imagine 
# the view from an observer in the interior of the sphere ... moving "stage right" 
# from this perspective appears as moving to the left from an observer outside of the sphere. 

# Exploring rotations in healpix
# -----------------------------------------
# Rotations are specified by a sequence of three right handed rotation angles in radians.
# The convention is z,y,z rotation where: 
#   - z points to the north pole (i.e. θ=0) 
#   - y points toward (θ,φ) = (π/2,3π/4) (i.e to the right in mollview)

# Lets explore this a bit. First lets look at an initial z rotation by π/4
θrot = HT.rotate_map_ZYZ(θfull, π/4, 0, 0)
φrot = HT.rotate_map_ZYZ(φfull, π/4, 0, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")


# Now set initial z rotation to zero, but the y rotation to π/4
θrot = HT.rotate_map_ZYZ(θfull, 0, π/4, 0)
φrot = HT.rotate_map_ZYZ(φfull, 0, π/4, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")


# This can be imagined by grabbing (with your right hand, thumb pointing away from the origin) a vector emitting from the right 
# side of mollview. Then twisting by π/4 in the direction of your fingers. 

# Just to double check our understanding, if we to a zyz rotation of (0,-π/4,0) the south pole 
# should move directly up the center of view, yeilding the red hot spot moving up in the polar coordinate map.  
θrot = HT.rotate_map_ZYZ(θfull, 0, -π/4, 0)
φrot = HT.rotate_map_ZYZ(φfull, 0, -π/4, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")



# Now we can move the red hot spot in the polar map, to "stage right" by a final z rotation 
# which is negative.
θrot = HT.rotate_map_ZYZ(θfull, 0, -π/4, -π/4)
φrot = HT.rotate_map_ZYZ(φfull, 0, -π/4, -π/4)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")


# Rings to rows 
# -----------------------------------------
# HealpixTransforms.jl has functions for extracting the equitorial belt, 
# organized in a retangular array with polar angle running down the rows 
# (increasing row index -> increasing polar angle) and 
# azimuthal running across the columnes (increasing column index -> increasing azimuth)



# ## Plot rings in rows
# -----------------------------------------

# Healpix is organized in `4 Nside - 1` azimuthal rings. These rings are equi-spaced in `cos(θ)`.
# The rings can be grouped into 5 "regions":
#   - `Nside - 1` rings in the North cap (ring index: `1,...,Nside-1`)  
#   - `Nside` rings in the North equitorial belt (ring index: `Nside,...,2*Nside-1`)  
#   - `1` equator ring (ring index: `2*Nside`)  
#   - `Nside` rings in the North equitorial belt (ring index: `2*Nside+1,...,3*Nside`)  
#   - `Nside - 1` rings in the South cap (ring index: `3*Nside + 1,...,4*Nside-1`)  
#
# For easy---no processing---views of the pixel values we can put each ring into a row of a matrix
# organized with polar angle running down the rows 
# (increasing row index -> increasing polar angle) and 
# azimuthal running across the columnes (increasing column index -> increasing azimuth)

θmat = HT.rings2rows(θfull, trn.nside)

#-

φmat = HT.rings2rows(φfull, trn.nside)

#-

fig,ax = subplots(2, figsize=(14,12))
for (i,f,s) ∈ zip((1,2),(θmat, φmat), ("polar coordinates", "azimuthal coordinates"))
	img = ax[i].imshow(f)
	fig.colorbar(img, ax=ax[i])
	ax[i].set_title(s)
	fig.tight_layout()
end

