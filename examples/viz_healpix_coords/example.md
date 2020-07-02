Visualizing the healpix coordinates
========================


```julia
using HealpixTransforms
using PyPlot

const HT = HealpixTransforms
```




    HealpixTransforms



Set fourier grid parameters
--------------------------


```julia
trn = let nside = 1024 # 512 # 1024 # 2048
    ℍ0(nside, iter=1)
end
```




    ℍ0(1024, 3071, 1)



Rotations in healpix
-----------------------------------------

I want to be able to move the SPT patch to the equitorial belt. Lets take a look at the rotation
functionality in healpix


```julia
θfull, φfull  = pix(trn)
HT.mollview(θfull, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φfull, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")
```


![png](example_files/example_6_0.png)


From the top plot we see that in mollview the northpole ($\theta = 0$) is shown at the top.

The bottom plot appears counter intuitive at first glance.
Azmimuth $\varphi = 0$ is right at the center of the view,
as expected, but $\varphi \uparrow$ corresponds to moving left rather than
right as would be typical. To understand it I think the best way is to imagine
the view from an observer in the interior of the sphere ... moving "stage right"
from this perspective appears as moving to the left from an observer outside of the sphere.

Exploring rotations in healpix
-----------------------------------------
Rotations are specified by a sequence of three right handed rotation angles in radians.
The convention is z,y,z rotation where:
  - z points to the north pole (i.e. θ=0)
  - y points toward (θ,φ) = (π/2,3π/4) (i.e to the right in mollview)

Lets explore this a bit. First lets look at an initial z rotation by π/4


```julia
θrot = HT.rotate_map_ZYZ(θfull, π/4, 0, 0)
φrot = HT.rotate_map_ZYZ(φfull, π/4, 0, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")
```


![png](example_files/example_11_0.png)


Now set initial z rotation to zero, but the y rotation to π/4


```julia
θrot = HT.rotate_map_ZYZ(θfull, 0, π/4, 0)
φrot = HT.rotate_map_ZYZ(φfull, 0, π/4, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")
```


![png](example_files/example_13_0.png)


This can be imagined by grabbing (with your right hand, thumb pointing away from the origin) a vector emitting from the right
side of mollview. Then twisting by π/4 in the direction of your fingers.

Just to double check our understanding, if we to a zyz rotation of (0,-π/4,0) the south pole
should move directly up the center of view, yeilding the red hot spot moving up in the polar coordinate map.


```julia
θrot = HT.rotate_map_ZYZ(θfull, 0, -π/4, 0)
φrot = HT.rotate_map_ZYZ(φfull, 0, -π/4, 0)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")
```


![png](example_files/example_16_0.png)


Now we can move the red hot spot in the polar map, to "stage right" by a final z rotation
which is negative.


```julia
θrot = HT.rotate_map_ZYZ(θfull, 0, -π/4, -π/4)
φrot = HT.rotate_map_ZYZ(φfull, 0, -π/4, -π/4)
HT.mollview(θrot, sub=(1,2,1), vmin=0, vmax=π, title="polar coordinates")
HT.mollview(φrot, sub=(1,2,2), vmin=0, vmax=2π, title="azimuthal coordinates")
```


![png](example_files/example_18_0.png)


Rings to rows
-----------------------------------------
HealpixTransforms.jl has functions for extracting the equitorial belt,
organized in a retangular array with polar angle running down the rows
(increasing row index -> increasing polar angle) and
azimuthal running across the columnes (increasing column index -> increasing azimuth)

## Plot rings in rows
-----------------------------------------

Healpix is organized in `4 Nside - 1` azimuthal rings. These rings are equi-spaced in `cos(θ)`.
The rings can be grouped into 5 "regions":
  - `Nside - 1` rings in the North cap (ring index: `1,...,Nside-1`)
  - `Nside` rings in the North equitorial belt (ring index: `Nside,...,2*Nside-1`)
  - `1` equator ring (ring index: `2*Nside`)
  - `Nside` rings in the North equitorial belt (ring index: `2*Nside+1,...,3*Nside`)
  - `Nside - 1` rings in the South cap (ring index: `3*Nside + 1,...,4*Nside-1`)

For easy---no processing---views of the pixel values we can put each ring into a row of a matrix
organized with polar angle running down the rows
(increasing row index -> increasing polar angle) and
azimuthal running across the columnes (increasing column index -> increasing azimuth)


```julia
θmat = HT.rings2rows(θfull, trn.nside)
```




    4095×4096 Array{Float64,2}:
     0.00079736  0.00079736  0.00079736  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00159472  0.00159472  0.00159472     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00239208  0.00239208  0.00239208     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00318944  0.00318944  0.00318944     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0039868   0.0039868   0.0039868      0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00478416  0.00478416  0.00478416  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00558153  0.00558153  0.00558153     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00637889  0.00637889  0.00637889     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00717625  0.00717625  0.00717625     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00797362  0.00797362  0.00797362     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00877099  0.00877099  0.00877099  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.00956836  0.00956836  0.00956836     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0103657   0.0103657   0.0103657      0.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋮                                   ⋱       ⋮                        ⋮
     3.13202     3.13202     3.13202        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13282     3.13282     3.13282        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13362     3.13362     3.13362     …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13442     3.13442     3.13442        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13521     3.13521     3.13521        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13601     3.13601     3.13601        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13681     3.13681     3.13681        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.13761     3.13761     3.13761     …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.1384      3.1384      3.1384         0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.1392      3.1392      3.1392         0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.14        3.14        3.14           0.0  0.0  0.0  0.0  0.0  0.0  0.0
     3.1408      3.1408      3.1408         0.0  0.0  0.0  0.0  0.0  0.0  0.0




```julia
φmat = HT.rings2rows(φfull, trn.nside)
```




    4095×4096 Array{Float64,2}:
     0.785398   2.35619   3.92699   …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.392699   1.1781    1.9635       0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.261799   0.785398  1.309        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.19635    0.589049  0.981748     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.15708    0.471239  0.785398     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.1309     0.392699  0.654498  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.1122     0.336599  0.560999     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0981748  0.294524  0.490874     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0872665  0.261799  0.436332     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0785398  0.235619  0.392699     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0713998  0.214199  0.356999  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0654498  0.19635   0.327249     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0604152  0.181246  0.302076     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋮                              ⋱       ⋮                        ⋮
     0.0654498  0.19635   0.327249     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0713998  0.214199  0.356999     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0785398  0.235619  0.392699  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0872665  0.261799  0.436332     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.0981748  0.294524  0.490874     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.1122     0.336599  0.560999     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.1309     0.392699  0.654498     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.15708    0.471239  0.785398  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.19635    0.589049  0.981748     0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.261799   0.785398  1.309        0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.392699   1.1781    1.9635       0.0  0.0  0.0  0.0  0.0  0.0  0.0
     0.785398   2.35619   3.92699      0.0  0.0  0.0  0.0  0.0  0.0  0.0




```julia
fig,ax = subplots(2, figsize=(14,12))
for (i,f,s) ∈ zip((1,2),(θmat, φmat), ("polar coordinates", "azimuthal coordinates"))
	img = ax[i].imshow(f)
	fig.colorbar(img, ax=ax[i])
	ax[i].set_title(s)
	fig.tight_layout()
end
```


![png](example_files/example_24_0.png)

