#MieHolo

MieHolo offers some matlab functions to generate and fit holograms
that produced by a single sphere under plane incident light.

##Usage
First add root folder and `mie/` to matlab path.

Typically there are 3 functions to be used:

* mie_holo
* fit1d
* fit2d

Example codes can be found  in `example/` folder.

###Generate holograms

`mie_holo` is used to generate holograms.
You should first generate a meshgrid to define a image plane. 
See the example below.
v and h are vertical and horizonal coordinates relative to the 
origin of the image plane.

To generate a 201*201 hologram with pixel size 0.2um, and a silica sphere
(n=1.456, radius=1) located at 30um above the focal plane in water:

```matlab
np = 1.456;
ns = 1.33;
wavelen = 0.66;
radius = 1;
pixsize = 0.2;
v_range = linspace(0, 200*pixsize, 201); %grid origin is in upper left conner
h_range = linspace(0, 200*pixsize, 201); %point values are in actual unit
v = 20;
h = 20;
z = 30;
[h_grid, v_grid] = meshgrid(h_range, v_range);
holo = mie_holo(np,ns,wavelen,radius,v_grid,h_grid,v,h,z);
imshow(holo,[]);
```
You can pass additional parameters to `mie_holo` to change incident light
properties, see `mie/mie_holo.m` for more details.

###Fit holograms

'fit1d' and 'fit2d' are used to fit experimental holograms to mie solutions
to extract characteristics of a sphere.

To use these functions, you need to run the configuration code to set parameters.
Copy `config_sample.m` and edit parameters in the copied file,
then pass hologram and parameters to the function.

See `example/fit_holo.m`

##References
1. Lee. Characterizing and tracking single colloidal particles with video holographic microscopy Opt. Express, OSA, 2007.
2. Bohren & Huffman. Absorption and Scattering by a Sphere Absorption and Scattering of Light by Small Particles.

##License
This code is licensed under GNU GPL V2 license.
