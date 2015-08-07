clear;
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
