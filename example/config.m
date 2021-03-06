%sample configuration for mie theory fitting
%all fitting bound are given as [lower, guess, upper].
%currently ang_h, ang_v, particle_index are not fitted,
%but they are included here for future changes.
%you can modify fit1d.m and fit2d.m yourself to choose what to fit.
%
%written by Shuyu Wei
%
%This code is licensed under GNU GPL V2.

spacing = 0.1375;   %width of a sigle pixel in hologram 
wavelen = 0.66;  %wave length of incident light in vacuum
radius = 1;      %particle radius, currently not used.
media_index = 1.33;       %media refractive index
particle_index = 1.4891;   %particle refractive index, 
polarization = NaN; %NaN for unpolarized
z_bound = [50,55,60]; %[lower, guess, upper];
alpha_bound = [0.6,0.8,1.0];
radius_bound = [1.05,1.15,1.2];
particle_index_bound = [0.5,1,1.5];
ang_h_bound = deg2rad([-12,0,12]);
ang_v_bound = deg2rad([-12,0,12]);
cfoptions = optimoptions(@lsqcurvefit,'Display','none','TolFun',1e-10,'TolX',1e-10);
param = struct( ...
        'spacing',spacing, ...
		'wavelen',wavelen, ...
		'radius',radius, ...
		'particle_index',particle_index, ...
		'media_index',media_index, ...
		'media_wavelen',wavelen/media_index, ...
		'cfoptions',cfoptions,...
		'polarization',polarization);

fit_bound = struct( ...
		'z',z_bound, ...
		'alpha',alpha_bound, ...
		'radius',radius_bound, ...
		'particle_index',particle_index_bound, ...
		'ang_h',ang_h_bound, ...
		'ang_v',ang_v_bound); 
