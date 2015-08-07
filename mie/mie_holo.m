function holo = mie_holo(np,ns,wavelen,radius,v_grid,h_grid,v,h,z,varargin)
%holo = mie_holo(np,ns,wavelen,radius,v_grid,h_grid,v,h,z,alpha,polar_ang,ang_v,ang_h)
%	np: partical refractive index
%	ns: solvent refractive index
%	wavelen: wavelength in vacuum
%	radius: partical radius
%	v_grid: grid of vertical coords
%	h_grid: grid of horizonal coords
%	h,v,z: partical position relative to grid coord zero
%	alpha: illumination parameter
%	polar_ang: angle of polarization in rad. 0 is vertical. 
%              Positive is rotating counterclockwise.
%	           NaN is unpolarized
%	ang_v/ang_h: rotation with axis v/h. in rad
%alpha, polar_ang, ang_v, ang_h are optional parameters
%assumption : incident light amplitiude is 1.
%
%written by Shuyu Wei
%
% This code is licensed under GNU GPL V2.
defaults = [1,0,0,0];

for ii = 1:length(varargin)
	defaults(ii) = varargin{ii};
end

alpha = defaults(1);
polar_ang = defaults(2);
ang_v = defaults(3);
ang_h = defaults(4);




wavelen_media = wavelen / ns;
k = 2 * pi / wavelen_media ;
sizeparam = k*radius; %size parameter
v_grid = v_grid - v;
h_grid = h_grid - h;
z_grid = z + tan(ang_v)*h_grid + tan(ang_h)*v_grid;
[phi_grid,theta_grid,r_grid] = cart2sph(v_grid,h_grid,z_grid);
theta_grid = 0.5*pi - theta_grid;

if isnan(polar_ang)
	phi_grid1 = phi_grid - 0.5*pi; %equivalent to rotating light
	[es_theta,es_phi,es_r] = mie_es(np/ns,sizeparam,k,theta_grid,phi_grid1,r_grid,false);
	[es_v,es_h,es_z] = sph2cartv(es_theta,es_phi,es_r,theta_grid,phi_grid1);
	es_v = es_v * alpha;
	es_h = es_h * alpha;
	es_z = es_z * alpha;
	holo1 = abs(es_v + exp(1i*k.*z_grid).*cos(ang_h)).^2 + abs(es_h.*cos(ang_v)).^2 + abs(es_z.*sin(ang_v).*sin(ang_h)).^2;
else
	holo1 = 0;
	phi_grid = phi_grid - polar_ang; %equivalent to rotating light
end

[es_theta,es_phi,es_r] = mie_es(np/ns,sizeparam,k,theta_grid,phi_grid,r_grid,false);
[es_v,es_h,es_z] = sph2cartv(es_theta,es_phi,es_r,theta_grid,phi_grid);
es_v = es_v * alpha;
es_h = es_h * alpha;
es_z = es_z * alpha;
holo = abs(es_v + exp(1i*k.*z_grid).*cos(ang_h)).^2 + abs(es_h.*cos(ang_v)).^2 + abs(es_z.*sin(ang_v).*sin(ang_h)).^2;

if isnan(polar_ang)
	holo = 0.5*(holo + holo1);
end
