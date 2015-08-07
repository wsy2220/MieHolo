function [es_theta,es_phi,es_r] = mie_es(m,x,k,theta,phi,r,isasm)
%[es_theta,es_phi,es_r] = mie_es(m,x,k,theta,phi,r,isasm)
%
%calculate scattered field of gaven parameters
%returns three vector component at each point
%assumption: the incident field E=1
%m, x, k, and coordinate system are defined in B&H book chapter 4.
%theta,phi,and r are 2D or 1D matrices in spherical coordinates 
%which define the positions to be calculated, and must have the same size
%
%written by Shuyu Wei
%
% This code is licensed under GNU GPL V2.
if nargin < 7
	isasm = false;
end
grid_size=size(r);
nmax = round(2+x+4*x^(1/3));
n = (1:nmax)';
en = (1i).^n .* (2*n+1) ./ (n.*(n+1)) ;
n = repmat(n,[1,grid_size]);
en = repmat(en,[1,grid_size]);
kr = k * r;
kr = reshape(kr,[1,grid_size]);
kr = repmat(kr,[nmax,1,1]);
if isasm
	zn = (-1i).^n .* exp(1i*kr) ./ (1i*kr);
	krzn = zn + kr .* (-1i).^n .* exp(1i*kr) ./ kr;
else
	c2s = sqrt(0.5*pi./kr);
	zn = c2s .* besselh(n+0.5,1,kr);
	krzn = zn + 0.5 * kr .* c2s .* (besselh(n-0.5,1,kr) - besselh(n+1.5,1,kr));
end
[an,bn] = mie_abcd(m,x);
an = repmat(an,[1,grid_size]);
bn = repmat(bn,[1,grid_size]);
cos_th = cos(theta);
sin_th = sin(theta);
cos_ph = cos(phi);
sin_ph = sin(phi);
[p,t] = mie_pt(cos_th,nmax);
%cos_th = reshape(cos_th,[1,grid_size]);
sin_th = reshape(sin_th,[1,grid_size]);
cos_ph = reshape(cos_ph,[1,grid_size]);
sin_ph = reshape(sin_ph,[1,grid_size]);
%cos_th = repmat(cos_th,[nmax,1,1]);
sin_th = repmat(sin_th,[nmax,1,1]);
cos_ph = repmat(cos_ph,[nmax,1,1]);
sin_ph = repmat(sin_ph,[nmax,1,1]);

mo_theta = cos_ph .* p .* zn;
mo_phi   = -sin_ph .* t .* zn;
%me_theta = -sin_ph .* p .* zn;
%me_phi   = -cos_ph .* t .* zn;
%no_r     = sin_ph .* n .* (n+1) .* sin_th .* p .* zn ./ kr;
%no_theta = sin_ph .* t .* krzn ./ kr;
%no_phi   = cos_ph .* p .* krzn ./kr;
ne_r     = cos_ph .* n .* (n+1) .* sin_th .* p .* zn ./kr;
ne_theta = cos_ph .* t .* krzn ./ kr;
ne_phi   = -sin_ph .* p .* krzn ./kr;

es_theta = en .* (1i * an .* ne_theta - bn .* mo_theta);
es_phi   = en .* (1i * an .* ne_phi - bn .* mo_phi);
es_r     = en .* (1i * an .* ne_r);

es_theta = shiftdim(sum(es_theta,1),1);
es_phi = shiftdim(sum(es_phi,1),1);
es_r = shiftdim(sum(es_r,1),1);
