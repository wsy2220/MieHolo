function [x,y,z]=sph2cartv(vtheta,vphi,vr,theta,phi)
%[x,y,z]=sph2cartv(vtheta,vphi,vr,theta,phi)
%Convert  spherical vectors to cartesian vectors
%all parameters must have the same size
%
%written by Shuyu Wei
%
% This code is licensed under GNU GPL V2.
x = vr .* sin(theta) .* cos(phi) ...
     + vtheta .* cos(theta) .* cos(phi) ...
	 - vphi .* sin(phi);
y = vr .* sin(theta) .* sin(phi) ...
     + vtheta .* cos(theta) .* sin(phi) ...
	 + vphi .* cos(phi);
z = vr .* cos(theta) - vtheta .* sin(theta);
