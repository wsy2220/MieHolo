function [track_3d, fitresult] = fit1d(im, param, fit_bound,rcenter_v,rcenter_h)
%track_3d = fit1d(im, param, fit_bound)
%track_3d is the fitted positon
%fitresult is the raw result by lsqcurvefit
%param and bound are defined in config file
%
%written by Shuyu Wei
%
%This code is licensed under GNU GPL V2.
if nargin~=5
	[rcenter_h, rcenter_v] = radialcenter(im);
end
center_h = rcenter_h * param.spacing;
center_v = rcenter_v * param.spacing;
[avg,count,r_range] = radial_avg(im,rcenter_v,rcenter_h);
lower = [fit_bound.radius(1), fit_bound.z(1), fit_bound.alpha(1)];
guess = [fit_bound.radius(2), fit_bound.z(2), fit_bound.alpha(2)];
upper = [fit_bound.radius(3), fit_bound.z(3), fit_bound.alpha(3)];
weight = sqrt(2*pi*(0:length(avg)-1).*count/count(1));
r_line = r_range * param.spacing;
fitobj = @(x,r_line)mie_holo(param.particle_index,param.media_index,param.wavelen,x(1),0,r_line,...
0,0,x(2),x(3),param.polarization) .* weight; 
fitresult = lsqcurvefit(fitobj,guess,r_line,avg.*weight,lower,upper,param.cfoptions);
center_z = fitresult(2);
track_3d = [center_v, center_h, center_z];
