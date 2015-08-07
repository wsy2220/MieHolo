function [track_3d,fitresult] = fit2d(im, param, fit_bound)
%track_3d is the fitted positon
%fitresult is the raw result by lsqcurvefit
%param and bound are defined in config file
%
%written by Shuyu Wei
%
%This code is licensed under GNU GPL V2.

[rcenter_h, rcenter_v] = radialcenter(im);
center_h = rcenter_h * param.spacing;
center_v = rcenter_v * param.spacing;
lower = [fit_bound.radius(1), center_v-1, center_h-1, fit_bound.z(1), fit_bound.alpha(1)];%,fit_bound.ang_v(1),fit_bound.ang_h(1)];
guess = [fit_bound.radius(2), center_v, center_h,     fit_bound.z(2), fit_bound.alpha(2)];%,fit_bound.ang_v(2),fit_bound.ang_h(2)];
upper = [fit_bound.radius(3), center_v+1, center_h+1, fit_bound.z(3), fit_bound.alpha(3)];%,fit_bound.ang_v(3),fit_bound.ang_h(3)];
[h_grid,v_grid] = meshgrid((1:size(im,2))*param.spacing,(1:size(im,1))*param.spacing);
vh_grid = cat(3,v_grid,h_grid);
fitobj = @(x,vh_grid)mie_holo(param.particle_index,param.media_index,param.wavelen,x(1),vh_grid(:,:,1),vh_grid(:,:,2),...
x(2),x(3),x(4),x(5),param.polarization) ; 
fitresult = lsqcurvefit(fitobj,guess,vh_grid,im,lower,upper,param.cfoptions);
track_3d = fitresult(2:4);
