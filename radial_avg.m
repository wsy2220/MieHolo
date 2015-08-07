function [avg,count,r_range] = radial_avg(im, center_v, center_h, th_num)
%[avg,count] = radial_avg(im, center_v, center_h, th_num)
%calculate radial average over theta
%all units are pixel
%
%written by Shuyu Wei
%
% This code is licensed under GNU GPL V2.

if nargin < 4
	th_num = sum(size(im))*2;
end
step = 1;
[v,h]=size(im);
center = [center_v,center_h];
sample_num = max([norm(center),norm(center-[v,1]),norm(center-[v,h]),norm(center-[1,h])]);
sample_num = fix(sample_num)-8;
r_range = 0:step:sample_num;
th_range = linspace(0,2*pi,th_num+1)';
th_range = th_range(1:end-1);
grid_sh = center_h + cos(th_range) * r_range;
grid_sv = center_v - sin(th_range) * r_range;
inside = (grid_sh>1) & (grid_sh<h) & (grid_sv>1) & (grid_sv<v);
index = find(inside);
grid_res = interp2(im,grid_sh(index),grid_sv(index),'spline');
count = sum(inside);
avg = zeros(size(inside));
avg(index) = grid_res;
avg = sum(avg)./count;
