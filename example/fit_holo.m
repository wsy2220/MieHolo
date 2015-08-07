clear;
load sample;
config;
position_1d = fit1d(imr, param, fit_bound)
position_2d = fit2d(imr, param, fit_bound)
