% Performs the fitting algorithm on "experimental" data, being located in
% "a.nameDest" folder. To apply on arbitrary experimental data, adjust
% section 1.) with correct parameters from your experiment
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
addpath('../scripts');
a = simulation;                   % load methods for calculation

%--------------------------------------------------------------------------
% 1.) Experimental parameters + fitting parameter margins
%--------------------------------------------------------------------------
a.nameDest = 'exp_data';          % folder that contains experimental data
a.psize    = 0.38e-6;             % [m] px size of detector
a.r        = 25;                  % [m] source-to-grating distance
a.a        = 6.84e-6;             % [m] grating period
a.z        = linspace(0,1.1,55);  % [m] experimental propagation distances
a.usewin   = 1;                   % apply Tukey/Hanning window function

a.h        = 3.39e-6;             % [m] height of grating structure
a.periods  = 32;                  % grating-size (in terms of periods)
a.N        = 2^13;                % number of particles (pixels)
a.plotper  = 14;      % number of periods to be plotted --> sets a.x1, a.x2

dc_min     = 0.5;                 % [1] duty cycle lower limit
dc_max     = 0.54;                % [1] duty cycle upper limit
alpha_min  = 0;                   % [°] grating angle lower limit
alpha_max  = 4;                   % [°] grating angle upper limit
sigma_min  = 0e-6;                % [m] source size lower limit
sigma_max  = 200e-6;              % [m] source size upper limit
E_min      = 14-0.5;              % [keV] x-ray energy lower limit
E_max      = 14+0.5;              % [keV] x-ray energy upper limit
s          = 1.5;                 % interval stretching factor (from paper)

n_max  = 3;                       % number of intervals to be nested
k_max  = 7;                       % number of iteration steps

%--------------------------------------------------------------------------
% 2.) Calculate Fourier coefficients from experimental data
%--------------------------------------------------------------------------
for jj=1:length(a.z)
	img     = a.loadSmallImg(jj);     % load experimental Talbot imgs
    [F_Hi F_Vi] = a.visCalc(img,jj);  % extract first Fourier coefficients
    F_exp(jj,1) = F_Hi;
    F_exp(jj,2) = F_Vi;
end

%--------------------------------------------------------------------------
% 3.) Run fitting algorithm with the above parameters
%--------------------------------------------------------------------------
fitfunc('results',a,F_exp,dc_min,dc_max,alpha_min, ...
                  alpha_max,sigma_min,sigma_max,E_min,E_max,k_max,n_max,s);