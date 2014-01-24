% Creates synthetic Talbot images in folder "exp_data"for demonstrating
% the results of the proposed experimental setup.
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');

a = simulation;                   % load methods for calculation
if isdir('exp_data')
    rmdir('exp_data','s');        % delete previously created files
end
a.nameDest = 'exp_data';          % destination folder for saving images

%--------------------------------------------------------------------------
% 1.) Parameters for synthetic data
%--------------------------------------------------------------------------
a.a       = 6.84e-6;              % [m] grating period
a.r       = 25;                   % [m] radius of incidient wave curvature
a.E       = 14;                   % [keV] x-ray energy
a.srcsz   = [133e-6 52e-6];       % [m] TOMCAT src sizes
a.periods = 48;                   % grating-size (in terms of periods)
a.N       = 2^12;                 % number of particles (pixels)
a.psize   = 0.38e-6;              % [m] px size of detector

% Grating parameters
a.gHeight = 3.39e-6;              % [m] height of grating structure
a.gAngle  = 1.5;                  % [°] angle of grating's bumps slope
a.duty    = 0.52;                 % [1] duty cycle

z = linspace(0,1.1,55);           % [m] propagation distances in meters

%--------------------------------------------------------------------------
% 2.) Construct grating
%--------------------------------------------------------------------------
gra = a.talbotGrid2D;

%--------------------------------------------------------------------------
% 3.) Wave field @ grating
%--------------------------------------------------------------------------
a.calcRefracAbsorb(17); % calculate delta & beta for grating with rho = 17
f = a.waveFieldGrat(gra);

%--------------------------------------------------------------------------
% 4.) Propagation along z-axis
%--------------------------------------------------------------------------
[x y] = meshgrid(round(-10/2):round(10/2), round(-10/2):round(10/2));
gauss = exp( (-x.^2./2) - (y.^2./2) );
gauss = gauss./sum(gauss(:));

for ii=1:length(z)
	calc = a.waveFieldPropMutual2D(z(ii),f); % wave propagation
	crop = a.scale2Det(calc);                % scale to detector pixel size
	crop = conv2(crop,gauss,'same');         % simulate detector's PSF
	crop = crop./4;                          % correct range in Talbot imgs
	crop = crop(70:795,70:795);              % correct for border areas
	a.Save2img(crop,ii);                     % save images
end