% Creates synthetic Talbot images in folder "exp_data_refl_beam" for
% demonstrating the results of the multilayer (ML) characterization
% experiment according to the setup of Rack et al., JSR 17, 496-510 (2010).
%--------------------------------------------------------------------------
% Date: 2014-01-24
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');

a = simulation;                      % load methods for calculation
if isdir('exp_data_refl_beam')
    rmdir('exp_data_refl_beam','s'); % delete previously created files
end
a.nameDest = 'exp_data_refl_beam';   % destination folder for saving images

%--------------------------------------------------------------------------
% 1.) Parameters for synthetic data
%--------------------------------------------------------------------------
a.a       = 6.84e-6;              % [m] grating period
a.R       = 25.0;                 % [m] radius of incidient wave curvature
a.R2      = 5.0;                  % [m] radius of vertical wave curvature
a.E       = 18;                   % [keV] x-ray energy
a.sigma   = [154e-6 43e-6];       % [m] TOMCAT src sizes
a.periods = 48;                   % grating-size (in terms of periods)
a.N       = 2^12;                 % number of particles (pixels)
a.psize   = 0.76e-6;              % [m] px size of detector

% Grating parameters
a.h       = 3.39e-6;              % [m] height of grating structure
a.alpha   = 1.5;                  % [°] angle of grating's bumps slope
a.dc      = 0.52;                 % [1] duty cycle

a.z = linspace(0,1.1,56);         % [m] propagation distances in meters

%--------------------------------------------------------------------------
% 2.) Construct grating
%--------------------------------------------------------------------------
gra = a.talbotGrid2D;

%--------------------------------------------------------------------------
% 3.) Wave field @ grating
%--------------------------------------------------------------------------
a.calcRefracAbsorb('Au',17); % calculate delta,beta for grating with rho=17
f = a.waveFieldGrat(gra);

%--------------------------------------------------------------------------
% 4.) Propagation along z-axis
%--------------------------------------------------------------------------
[x y] = meshgrid(round(-10/2):round(10/2), round(-10/2):round(10/2));
gauss = exp( (-x.^2./2) - (y.^2./2) );
gauss = gauss./sum(gauss(:));

for ii=1:length(a.z)
	calc = a.waveFieldPropMutual2D(a.z(ii),f); % wave propagation
	crop = a.scale2Det(calc);                % scale to detector pixel size
	crop = conv2(crop,gauss,'same');         % simulate detector's PSF
	crop = crop./4;                          % correct range in Talbot imgs
	crop = crop(33:400,33:400);              % correct for border areas
	a.Save2img(crop,ii);                     % save images
end