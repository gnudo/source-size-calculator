% creates synthetic Talbot images in folder source-size-calculator/exp_data
% for demonstrating the results of the proposed experimental setup
%--------------------------------------------------------------------------
% Date: 2013-12-19
% Author: Goran Lovric
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');

a = simulation;               % load methods for calculation
b = analysis;                 % load method for saving images
if isdir('../exp_data')
    rmdir('../exp_data','s')  % delete previously created files
end
b.nameDest = '../exp_data';   % destination folder for saving images

%--------------------------------------------------------------------------
% 1.) CONSTANTS
%--------------------------------------------------------------------------
a.a       = 6.84e-6;               % [m] grating period
a.r       = 25;                    % [m] radius of incidient wave curvature
a.E       = 14;                    % [keV]
a.srcsz   = [133e-6 52e-6];        % TOMCAT src sizes

a.periods = 8;                     % grating-size (in terms of periods)
a.N       = 2^11;                  % number of particles --> 2^n !!!
a.padding = 1*a.N;                 % total number (with zero-padding)

% Grating parameters
a.gHeight = 3.39e-6;               % height of grating structure
a.gAngle  = 1.5;                   % angle of bump's slope
a.duty    = 0.518;                 % duty cycle

z = linspace(0,1.1,110);           % [m] propagation distance in meters

%--------------------------------------------------------------------------
% 2.) GRID creation
%--------------------------------------------------------------------------
gra = a.talbotGrid2D;

%--------------------------------------------------------------------------
% 3.) Wave field @ Grid
%--------------------------------------------------------------------------
f = a.waveFieldGrat(gra);

%--------------------------------------------------------------------------
% 4.) Propagation along z-axis
%--------------------------------------------------------------------------

for ii=1:length(z)
    calc = a.waveFieldPropMutual2D(z(ii),f);
    crop = a.scale2Det(calc,0.38e-6);     % scale to detector resolution
    crop = crop + rand(size(crop)).*(max(crop(:))-min(crop(:))).*0.1;
                                          % add here 10% noise
    b.Save2img(crop,ii);
end