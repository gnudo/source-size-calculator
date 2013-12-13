% sends the used parameters to the fitting algorithm and contains mostly
% experiment-relevant information...(will be changed in the future, and probably
% obsolete)
%--------------------------------------------------------------------------
% Author: Goran Lovric
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
a = simulation;               % load classes for calculation
b = analysis;
sample = {'G_14' 'G_21_8' 'GWSi_21_8' 'GSi111_22' 'G_15' 'G_16' 'G_17' 'ML2_Si111_4x_'};
enes = [14 21.8 21.8 22 15 16 17 18];

for iii = 8:8
iii
% Experiment-specific settings
a.E    = enes(iii);
load(['../experiment/' sample{iii} 'kev_exp']);
if iii == 1
    z   = linspace(0.002,1.102,111);
    a.r = 26.3;                 % [m] radius of incidient wave curvature
elseif iii == 2 | iii == 3 |  iii == 4
    z1  = linspace(0.002,1.102,56);     % [m] propagation distance in meters
    z2  = linspace(1.002,2.102,56);     % [m] propagation distance in meters
    z   = [z1(1:end-6) z2];
    a.r = 26.3;                 % [m] radius of incidient wave curvature
else
    z   = linspace(0.011,1.111,111);
    a.r = 26.55;                 % [m] radius of incidient wave curvature
end

% Fitting parameters
a.a       = 6.84e-6;               % [m] grating period
a.periods = 32;                    % grating-size (in terms of periods)
a.N       = 2^13;                  % number of particles --> 2^n !!!
a.padding = 1*a.N;                 % total number (with zero-padding)
a.gHeight = 3.39e-6;               % height of grating structure

% construct fourCoeffexp from experimental data
if iii == 2 | iii == 3 | iii == 4
    horz_tmp = [horz(1,1:end-6) horz(2,:)];
    vert_tmp = [vert(1,1:end-6) vert(2,:)];
    fourCoeffexp = [horz_tmp.' vert_tmp.'];
else
    fourCoeffexp = [horz.' vert.'];
end

% starting values
stepsz  = 3;
iterat  = 7;
dutDN   = 0.5;
dutUP   = 0.54;
angDN   = 0;
angUP   = 4.2;
src_DN = 0e-6;
src_UP = 200e-6;
range_sz = 1.5;
ene_DN = a.E-0.5;
ene_UP = a.E+0.5;

fitpar = fitfunc(a,b,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP, ...
                 src_DN,ene_UP,ene_DN,iterat,stepsz,range_sz,z, ...
                 [sample{iii} 'kev_simFIT']);
end