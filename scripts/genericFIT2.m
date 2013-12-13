% generic fitting procedure for multilayer evaluation
% TODO: obsolet (remove redunacy when changing genericFIT.m)
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
z   = linspace(0.011,1.111,111);
a.r = 26.55;                 % [m] radius of incidient wave curvature

% Fitting parameters
a.a       = 6.84e-6;               % [m] grating period
a.periods = 32;                    % grating-size (in terms of periods)
a.N       = 2^13;                  % number of particles --> 2^n !!!
a.padding = 1*a.N;                 % total number (with zero-padding)
a.gHeight = 3.39e-6;               % height of grating structure

% construct fourCoeffexp from experimental data
fourCoeffexp = [horz.' vert.'];


% starting values
stepsz  = 3;
iterat  = 7;
dutDN   = 0.5;
dutUP   = 0.54;
angDN   = 0;
angUP   = 4.2;
src_DN = 1e-6;
src_UP = 200e-6;
range_sz = 1.5;
RR_DNh = 1;
RR_UPh = 25;
RR_DNv = RR_DNh;
RR_UPv = RR_UPh;

% change: instead of varying energy, the source-to-sample distance is
% varied!
fitpar = fitfunc2(a,b,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP, ...
                 src_DN,RR_UPh,RR_DNh,RR_UPv,RR_DNv,iterat,stepsz,range_sz,z, ...
                 [sample{iii} 'kev_simFIT']);
end