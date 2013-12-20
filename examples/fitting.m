% performs the fitting algorithm on the synthetic "experimental" data
%--------------------------------------------------------------------------
% Date: 2013-12-19
% Author: Goran Lovric
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
addpath('../scripts');

a = simulation;               % load methods for calculation
b = analysis;                 % load method for saving images

%--------------------------------------------------------------------------
% 1.) Experimental parameters + fitting parameter margins
%--------------------------------------------------------------------------
a.E       = 14;
a.r       = 25;
a.a       = 6.84e-6;          % [m] grating period
a.periods = 32;               % grating-size (in terms of periods)
a.N       = 2^13;             % number of particles --> 2^n !!!
a.padding = 1*a.N;            % total number (with zero-padding)
a.gHeight = 3.39e-6;          % height of grating structure

dutDN   = 0.5;
dutUP   = 0.54;
angDN   = 0;
angUP   = 4.2;
src_DN = 0e-6;
src_UP = 200e-6;
range_sz = 1.5;
ene_DN = a.E-0.5;
ene_UP = a.E+0.5;

stepsz  = 3;
iterat  = 4;

%--------------------------------------------------------------------------
% analyze experimental data
%--------------------------------------------------------------------------
%b.nameData = 'exp_data';     % set <b.src>-folder with img-tifs
b.nameDest = ['../' 'exp_data'];%b.nameData];
points     = 55;     % number of measured distances
b.z        = linspace(0,1.1,points);
b.usewin = 1;

for jj=1:points
	img = b.loadSmallImg(jj);       % loads small images created
    %img = img(54:520,54:520);
    [Fv Fh] = b.visCalc(img,1,jj);
    vert(jj) = Fv;
    horz(jj) = Fh;
end

fourCoeffexp = [horz.' vert.'];
% 
% fitpar = fitfunc(a,b,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP, ...
%                  src_DN,ene_UP,ene_DN,iterat,stepsz,range_sz,b.z, ...
%                  ['synth_' 'kev_simFIT']);
load('synth_kev_simFIT.mat')

fig1 = figure;
    set(fig1,'Position',[80 680 800 248]);
subplot(1,2,1)
    plot(horz./mean(horz))
    hold on
    plot(horzsim(:,end)./mean(horzsim(:,end)),'r')
    hold off
subplot(1,2,2)
    plot(vert./mean(vert))
    hold on
    plot(vertsim(:,end)./mean(vertsim(:,end)),'r')
    hold off