% Performs the fitting algorithm on "experimental" data, being located in
% "b.nameDest" folder. To apply on arbitrary experimental data, adjust
% section 1.) with correct parameters
% TODO: - reduce redundancies and unite "simulation" and "analysis" classes
%       - move object creation and class properties settings to fitfunc
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
addpath('../scripts');

a = simulation;                   % load methods for calculation
b = analysis;                     % load method for saving images

%--------------------------------------------------------------------------
% 1.) Experimental parameters + fitting parameter margins
%--------------------------------------------------------------------------
b.nameDest = 'exp_data';          % folder that contains experimental data
b.psize    = 0.38e-6;             % [m] px size of detector
b.r        = 25;                  % [m] source-to-grating distance
b.a        = 6.84e-6;             % [m] grating period
b.z        = linspace(0,1.1,55);  % [m] experimental propagation distances
b.usewin   = 1;                   % apply Tukey/Hanning window function

a.E        = 14;                  % [keV] x-ray energy
a.r        = b.r;                 % [m] radius of incidient wave curvature
a.a        = b.a;                 % [m] grating period
a.gHeight  = 3.39e-6;             % [m] height of grating structure
a.periods  = 32;                  % grating-size (in terms of periods)
a.N        = 2^13;                % number of particles (pixels)

dutDN  = 0.5;                     % [1] duty cycle lower limit
dutUP  = 0.54;                    % [1] duty cycle upper limit
angDN  = 0;                       % [°] grating angle lower limit
angUP  = 4;                       % [°] grating angle upper limit
src_DN = 0e-6;                    % [m] source size lower limit
src_UP = 200e-6;                  % [m] source size upper limit
ene_DN = a.E-0.5;                 % [keV] x-ray energy lower limit
ene_UP = a.E+0.5;                 % [keV] x-ray energy upper limit
s      = 1.5;                     % interval stretching factor (from paper)

n_max  = 3;                       % number of intervals to be nested
k_max  = 7;                       % number of iteration steps

%--------------------------------------------------------------------------
% 2.) Calculate Fourier coefficients from experimental data
%--------------------------------------------------------------------------
for jj=1:length(b.z)
	img     = b.loadSmallImg(jj); % load experimental Talbot imgs
    [Fv Fh] = b.visCalc(img,jj);  % extract first Fourier coefficients
    Fvv(jj) = Fv;
    Fhh(jj) = Fh;
end
fourCoeffexp = [Fhh.' Fvv.'];

%--------------------------------------------------------------------------
% 3.) Run fitting algorithm with the above parameters (or load results)
%--------------------------------------------------------------------------
if ~exist('results.mat', 'file')
    fitpar = fitfunc(a,b,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP, ...
                     src_DN,ene_UP,ene_DN,k_max,n_max,s,b.z, 'results');
end
load('results.mat')

%--------------------------------------------------------------------------
% 4.) Print and plot the results
%--------------------------------------------------------------------------

% Source sizes
clc;
disp(['Hor. source size = ' num2str(m_srH(end)*1e6) ' micrometer']);
disp(['Ver. source size = ' num2str(m_srV(end)*1e6) ' micrometer']);
fprintf('\n');
% Duty cycles
disp(['Hor. duty cycle = ' num2str(m_dutH(end))]);
disp(['Ver. duty cycle = ' num2str(m_dutV(end))]);
fprintf('\n');
% Angle
disp(['Hor. angle = ' num2str(m_angH(end))]);
disp(['Ver. angle = ' num2str(m_angV(end))]);

% Plot Fourier coefficients
fig1 = figure;
    set(fig1,'Position',[80 680 800 248]);
subplot(1,2,1)
    plot(b.z,Fhh./mean(Fhh),'o')
    hold on;
    plot(b.z,horzsim(:,end)./mean(horzsim(:,end)),'k')
    hold off;
    xlim([0 b.z(end)]);
    legend('experimental values','best fit')
subplot(1,2,2)
    plot(b.z,Fvv./mean(Fvv),'ro')
    hold on;
    plot(b.z,vertsim(:,end)./mean(vertsim(:,end)),'k')
    hold off;
    xlim([0 b.z(end)]);
    legend('experimental values','best fit')