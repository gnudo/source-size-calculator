% Performs the fitting algorithm on "experimental" data, being located in
% "b.nameDest" folder. To apply on arbitrary experimental data, adjust
% section 1.) with correct parameters
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
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

a.gHeight  = 3.39e-6;             % [m] height of grating structure %!!!!!!!!!!!!!!!!!!!!
a.periods  = 32;                  % grating-size (in terms of periods)
a.N        = 2^13;                % number of particles (pixels)

dutDN  = 0.5;                     % [1] duty cycle lower limit
dutUP  = 0.54;                    % [1] duty cycle upper limit
angDN  = 0;                       % [°] grating angle lower limit
angUP  = 4;                       % [°] grating angle upper limit
src_DN = 0e-6;                    % [m] source size lower limit
src_UP = 200e-6;                  % [m] source size upper limit
ene_DN = 14-0.5;                  % [keV] x-ray energy lower limit
ene_UP = 14+0.5;                  % [keV] x-ray energy upper limit
s      = 1.5;                     % interval stretching factor (from paper)

n_max  = 3;                       % number of intervals to be nested
k_max  = 7;                       % number of iteration steps

%--------------------------------------------------------------------------
% 2.) Calculate Fourier coefficients from experimental data
%--------------------------------------------------------------------------
for jj=1:length(a.z)
	img     = a.loadSmallImg(jj); % load experimental Talbot imgs
    [F_Hi F_Vi] = a.visCalc(img,jj);  % extract first Fourier coefficients
    F_exp(jj,1) = F_Hi;
    F_exp(jj,2) = F_Vi;
end

%--------------------------------------------------------------------------
% 3.) Run fitting algorithm with the above parameters (or load results)
%--------------------------------------------------------------------------
if ~exist('results.mat', 'file')
    fitpar = fitfunc('results', a,F_exp,dutUP,dutDN,angUP,angDN, ...
                     src_UP,src_DN,ene_UP,ene_DN,k_max,n_max,s);
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
    plot(a.z,F_exp(:,1)./mean(F_exp(:,1)),'o')
    hold on;
    plot(a.z,horzsim(:,end)./mean(horzsim(:,end)),'k')
    hold off;
    xlim([0 a.z(end)]);
    legend('experimental values','best fit')
subplot(1,2,2)
    plot(a.z,F_exp(:,2)./mean(F_exp(:,2)),'ro')
    hold on;
    plot(a.z,vertsim(:,end)./mean(vertsim(:,end)),'k')
    hold off;
    xlim([0 a.z(end)]);
    legend('experimental values','best fit')