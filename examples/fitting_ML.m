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
addpath('../functions');
a = simulation;                   % load methods for calculation

%--------------------------------------------------------------------------
% 1.) Experimental parameters + fitting parameter margins
%--------------------------------------------------------------------------
a.nameDest = 'exp_data_refl_beam';% folder that contains experimental data
a.psize    = 0.76e-6;             % [m] px size of detector
a.R        = 25;                  % [m] source-to-grating distance
a.a        = 6.84e-6;             % [m] grating period
a.z        = linspace(0,1.1,56);  % [m] experimental propagation distances
a.usewin   = 1;                   % apply Tukey/Hanning window function
a.E        = 18;                  % [keV] X-ray energy

a.h        = 3.39e-6;             % [m] height of grating structure
a.periods  = 32;                  % grating-size (in terms of periods)
a.N        = 2^13;                % number of particles (pixels)

dc_min     = 0.5;                 % [1] duty cycle lower limit
dc_max     = 0.54;                % [1] duty cycle upper limit
alpha_min  = 0;                   % [°] grating angle lower limit
alpha_max  = 4;                   % [°] grating angle upper limit
sigma_min  = 0e-6;                % [m] source size lower limit
sigma_max  = 200e-6;              % [m] source size upper limit
R_min      = 1;                   % [m] curvature radius lower limit
R_max      = 26;                  % [m] curvature radius upper limit
s          = 1.5;                 % interval stretching factor (from paper)

n_max      = 3;                   % number of intervals to be nested
k_max      = 7;                   % number of iteration steps
resolu     = 512;                 % resolution for Talbot images simulation

%--------------------------------------------------------------------------
% 2.) Calculate Fourier coefficients from experimental data
%--------------------------------------------------------------------------
for jj=1:length(a.z)
	img         = a.loadSmallImg(jj);       % load experimental Talbot imgs
    [F_Hi F_Vi] = a.FourierAnalysis2D(img,jj); % extract 1st F-coefficients
    F_exp(jj,1) = F_Hi;
    F_exp(jj,2) = F_Vi;
end

%--------------------------------------------------------------------------
% 3.) Run fitting algorithm with the above parameters (or load results)
%--------------------------------------------------------------------------
if ~exist('results_ML.mat', 'file')
    fit_algorithm_ML('results_ML',a,'Au',F_exp,dc_min,dc_max,alpha_min, ...
           alpha_max,sigma_min,sigma_max,R_min,R_max,k_max,n_max,s,resolu);
end
load('results_ML.mat')

%--------------------------------------------------------------------------
% 4.) Simulate F-coefficients from the loaded parameters (for plotting)
%--------------------------------------------------------------------------
a.sigma = src_H.val;
a.dc    = duty_H.val;
a.alpha = ang_H.val;
a.RR    = rad_H.val;
F_simH  = a.calculateFsim(resolu);

a.sigma = src_V.val;
a.dc    = duty_V.val;
a.alpha = ang_V.val;
a.RR    = rad_V.val;
F_simV  = a.calculateFsim(resolu);

%--------------------------------------------------------------------------
% 5.) Print and plot the results
%--------------------------------------------------------------------------

% Source sizes
clc;
disp(['Hor. source size = (' num2str(round(src_H.val*1e6)), ...
    ' +- ' num2str(round(src_H.del*1e6)) ') micrometer']);
disp(['Ver. source size = (' num2str(round(src_V.val*1e6)), ...
    ' +- ' num2str(round(src_V.del*1e6)) ') micrometer']);
fprintf('\n');
% Duty cycles
disp(['Hor. duty cycle = ' num2str(duty_H.val)]);
disp(['Ver. duty cycle = ' num2str(duty_V.val)]);
fprintf('\n');
% Angle
disp(['Hor. angle = ' num2str(ang_H.val)]);
disp(['Ver. angle = ' num2str(ang_V.val)]);
fprintf('\n');
% Source-to sample distance
disp(['R = (' num2str(rad_H.val) ' +- ' num2str(rad_H.del) ') m']);
disp(['R = (' num2str(rad_V.val) ' +- ' num2str(rad_V.del) ') m']);

disp(['Note that these are not physical uncertainties, but just from', ...
     ' the fitting procedure!'])
% Plot Fourier coefficients
fig1 = figure;
    set(fig1,'Position',[80 680 800 248]);
subplot(1,2,1)
    plot(a.z,F_exp(:,1)./mean(F_exp(:,1)),'o')
    hold on;
    plot(a.z,F_simH./mean(F_simH),'k')
    hold off;
    xlim([0 a.z(end)]);
    legend('experimental values','best fit')
subplot(1,2,2)
    plot(a.z,F_exp(:,2)./mean(F_exp(:,2)),'ro')
    hold on;
    plot(a.z,F_simV./mean(F_simV),'k')
    hold off;
    xlim([0 a.z(end)]);
    legend('experimental values','best fit')