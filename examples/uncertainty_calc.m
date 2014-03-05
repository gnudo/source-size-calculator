% This is a wrapper function for calculating the energy and source size
% uncertainties with the Newton method. The fitting function has to be run
% first and the results have to be saved into <results.mat>.
%--------------------------------------------------------------------------
% Date: 2014-02-28
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
addpath('../functions');

load('results.mat');

del_z   = 0.002;                % [m] propagation distance uncertainty

% starting values for Newton searches
E_start = ene.val.*1.05;        % [keV] starting value of E
sigma_H_start = src_H.val*1.5;  % [m] starting value of simga_H
sigma_V_start = src_V.val*1.5;  % [m] starting value of simga_V

%--------------------------------------------------------------------------
% 1.) Calculate physical uncertainties
%--------------------------------------------------------------------------
[del_E del_sigma]= errorcalc(a,del_z,ene.val,E_start, ...
                            [src_H.val src_V.val], ...
                            [sigma_H_start sigma_V_start], ...
                            [duty_H.val duty_V.val],[ang_H.val ang_V.val]);

%--------------------------------------------------------------------------
% 2.) Add uncertainties from fitting procedure
%--------------------------------------------------------------------------                        
del_E_eff = del_E + ene.del;
del_sigma_eff(1) = del_sigma(1) + src_H.del;
del_sigma_eff(2) = del_sigma(2) + src_V.del;

%--------------------------------------------------------------------------
% 3.) Save all into results.mat
%--------------------------------------------------------------------------
save('results.mat','del_E_eff','del_sigma_eff','-append')