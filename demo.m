% Run full test to demonstrate algorithm chain for source size calculation.
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;

%--------------------------------------------------------------------------
% 0.) create "fake" experimental data (takes ~5min)
%--------------------------------------------------------------------------
%run examples/synthetic_data

%--------------------------------------------------------------------------
% 1.) Apply the fitting algorithm from G. Lovric et. al. (submitted)
%--------------------------------------------------------------------------
if ~exist('examples/results.mat', 'file')
    run examples/fitting
end

%--------------------------------------------------------------------------
% 2.) Calculate uncertainties
%--------------------------------------------------------------------------
run examples/uncertainty_calc

%--------------------------------------------------------------------------
% 3.) Print all results
%--------------------------------------------------------------------------
run examples/print_all