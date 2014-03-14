% Run full test to demonstrate algorithm chain for source size calculation
% from Lovric et al., Opt. Express 22, 2745 (2014): http://doi.org/q9m.
%--------------------------------------------------------------------------
% Date: 2013-12-20
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;

%--------------------------------------------------------------------------
% 0.) Create "fake" experimental data (takes ~5min)
%--------------------------------------------------------------------------
run examples/create_fake_data

%--------------------------------------------------------------------------
% 1.) Run fitting algorithm
%--------------------------------------------------------------------------
if ~exist('examples/results.mat', 'file')
    run examples/run_fit
end

%--------------------------------------------------------------------------
% 2.) Calculate uncertainties
%--------------------------------------------------------------------------
run examples/run_uncertainty_calc

%--------------------------------------------------------------------------
% 3.) Print all results
%--------------------------------------------------------------------------
run examples/show_results