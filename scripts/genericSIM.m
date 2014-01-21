% generic simulation function for calculating Fourier coefficients under
% given parameters.
% --> almost ready for deletion
%--------------------------------------------------------------------------
% Author: Goran Lovric
%--------------------------------------------------------------------------
function [horzsim,vertsim] = genericSIM(z,E,a,h,alpha,d_c, sigma,R,RR)
addpath('../classes');
addpath('../classes/xray-interaction-constants');
calc = simulation;               % load classes for calculation
%--------------------------------------------------------------------------
% 1.) PATHS (m-files & src-files) and CONSTANTS
%--------------------------------------------------------------------------
calc.a       = a;               % [m] grating period
calc.r       = R;                % [m] radius of incidient wave curvature
calc.rr      = RR;
calc.E       = E;                % [keV]
calc.periods = 32;                    % grating-size (in terms of periods)
calc.N       = 2^13;                  % number of particles --> 2^n !!!

% Grating parameters
calc.gHeight = h;               % height of grating structure  !!!! EDIT

for kk=1:2
calc.srcsz   = sigma(kk);             % TOMCAT horizontal src size
calc.duty    = d_c(kk);
calc.gAngle  = alpha(kk);

% GRID creation & Wave field @ Grid
calc.calcRefracAbsorb(17);
gra = calc.talbotGrid1D;
f = calc.waveFieldGrat(gra);

% Propagation along z-axis
calc.plotper = 14;      % number of periods to be plotted --> sets calc.x1, calc.x2
resolu    = 512;     % desired resolution

for ii=1:length(z)
    lineTalbot = calc.waveFieldPropMutual(z(ii),f); % Fresnel-pro + mutual coh.fnct
    crop = lineTalbot(calc.x1:calc.x2);                % crop to <plotper> periods
    crop = interp1(crop,linspace(1,length(crop),resolu));
    pwav(:,ii) = crop;
end

%--------------------------------------------------------------------------
% 2.) Visibility calcs
%--------------------------------------------------------------------------
M   = ( calc.r + z) ./ calc.r;
per = M .* size(pwav,1)/calc.plotper;        % period in [px]

for ii=1:length(z)
    vec = pwav(:,ii);
    fourCoeff(ii,kk) = calc.vis1D (vec,per(ii));
end
end
%--------------------------------------------------------------------------
% 3.) Output coefficients
%--------------------------------------------------------------------------
horzsim = fourCoeff(:,1);
vertsim = fourCoeff(:,2);
end