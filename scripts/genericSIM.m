% generic simulation script for calculation under given parameters
% TODO: like in all others (probably script will be obsolete)
%--------------------------------------------------------------------------
% Author: Goran Lovric
%--------------------------------------------------------------------------
function [horzsim,vertsim] = genericSIM(name,radius,energy,z,srcsz,duties,angles)
%addpath('../classes');
a = simulation;               % load classes for calculation
%--------------------------------------------------------------------------
% 1.) PATHS (m-files & src-files) and CONSTANTS
%--------------------------------------------------------------------------
a.a       = 6.84e-6;               % [m] grating period
a.r       = radius;                % [m] radius of incidient wave curvature
a.rr = radius;%4.6;
a.E       = energy;                % [keV]
a.periods = 32;                    % grating-size (in terms of periods)
a.N       = 2^13;                  % number of particles --> 2^n !!!
a.padding = 1*a.N;                 % total number (with zero-padding)

% Grating parameters
a.gHeight = 3.39e-6;               % height of grating structure  !!!! EDIT

for kk=1:2
a.srcsz   = srcsz(kk);             % TOMCAT horizontal src size
a.duty    = duties(kk);
a.gAngle  = angles(kk);

% GRID creation & Wave field @ Grid
gra = a.talbotGrid1D;
f = a.waveFieldGrat(gra);

% Propagation along z-axis
a.plotper = 14;      % number of periods to be plotted --> sets a.x1, a.x2
resolu    = 512;     % desired resolution

for ii=1:length(z)
    calc = a.waveFieldPropMutual(z(ii),f); % Fresnel-pro + mutual coh.fnct
    crop = calc(a.x1:a.x2);                % crop to <plotper> periods
    crop = interp1(crop,linspace(1,length(crop),resolu));
    pwav(:,ii) = crop;
end

%--------------------------------------------------------------------------
% 2.) Visibility calcs
%--------------------------------------------------------------------------
b   = analysis;
M   = ( a.r + z) ./ a.r;
per = M .* size(pwav,1)/a.plotper;        % period in [px]

for ii=1:length(z)
    vec = pwav(:,ii);
    fourCoeff(ii,kk) = b.vis1D (vec,per(ii));
end
end
%--------------------------------------------------------------------------
% 3.) Output coefficients
%--------------------------------------------------------------------------
horzsim = fourCoeff(:,1);
vertsim = fourCoeff(:,2);
end