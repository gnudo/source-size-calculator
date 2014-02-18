% simple calculation of Talbot-carpet for "paralell" (set r > 1e6) and
% "cone-beam" impinging wave-front.
%--------------------------------------------------------------------------
% Name: Talbot-simulation
% Date: 2013-12-13
% Author: Goran Lovric
%--------------------------------------------------------------------------
close all;clc;clear;
addpath('../classes');
addpath('../classes/xray-interaction-constants');
a = simulation;               % load classes for calculation
%--------------------------------------------------------------------------
% 1.) PATHS (m-files & src-files) and CONSTANTS
%--------------------------------------------------------------------------
a.a       = 6.84e-6;               % [m] grating period
a.r       = 25;%1e11;              % [m] radius of incidient wave curvature
a.E       = 21;                    % [keV]
a.srcsz   = 53e-6;                 % TOMCAT horizontal src size

a.periods = 8;                     % grating-size (in terms of periods)
a.N       = 2^13;                  % number of particles --> 2^n !!!
a.padding = 1*a.N;                 % total number (with zero-padding)
a.psize   = 0.38e-6;               % [m] px size of detector

% Grating parameters
a.gHeight = 3.39e-6;                % height of grating structure
a.gAngle  = 4.2;                    % angle of bump's slope
a.duty    = sqrt(0.285);            % duty cycle

a.z = linspace(0,a.D_def*2,501);  % [m] propagation distance in meters

%--------------------------------------------------------------------------
% 2.) GRID creation
%--------------------------------------------------------------------------
gra = a.talbotGrid1D;
%--------------------------------------------------------------------------
% 3.) Wave field @ Grid
%--------------------------------------------------------------------------
a.calcRefracAbsorb(17);
f = a.waveFieldGrat(gra);
a.plotWaveAtGrating(f,gra);

%--------------------------------------------------------------------------
% 4.) Propagation along z-axis
%--------------------------------------------------------------------------
a.plotper = 4;      % number of periods to be plotted --> sets a.x1, a.x2

for ii=1:length(a.z)
    ii
    calc = a.waveFieldPropMutual(a.z(ii),f);
    crop = calc(a.x1:a.x2);     % crop to <plotper> periods
    crop = a.scale2Det(crop);   % scale to detector resolution
    pwav(:,ii) = crop;
    %--- Plot ---
    if (mod(ii,100) == 0)
        fig2=figure(2);
        set(fig2,'Position',[100 100 1024 400],'Color','white')
        sb1 = subplot (1,2,1);
        plot(calc);
        title('vertical mean')
    end
end

sb2 = subplot (1,2,2);
figure(2)
imagesc(a.z*1000,[],pwav);
colormap(gray);
line(a.D_R.*1e3.*[1 1],[0 a.padding],'Color','r','LineWidth',2);
line(a.D_T.*1e3.*[1 1],[0 a.padding],'Color','r','LineWidth',2);
line(a.D_defr.*1e3.*[1 1],[0 a.padding],'Color','b','LineWidth',2);
line(a.D_def.*1e3.*[1 1],[0 a.padding],'Color','b','LineWidth',2);
xlabel('propagation distance [mm]');
legend('D_R','D_T','D_R (defocussed)','D_T (defocussed)');