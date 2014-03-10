% Print and plot results for source sizes and uncertainties
%--------------------------------------------------------------------------
% Date: 2014-02-28
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
close all;clc;clear;

load('results.mat')

%--------------------------------------------------------------------------
% 1.) Simulate F-coefficients from the loaded parameters (for plotting)
%--------------------------------------------------------------------------
[F_simH F_simV] = a.calculateFsim2D(ene.val,[src_H.val src_V.val],...
                            [duty_H.val duty_V.val],[ang_H.val ang_V.val]);

%--------------------------------------------------------------------------
% 2.) Print and plot the results
%--------------------------------------------------------------------------

% Source sizes
clc;
disp(['Hor. source size = (' num2str(round(src_H.val*1e6)), ...
    ' +- ' num2str(round(del_sigma_eff(1)*1e6)) ') micrometer']);
disp(['Ver. source size = (' num2str(round(src_V.val*1e6)), ...
    ' +- ' num2str(round(del_sigma_eff(2)*1e6)) ') micrometer']);
fprintf('\n');
% Duty cycles
disp(['Hor. duty cycle = ' num2str(duty_H.val)]);
disp(['Ver. duty cycle = ' num2str(duty_V.val)]);
fprintf('\n');
% Angle
disp(['Hor. angle = ' num2str(ang_H.val)]);
disp(['Ver. angle = ' num2str(ang_V.val)]);
fprintf('\n');
% Energy
disp(['X-ray energy = ' num2str(ene.val) ' +- ' num2str(del_E_eff) ' keV']);

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