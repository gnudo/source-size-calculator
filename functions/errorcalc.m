% Function for calculating the energy and source size uncertainties
% according to Eqs. (19) and (20) from Lovric et al., Opt. Express 22,
% 2745 (2014). For solving the equations, we make use of the Newton method
% (http://en.wikipedia.org/wiki/Newton_method).
%--------------------------------------------------------------------------
% Date: 2014-02-28
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
function [del_E,del_sigma] = errorcalc(a,del_z,E,E_start,sigma, ...
                                                      sigma_start,dc,alpha)
% (0) Constants that are used below
eps = 0.01;  % precision of Newton method
f_H = 1;     % difference between two p after consecutive Newton iterations
f_V = f_H;
ff  = f_H;

% (1) First we calculate the left side of Eq. (19)
[F_simH,F_simV]   = a.calcFcoeff2D(E,[sigma(1) sigma(2)],...
                            [dc(1) dc(2)],[alpha(1) alpha(2)]);

a.z = a.z + del_z;
[F_simH2,F_simV2] = a.calcFcoeff2D(E,[sigma(1) sigma(2)],...
                            [dc(1) dc(2)],[alpha(1) alpha(2)]);
                        
p_H = a.weightedLSE(F_simH,F_simH2);
p_V = a.weightedLSE(F_simV,F_simV2);

% (2) Now we search for del_sigma of the right side from Eq. (19)
a.z = a.z - del_z;
sigma_H_start = sigma_start(1);
sigma_V_start = sigma_start(2);

while f_H > eps || f_V > eps
d_sigma_H = (sigma_H_start - sigma(1)) .* [1 1.05];
d_sigma_V = (sigma_V_start - sigma(2)) .* [1 1.05];

for ii=1:2
    [F_simH2,F_simV2] = a.calcFcoeff2D(E,[sigma(1)+d_sigma_H(ii) ...
                                       sigma(2)+d_sigma_V(ii)], ...
                            [dc(1) dc(2)],[alpha(1) alpha(2)]);
	del_p_H(ii) = a.weightedLSE(F_simH2,F_simH) - p_H;
    del_p_V(ii) = a.weightedLSE(F_simV2,F_simV) - p_V;
end
f_H = del_p_H(1);
f_V = del_p_V(1);

f_H_prime = (del_p_H(2)-del_p_H(1))/(d_sigma_H(2)-d_sigma_H(1));
f_V_prime = (del_p_V(2)-del_p_V(1))/(d_sigma_V(2)-d_sigma_V(1));
sigma_H_start = sigma_H_start - del_p_H(1) ./ f_H_prime; % Newton iteration
sigma_V_start = sigma_V_start - del_p_V(1) ./ f_V_prime; % Newton iteration
end

del_sigma = [d_sigma_H(1) d_sigma_V(1)];

% (3) We repeat the same for Eq. (20)
while ff > eps
d_E = (E_start - E) .* [1 1.05];

for ii=1:2
    [F_simH2,F_simV2] = a.calcFcoeff2D(E+d_E(ii),[sigma(1) sigma(2)], ...
                            [dc(1) dc(2)],[alpha(1) alpha(2)]);
	del_p(ii) = a.weightedLSE(F_simH2,F_simH) - max([p_H p_V]);
end
ff = del_p(1);

ff_prime = (del_p(2)-del_p(1))/(d_E(2)-d_E(1));
E_start = E_start - del_p(1) ./ ff_prime;                % Newton iteration
end

del_E = d_E(1)*3;  % according to Eq. (26)
end