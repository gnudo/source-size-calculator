% Implementation of the fitting algorithm according to Fig. 3 from paper
% Lovric et al., Opt. Express 22, 2745 (2014): http://doi.org/q9m.
%--------------------------------------------------------------------------
% Date: 2014-01-24
% Author: Goran Lovric
% License: GPL 3 (see LICENSE file in root folder)
%--------------------------------------------------------------------------
function fit_algorithm(name,a,material,F_exp,dc_min,dc_max,alpha_min, ...
            alpha_max,sigma_min,sigma_max,E_min,E_max,k_max,n_max,s,resolu)
% number of iterations
N_max = 2 * n_max^4 * k_max * length(a.z)  % Eq. (18) from paper
counter_tot = n_max^4 *k_max;
counter = 0;  % help variable to see progress of fitting procedure

p_start_H  = inf;
p_start_V  = inf;
p_start_HV = inf;

% Create interval objects for horizontal an vertial directions
duty_H  = interval(dc_min,dc_max);
duty_V  = interval(dc_min,dc_max);
ang_H   = interval(alpha_min,alpha_max);
ang_V   = interval(alpha_min,alpha_max);
src_H   = interval(sigma_min,sigma_max);
src_V   = interval(sigma_min,sigma_max);
ene     = interval(E_min,E_max);

for k=1:k_max
% "Nest intervals" subroutine (for each fitting parameter)
dcH     = duty_H.nestIntervals(n_max,s);
dcV     = duty_V.nestIntervals(n_max,s);
alphaH  = ang_H.nestIntervals(n_max,s);
alphaV  = ang_V.nestIntervals(n_max,s);
sigmaH  = src_H.nestIntervals(n_max,s);
sigmaV  = src_V.nestIntervals(n_max,s);
E       = ene.nestIntervals(n_max,s);

for n_E=1:n_max
  a.E = E(n_E);
  a.calcRefracAbsorb(material,17);% calculate delta & beta
  for n_sigma=1:n_max
    sigma = [sigmaH(n_sigma) sigmaV(n_sigma)];
    for n_alpha=1:n_max
      alpha = [alphaH(n_alpha) alphaV(n_alpha)];
      for n_dc=1:n_max
        dc = [dcH(n_dc) dcV(n_dc)];
        counter = counter+1;
        disp([num2str(counter) ' of ' num2str(counter_tot)]);
      
        % "Fourier Analysis" subroutine for simulated data
        [F_sim_x,F_sim_y] = a.calculateFsim2D(a.E,sigma,dc,alpha,resolu);
        F_sim = [F_sim_x F_sim_y];
      
        % "Weighted LSE" subroutine
        p_H = a.weightedLSE(F_sim(:,1),F_exp(:,1)); % horizontal F-coeff
        p_V = a.weightedLSE(F_sim(:,2),F_exp(:,2)); % vertical F-coeff
        
        if p_H < p_start_H
          p_start_H = p_H;
          duty_H.setNewMinMax(dc(1));
          ang_H.setNewMinMax(alpha(1));
          src_H.setNewMinMax(sigma(1));
        end % --> IF
        if p_V < p_start_V
          p_start_V = p_V;
          duty_V.setNewMinMax(dc(2));
          ang_V.setNewMinMax(alpha(2));
          src_V.setNewMinMax(sigma(2));
        end % --> IF 
      end % --> n_dc
    end % --> n_alpha
  end % --> n_sigma
  % A new energy is only accepted if it minimizes p for both horizontal and
  % vertical Fourier coefficients (i.e.: we assume E must be the same for
  % both directions).
  if (p_H+p_V) < p_start_HV
    p_start_HV = p_H+p_V;
    ene.setNewMinMax(a.E);
  end % --> IF 
end % --> n_E
end % --> k
%--------------------------------------------------------------------------
% 3.) Save to file
%--------------------------------------------------------------------------
save(name,'a','F_exp','duty_H','duty_V','ang_H','ang_V','src_H', ...
     'src_V','ene','resolu');
end