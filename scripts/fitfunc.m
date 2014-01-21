% function containing the actual fit algorithm
% TODO: - clean code
%       - unify nomenclature
function fitpar = fitfunc(name,a,F_exp,dutUP,dutDN,angUP,angDN,src_UP,src_DN,ene_UP,ene_DN,k_max,n_max,s)
% number of iterations
N_max = 2 * n_max^4 * k_max * length(a.z)     % Eq. (18) from paper
COUNTER_NUM = n_max^4 *k_max

a.plotper = 14;      % number of periods to be plotted --> sets a.x1, a.x2
a.rr = a.r;

counter = 0;
countH = 0;
countV = 0;
countVH = 0;

% Create interval objects for horizontal an vertial directions
dc_H  = interval(dutDN,dutUP);
dc_V  = interval(dutDN,dutUP);
ang_H = interval(angDN,angUP);
ang_V = interval(angDN,angUP);
src_H = interval(src_DN,src_UP);
src_V = interval(src_DN,src_UP);
E     = interval(ene_DN,ene_UP);

p_start_H = 1e11;
p_start_V = 1e11;

for k=1:k_max
% Nest intervals for each parameter "object"
dutiesH = dc_H.nestIntervals(n_max,s);
dutiesV = dc_V.nestIntervals(n_max,s);
anglesH = ang_H.nestIntervals(n_max,s);
anglesV = ang_V.nestIntervals(n_max,s);
srcH    = src_H.nestIntervals(n_max,s);
srcV    = src_V.nestIntervals(n_max,s);
ene     = E.nestIntervals(n_max,s);

for n_E=1:n_max
 a.E = ene(n_E);
 a.calcRefracAbsorb(17); % calculate delta & beta for grating with rho = 17
for n_sigma=1:n_max
  sigma = [srcH(n_sigma) srcV(n_sigma)];
  for n_alpha=1:n_max
    angles = [anglesH(n_alpha) anglesV(n_alpha)];
	for n_dc=1:n_max
      dc = [dutiesH(n_dc) dutiesV(n_dc)];
      counter = counter+1
      
      % Loop through horizontal and vertical Fourier coefficients
      for ii = 1:2
          a.srcsz   = sigma(ii);
          a.duty    = dc(ii);
          a.gAngle  = angles(ii);
          F_sim(:,ii) = a.calcFcoeff;
      end
      
      % Least-Square-Fit
      p_H = a.modifiedLSE(F_sim(:,1),F_exp(:,1));
      p_V = a.modifiedLSE(F_sim(:,2),F_exp(:,2));
      fitpar(counter,:) = [p_H p_V];

      p_start_HV = p_start_H + p_start_V;
      if p_H < p_start_H
        p_start_H = p_H;
        countH = countH+1;
        fitH(countH) = p_start_H;
        dc_H.setNewMinMax(dc(1));
        ang_H.setNewMinMax(angles(1));
        src_H.setNewMinMax(sigma(1));
        horzsim(:,countH) = F_sim(:,1);
        eneH(countH) = a.E;
      end % --> IF
      if p_V < p_start_V
        p_start_V = p_V;
        countV = countV+1;
        fitV(countV) = p_start_V;
        dc_V.setNewMinMax(dc(2));
        ang_V.setNewMinMax(angles(2));
        src_V.setNewMinMax(sigma(2));
        vertsim(:,countV) = F_sim(:,2);
        eneV(countV) = a.E;
      end % --> IF 
      if (p_H+p_V) < p_start_HV
        countVH = countVH+1;
        m_ene(countVH) = a.E;
        E.setNewMinMax(a.E);
      end % --> IF 
    end % --> n_dc
  end % --> n_alpha
end % --> n_sigma
end % --> n_E
end % --> k

m_dutH = (dc_H.DN +dc_H.UP)/2; %(dutHDN + dutHUP) / 2;
m_angH = (ang_H.DN + ang_H.UP)/2; %(angHDN + angHUP) / 2;
m_srH  = (src_H.DN + src_H.UP)/2; %(srcH_DN + srcH_UP) / 2;
m_dutV = (dc_V.DN +dc_V.UP)/2; %(dutVDN + dutVUP) / 2;
m_angV = (ang_V.DN + ang_V.UP)/2; %(angVDN + angVUP) / 2;
m_srV  = (src_V.DN + src_V.UP)/2; %(srcV_DN + srcV_UP) / 2;
m_ene = (E.DN + E.UP)/2; %(ene_DN + ene_UP) /2;
d_dutH = dc_H.del;%d_dutiesH;
d_dutV = dc_V.del;%d_dutiesV;
d_angH = ang_H.del;%d_anglesH;
d_angV = ang_V.del;%d_anglesV;
d_srH = src_H.del;%d_srcH;
d_srV = src_V.del;%d_srcV;
d_ene = E.del;
z = a.z;
%--------------------------------------------------------------------------
% 3.) Save to file
%--------------------------------------------------------------------------
save(name,'horzsim','vertsim','z','m_dutH','m_angH','m_srH','m_dutV', ...
          'm_angV','m_srV','d_dutH','d_dutV','d_angH','d_angV','d_srH', ...
          'd_srV','fitpar','fitH','fitV','m_ene','d_ene','eneH','eneV');
end