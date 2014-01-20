% function containing the actual fit algorithm
% TODO: - clean code
%       - unify nomenclature
%       - rebase
%       - use generic simulation+FourierCalc script (to be created)
function fitpar = fitfunc(a,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP,src_DN,ene_UP,ene_DN,k_max,n_max,s,z,name)
% number of iterations
N_max = 2 * n_max^4 * k_max * length(z)     % Eq. (18) from paper
COUNTER_NUM = n_max^4 *k_max

a.plotper = 14;      % number of periods to be plotted --> sets a.x1, a.x2
resolu    = 512;     % desired resolution
M   = ( a.r + z) ./ a.r;

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

lsHorz_tmp = 1e11;
lsVert_tmp = 1e11;

for k=1:k_max
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
  srcsz = [srcH(n_sigma) srcV(n_sigma)];
  for n_alpha=1:n_max
    angles = [anglesH(n_alpha) anglesV(n_alpha)];
	for n_dc=1:n_max
      duties = [dutiesH(n_dc) dutiesV(n_dc)];
      counter = counter+1
      %--------------------------------------------------------------------
      for kk=1:2
        a.gAngle  = angles(kk);
        a.srcsz   = srcsz(kk);                 % TOMCAT horizontal src size
        a.duty    = duties(kk);
        
        gra = a.talbotGrid1D;
        f = a.waveFieldGrat(gra);
        for ii=1:length(z)
          calc = a.waveFieldPropMutual(z(ii),f);
          crop = calc(a.x1:a.x2);
          crop = interp1(crop,linspace(1,length(crop),resolu));
          pwav(:,ii) = crop;
          % --- visibility calcs
          per = M .* size(pwav,1)/a.plotper;        % period in [px]
          vec = pwav(:,ii);
          fourCoeff(ii,kk) = a.vis1D (vec,per(ii));
        end
      end

      % Least-Square-Fit
      lsHorz = a.weightedLSQ(fourCoeff(:,1),fourCoeffexp(:,1));
      lsVert = a.weightedLSQ(fourCoeff(:,2),fourCoeffexp(:,2));
      fitpar(counter,:) = [lsHorz lsVert];

      sumVH = lsHorz_tmp + lsVert_tmp;
      if lsHorz < lsHorz_tmp
        lsHorz_tmp = lsHorz;
        countH = countH+1;
        fitH(countH) = lsHorz_tmp;
        dc_H.setNewMinMax(duties(1));
        ang_H.setNewMinMax(angles(1));
        src_H.setNewMinMax(srcsz(1));
        horzsim(:,countH) = fourCoeff(:,1);
        eneH(countH) = a.E;
      end % --> IF
      if lsVert < lsVert_tmp
        lsVert_tmp = lsVert;
        countV = countV+1;
        fitV(countV) = lsVert_tmp;
        dc_V.setNewMinMax(duties(2));
        ang_V.setNewMinMax(angles(2));
        src_V.setNewMinMax(srcsz(2));
        vertsim(:,countV) = fourCoeff(:,2);
        eneV(countV) = a.E;
      end % --> IF 
      if (lsHorz+lsVert) < sumVH
        countVH = countVH+1;
        m_ene(countVH) = a.E;
        E.setNewMinMax(a.E);
      end % --> IF 
    end % --> n_dc
  end % --> n_alpha
end % --> n_sigma
end % --> n_E
end % --> iter

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
%--------------------------------------------------------------------------
% 3.) Save to file
%--------------------------------------------------------------------------
save(name,'horzsim','vertsim','z','m_dutH','m_angH','m_srH','m_dutV', ...
          'm_angV','m_srV','d_dutH','d_dutV','d_angH','d_angV','d_srH', ...
          'd_srV','fitpar','fitH','fitV','m_ene','d_ene','eneH','eneV');
end