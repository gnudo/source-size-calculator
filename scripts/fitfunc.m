% function containing the actual fit algorithm
% TODO: - clean code
%       - unify nomenclature
%       - include into class
%       - rebase
%       - ...
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

a.n_max = n_max;
a.k_max = k_max;
a.s = s;

dutHDN  = dutDN;
dutVDN  = dutDN;
dutHUP  = dutUP;
dutVUP  = dutUP;
angHDN  = angDN;
angVDN  = angDN;
angHUP  = angUP;
angVUP  = angUP;
srcH_DN = src_DN;
srcV_DN = src_DN;
srcH_UP = src_UP;
srcV_UP = src_UP;

lsHorz_tmp = 1e11;
lsVert_tmp = 1e11;

for ii=1:a.k_max
[dutiesH d_dutiesH] = a.nestIntervals(dutHDN,dutHUP);
[dutiesV d_dutiesV] = a.nestIntervals(dutVDN,dutVUP);
[anglesH d_anglesH] = a.nestIntervals(angHDN,angHUP);
[anglesV d_anglesV] = a.nestIntervals(angVDN,angVUP);
[srcH d_srcH]       = a.nestIntervals(srcH_DN,srcH_UP);
[srcV d_srcV]       = a.nestIntervals(srcV_DN,srcV_UP);
[ene d_ene]         = a.nestIntervals(ene_DN,ene_UP);

for lene=1:a.n_max
 a.E = ene(lene);
 a.calcRefracAbsorb(17); % calculate delta & beta for grating with rho = 17
for lsrc=1:a.n_max
  srcsz = [srcH(lsrc) srcV(lsrc)];
  for lang=1:a.n_max
    angles = [anglesH(lang) anglesV(lang)];
	for ldut=1:a.n_max
      duties = [dutiesH(ldut) dutiesV(ldut)];
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
        [dutHDN dutHUP]   = a.setNewMinMax(duties(1),d_dutiesH);
        [angHDN angHUP]   = a.setNewMinMax(angles(1),d_anglesH);
        [srcH_DN srcH_UP] = a.setNewMinMax(srcsz(1),d_srcH);
        horzsim(:,countH) = fourCoeff(:,1);
        eneH(countH) = a.E;
      end % --> IF
      if lsVert < lsVert_tmp
        lsVert_tmp = lsVert;
        countV = countV+1;
        fitV(countV) = lsVert_tmp;
        [dutVDN dutVUP]   = a.setNewMinMax(duties(2),d_dutiesV);
        [angVDN angVUP]   = a.setNewMinMax(angles(2),d_anglesV);
        [srcV_DN srcV_UP] = a.setNewMinMax(srcsz(2),d_srcV);
        vertsim(:,countV) = fourCoeff(:,2);
        eneV(countV) = a.E;
      end % --> IF 
      if (lsHorz+lsVert) < sumVH
        countVH = countVH+1;
        m_ene(countVH) = a.E;
      end % --> IF 
    end % --> ldut
  end % --> lang
end % --> lsrc
end % --> ene
end % --> iter

% fitpar
% m_dutH
% m_angH
% m_srH
% m_dutV
% m_angV
% m_srV
% m_ene
m_dutH = (dutHDN + dutHUP) / 2;
m_angH = (angHDN + angHUP) / 2;
m_srH  = (srcH_DN + srcH_UP) / 2;
m_dutV = (dutVDN + dutVUP) / 2;
m_angV = (angVDN + angVUP) / 2;
m_srV  = (srcV_DN + srcV_UP) / 2;
m_ene = (ene_DN + ene_UP) /2;
d_dutH = d_dutiesH;
d_dutV = d_dutiesV;
d_angH = d_anglesH;
d_angV = d_anglesV;
d_srH = d_srcH;
d_srV = d_srcV;
%--------------------------------------------------------------------------
% 3.) Save to file
%--------------------------------------------------------------------------
save(name,'horzsim','vertsim','z','m_dutH','m_angH','m_srH','m_dutV', ...
          'm_angV','m_srV','d_dutH','d_dutV','d_angH','d_angV','d_srH', ...
          'd_srV','fitpar','fitH','fitV','m_ene','d_ene','eneH','eneV');
end