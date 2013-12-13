% function containing the actual fit algorithm
% TODO: - clean code
%       - unify nomenclature
%       - include into class
%       - rebase
%       - ...
function fitpar = fitfunc(a,b,fourCoeffexp,dutUP,dutDN,angUP,angDN,src_UP,src_DN,ene_UP,ene_DN,iterat,stepsz,range_sz,z,name)
% number of iterations
NUM_ITERATIONS = 2 * stepsz^4 * iterat * length(z)
COUNTER_NUM = stepsz^4 *iterat

a.plotper = 14;      % number of periods to be plotted --> sets a.x1, a.x2
resolu    = 512;     % desired resolution
M   = ( a.r + z) ./ a.r;

counter = 0;
countH = 0;
countV = 0;
countVH = 0;

d_dutH = (dutUP-dutDN)./2;
m_dutH = (dutUP+dutDN)./2;
d_dutV = d_dutH;
m_dutV = m_dutH;
d_angH = (angUP-angDN)./2;
m_angH = (angUP+angDN)./2;
d_angV = d_angH;
m_angV = m_angH;
d_srH = (src_UP-src_DN)./2;
m_srH = (src_UP+src_DN)./2;
d_srV = d_srH;
m_srV = m_srH;
d_ene = (ene_UP-ene_DN)./2;
m_ene = (ene_UP+ene_DN)./2;
lsHorz_tmp = 1e11;
lsVert_tmp = 1e11;

for ii=1:iterat
dutiesH = linspace(m_dutH(end)-d_dutH,m_dutH(end)+d_dutH,stepsz+2);
dutiesV = linspace(m_dutV(end)-d_dutV,m_dutV(end)+d_dutV,stepsz+2);
anglesH = linspace(m_angH(end)-d_angH,m_angH(end)+d_angH,stepsz+2);
anglesV = linspace(m_angV(end)-d_angV,m_angV(end)+d_angV,stepsz+2);
srcH   = linspace(m_srH(end)-d_srH,m_srH(end)+d_srH,stepsz+2);
srcV   = linspace(m_srV(end)-d_srV,m_srV(end)+d_srV,stepsz+2);
ene     = linspace(m_ene(end)-d_ene,m_ene(end)+d_ene,stepsz+2);

% divide range by 3 for next round
d_dutH = (d_dutH/stepsz).*range_sz;
d_dutV = (d_dutV/stepsz).*range_sz;
d_angH = (d_angH/stepsz).*range_sz;
d_angV = (d_angV/stepsz).*range_sz;
d_srH  = (d_srH/stepsz).*range_sz;
d_srV  = (d_srV/stepsz).*range_sz;
d_ene  = (d_ene/stepsz).*range_sz;

for lene=2:stepsz+1
 a.E = ene(lene);
 a.gHeight = 3.39e-6;
for lsrc=2:stepsz+1
  srcsz = [srcH(lsrc) srcV(lsrc)];
  for lang=2:stepsz+1
    angles = [anglesH(lang) anglesV(lang)];
    for ldut=2:stepsz+1
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
          fourCoeff(ii,kk) = b.vis1D (vec,per(ii));
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
        m_dutH(countH) = duties(1);
        m_angH(countH) = angles(1);
        m_srH(countH)  = srcsz(1);
        horzsim(:,countH) = fourCoeff(:,1);
        eneH(countH) = a.E;
      end % --> IF
      if lsVert < lsVert_tmp
        lsVert_tmp = lsVert;
        countV = countV+1;
        fitV(countV) = lsVert_tmp;
        m_dutV(countV) = duties(2);
        m_angV(countV) = angles(2);
        m_srV(countV)  = srcsz(2);
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

fitpar
m_dutH
m_angH
m_srH
m_dutV
m_angV
m_srV
m_ene

%--------------------------------------------------------------------------
% 3.) Save to file
%--------------------------------------------------------------------------
save(name,'horzsim','vertsim','z','m_dutH','m_angH','m_srH','m_dutV', ...
          'm_angV','m_srV','d_dutH','d_dutV','d_angH','d_angV','d_srH', ...
          'd_srV','fitpar','fitH','fitV','m_ene','d_ene','eneH','eneV');
end