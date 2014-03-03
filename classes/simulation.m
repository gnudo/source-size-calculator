classdef simulation < handle
% contains methods for simulating a spherical wave propagation after a
% grating in 1D and 2D and analysis tools for loading experimental data
% and visibility processing for both experimental and simulation data
    properties (Constant)
        % natural constants
        c   = 299792458;              % [m/s]
        eV  = 1.6021764874e-19;       % [J] transfer to Joule
        h   = 6.6260689633e-34;       % [Js]
    end
    properties
        % constants from calling script
        a          % [m] grating period
        absorb     % [1] grating absorption (values between: 0-1)
        dc         % [1] duty cycle
        E          % [keV]
        alpha      % [°] angle of grating's pillar slopes
        h          % [m] height of grating's pillars
        N          % [1] number of particles --> 2^n !!!
        padding    % [1] total number (with zero-padding)
        periods    % [1] grating-size (in terms of periods)
        phShift    % [1] grating phase shift
        plotper    % [1] number of periods for plotting
        r          % [m] radius of incidient wave curvature
        srcsz      % [m] source size of beam for simulation
        rr         % [m] radius of incident wave curvature used for ML fit
        nameDest   % destination directory (full path)
        psize      % [m] px size of detector
        usewin     % switch whether to use window function
        z          % propagation distances
    end
    properties
        % calculated values from parameters
        D_def      % [m] defocussed D_T
        D_defr     % [m] defocussed D_R
        D_T        % [m] Talbot distance
        D_R        % [m] Replication distance
        du         % [1/m] distance in k-space
        dy0        % [m] distance between sampling points (from parameters)
        k          % wave vector
        lambda     % [m] wave length of impinging wave-front
        pxperiod   % [px] period of grating in px-values
        u          % [1/m] grating coordinates in k-space (zero-padded)
        x1         % [px] px-value for plotting-ROI
        x2         % [px] px-value for plotting-ROI
        y0         % [m] grating coordinates (zero-padded)
        per         % [px] number of px per period
        per_approx  % [px] approx. period due to beam divergence
        M           % magnification due to beam divergence
    end
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
methods
    function set.E(obj,value)
        % method that is called when energy E [keV] is "set". it
        % calculates wave length [m] of impinging wavefront, Talbot
        % distances [m] (defocussed, replication etc.) and the wavevector
        obj.E = value;      % set energy
        obj.lambda = obj.c.*obj.h ./ (obj.E.*obj.eV.*1000);
        obj.D_T    = 2.*(obj.a.^2)/obj.lambda;        % Talbot D
        obj.D_R    = obj.D_T/2;                       % replication D
        obj.D_def  = (obj.r*obj.D_T)/(obj.r-obj.D_T); % defocussed D_T
        obj.D_defr = (obj.r*obj.D_R)/(obj.r-obj.D_R); % defocussed D_R
        obj.k      = 2.*pi./obj.lambda;               % wave vector
    end
    function set.plotper(obj,value)
        % when setting "plotper" (number of periods to be plotted) the x1
        % and x2 coordinates are calculated to select the ROI from the
        % original img
        if rem(value,int32(value)) ~= 0
            error('number of period must be an INTEGER');
        end
        if isempty(obj.padding)
            obj.padding = obj.N;
        end
        obj.plotper = value;
        obj.pxperiod=obj.N/obj.periods;   % px per period
        middle = obj.padding/2;
        obj.x1 = middle - round( (value/2).*obj.pxperiod ) + 1;
        obj.x2 = middle + round( (value/2).*obj.pxperiod );
    end
    function calcRefracAbsorb(obj,density)
        % when setting "h" set correct phase-shift introduced by
        % grating and absorption level. here we make use of external
        % "xray-interaction-constants" library provided by Zhang Jiang
        %obj.h = value;
        result=refrac('Au',obj.E,density);
        delta = result.dispersion;
        bet   = result.absorption;
        obj.phShift = obj.k*delta*obj.h;
        obj.absorb  = (-bet.*obj.k.*obj.h);
    end
    function set.nameDest(obj,value)
        % set destination directory and mkdir
        mkdir(value);
        obj.nameDest = value;
    end
    function set.z(obj,value)
        % method that is called when propagation distances are set.
        % magnification M as well as the period in [px] is calculated
        obj.z   = value;
        obj.M   = ( obj.r + obj.z ) ./ obj.r;
        obj.per = obj.a/obj.psize;
        obj.per_approx = obj.per.*obj.M;
    end
    function G = talbotGrid1D (obj)
        % construction of 1D grid from all parameters. if "alpha" is
        % non-zero, then the grating bumps are assumed to be trapezoidal
        % and the grating is constructed.
        sca    = obj.N./obj.periods;
        xx     = round(sca.*obj.dc);
        ha_xx  = round(xx/2);
        if obj.alpha ~= 0 | isempty(obj.alpha) ~= 0
            angle = obj.alpha .* (pi./180);
            obj.dy0 = obj.periods.*(obj.a/obj.N); % TODO: remove redundancy with (waveFieldGrat)
            dmax   = ((sca-xx)*obj.dy0/2)/obj.h; % max length of bump's slope
            sinval = sin(angle);
            if sinval >= dmax
                error('angle of grating bumps is too big');
            end
            % length of slopes
            n1 = round(2*obj.h*sinval/obj.dy0);
            n2 = round(n1/2);
            s1 = linspace(1,0,n2);      % slope 1
            s2 = linspace(0,1,n1-n2);   % slope 2
            gra   = [ones(1,ha_xx) s1 zeros(1,sca-xx-n1) s2 ones(1,xx-ha_xx)];
        else
            gra   = [ones(1,ha_xx) zeros(1,sca-xx) ones(1,xx-ha_xx)];
        end

        % check whether total grating length will be the same as N,
        % otherwise correct for missing values
        if (length(gra).*obj.periods) == obj.N
            G  = repmat(gra,1,obj.periods);             % normal case
        elseif (length(gra).*obj.periods) < obj.N
            G = [];
            cornum = obj.N - length(gra).*obj.periods;
            for ii = 1:cornum
                G = [G gra 1];
            end
            for ii = 1:(obj.periods-cornum)
                G = [G gra];
            end
        elseif (length(gra).*obj.periods) > obj.N
            G = [];
            cornum = length(gra).*obj.periods - obj.N;
            for ii = 1:cornum
                G = [G gra(1:end-1)];
            end
            for ii = 1:(obj.periods-cornum)
                G = [G gra];
            end
        else
            disp('There was some failure in creating the grating.');
        end
        
        % ZEROPADDING grating (only if padding > N)mod
        if isempty(obj.padding)
            obj.padding = obj.N; % will become obsolete (substitute with N)
        end
        if (obj.padding ~= obj.N)
            gra0 = zeros(1,obj.padding);
            gra0(1,(obj.N/2+1):(3*obj.N/2))=G;
            G=gra0;
        end
    end
    function G = talbotGrid2D (obj)
        % construction of 2D grid by calling talbotGrid1D and creating
        % matrix
        grat1 = obj.talbotGrid1D;
        G = ones(length(grat1),length(grat1));
        for ii = 1:length(grat1)
            for jj = 1:length(grat1)
                G(ii,jj) = G(ii,jj) .* grat1(ii) .* grat1(jj);
            end
        end
    end
    function U = waveFieldGrat (obj, G)
        % calculates the wavefield of the impinging spherical wave after
        % passing the grating. (1) we consider the phase shift along the
        % grating (induced by the grating). (2) we consider the absorption
        % where the grating is. (3) we calculate the phase shift along the
        % grating due to the beam divergence.
        
        % Coordinate-specific values
        obj.dy0 = obj.periods.*(obj.a/obj.N);
        obj.y0  = [-(obj.padding/2):(obj.padding/2-1)].*obj.dy0;
        obj.du  = 1./ (obj.padding.*obj.dy0);
        obj.u   = [-(obj.padding/2):(obj.padding/2-1)].*obj.du;

        % Check whether grating is 1D or 2D
        if size(G,1) > 1 & size(G,2) > 1
            [XX YY] = meshgrid(obj.y0,obj.y0);
        else
            XX = obj.y0;
            YY = 0;
        end
        
        % Check whether curvature of impinging wavefront was differently
        % set
        if isempty(obj.rr)
            obj.rr = obj.r;
        end
        
        % Check whether "phase shift" and "absorption" level are set. Run
        % e.g. obj.calcRefracAbsorb(17) before obj.waveFieldGrat.
        if isempty(obj.phShift) | isempty(obj.absorb)
            error('obj.phShift and obj.absorb must be set before running obj.waveFieldGrat(G)!');
        end 
        
        % Grating Operators
        ds = sqrt(obj.rr.^2+XX.^2+YY.^2);    % beam diverg. phase shift
        U  = exp(i.*obj.k.*ds);              % (1) impinging spherical wave
        U  = U.*((exp(obj.absorb.*G)));      % (2) absorption
        U  = U.*exp(-i.*obj.phShift.*G);     % (3) phase-shift
        if (obj.padding ~= obj.N)
            intermed = zeros(1,obj.padding);
            intermed(1,(obj.N/2+1):(3*obj.N/2)) = U(1,(obj.N/2+1):(3*obj.N/2));
            U=intermed;
        end
    end
    function U = waveFieldProp (obj, z, U0)
        % calculates the itensity of a propagated wave field U0 along the
        % z-axis
        if isempty(U0)
            U0 = obj.waveFieldGrat(obj.talbotGrid1D);
        end
        ff = fftshift(fft(ifftshift(U0)));          % FFT of wave field
        H = exp( -i.*pi.*obj.lambda.*z.*obj.u.^2);  % Fresnel kernel
        C = exp(i.*obj.k.*z)./(2*obj.r);            % const
        C = C/abs(C);                               % normalized amplitude
        U = C .* fftshift(ifft(ifftshift(ff.*H)));  % convolution
        U = abs(U).^2;                              % intensity
    end
    function U = waveFieldPropMutual (obj, z, U0)
        % calculates the itensity of a propagated wave field U0 along the
        % z-axis and taking account of the mutual coherence (finite source
        % size) --> principle from Weitkamp-paper with Gaussian src
        w    = (z + (z==0)*1e-9) * (obj.srcsz/obj.r); % if z=0, set z=1e-9
        sigm_sq = w^2/(8*log(2));
        srcgauss = exp( -(1/2).* (obj.y0.^2)./ sigm_sq);
        srcgauss = srcgauss./sum(srcgauss);           % normalized Gauss
        
        U   = obj.waveFieldProp(z,U0);
        gam = fftshift(fft(ifftshift(srcgauss)));     % damping factor
        Uf  = fftshift(fft(ifftshift(U)));            % FFT of intensity
        U   = abs(fftshift(ifft(ifftshift(Uf.*gam))));
    end
    function U = waveFieldPropMutual2D (obj, z, U0)
        % calculates the itensity of a propagated wave field U0 along the
        % z-axis and taking account of the mutual coherence (given by
        % finite source size) --> principle from Weitkamp-paper with
        % Gaussian src --> all in 2D
        
        % Gaussian of Source sizes in both directions
        w    = (z + (z==0).*1e-9) * (obj.srcsz/obj.r); % if z=0, set z=1e-9
        sigm_sq = w.^2/(8.*log(2));
        
        [XX YY] = meshgrid(obj.y0,obj.y0);
        
        gaussX = exp(-(1/2).*(XX.^2)./sigm_sq(1));%./sqrt(2.*pi.*sigm_sq(1));
        gaussY = exp(-(1/2).*(YY.^2)./sigm_sq(2));%./sqrt(2.*pi.*sigm_sq(2));
        srcgauss = gaussX .* gaussY;
        srcgauss = srcgauss./sum(sum(srcgauss));  % normalized Gauss
        
        % if no wave-field at z=0 is given, it is calculated
        if isempty(U0)
            U0 = obj.waveFieldGrat(obj.talbotGrid2D);
        end

        [UU VV] = meshgrid(obj.u,obj.u);
        
        ff  = fftshift(fft2(ifftshift(U0)));         % FFT of wave field
        H   = exp( -i.*pi.*obj.lambda.*z.*(UU.^2+VV.^2)); % Fresnel kernel
        C   = exp(i.*obj.k.*z)./(2*obj.r);           % const
        C   = C/abs(C);                              % normalized amplitude
        U   = C .* fftshift(ifft2(ifftshift(ff.*H)));% convolution
        U   = abs(U).^2;                             % intensity
        gam = fftshift(fft2(ifftshift(srcgauss)));   % damping factor
        Uf  = fftshift(fft2(ifftshift(U)));          % FFT of intensity
        U   = abs(fftshift(ifft2(ifftshift(Uf.*gam))));
    end
    function F = calcFcoeff(obj)
        % calculates 1st Fourier coefficients in 1D from all class
        % properties

        % GRID creation & Wave field @ Grid
        obj.calcRefracAbsorb(17);
        gra = obj.talbotGrid1D;
        f = obj.waveFieldGrat(gra);
        
        % Propagation along z-axis
        obj.plotper = 14;      % number of periods to be plotted --> sets obj.x1, obj.x2
        resolu    = 512;     % desired resolution

        for ii=1:length(obj.z)
            lineTalbot = obj.waveFieldPropMutual(obj.z(ii),f); % Fresnel-pro + mutual coh.fnct
            crop = lineTalbot(obj.x1:obj.x2);                % crop to <plotper> periods
            crop = interp1(crop,linspace(1,length(crop),resolu));
            pwav(:,ii) = crop;
        end
        
        % Visibility calculation
        per = obj.M .* size(pwav,1)/obj.plotper;        % period in [px]

        for jj=1:length(obj.z)
            vec = pwav(:,jj);
            F(jj,1) = obj.vis1D (vec,per(jj));
        end 
    end
    function [F_H F_V] = calcFcoeff2D(obj,E,sigma,dc,alpha)
        % calculates 1st Fourier coefficients for all given parameters for
        % both the horizontal and vertical direction
        for ii=1:2
            obj.srcsz  = sigma(ii);
            obj.dc     = dc(ii);
            obj.alpha  = alpha(ii);
            obj.E      = E;
            F(:,ii)    = obj.calcFcoeff;
        end
        F_H = F(:,1);
        F_V = F(:,2);
    end
    function U = scale2Det (obj, U0)
        % scales Talbot-carpet to the pixel size (psize) of the detector.
        % if the images are not cropped, then obj.plotper is set with
        % obj.periods
        if isempty(obj.plotper)
            obj.plotper = obj.periods;
        end
        l     = obj.plotper.*obj.a;                  % [m] size of FOV
        res_o = length(obj.x1:obj.x2);               % [1] orig. resoultion
        res_n = round(obj.plotper.*obj.a./obj.psize);% [1] new resolution

        x_o   = linspace(0,l,res_o);      % original resolution from simu
        x_n   = linspace(0,l,res_n);      % resolution from detector
        
        if size(U0,1) > 1 & size(U0,2) > 1
            x = linspace(1,res_o,res_n);
            U = interp2(U0,x,x.');
        else
            U = interp1(x_o,U0,x_n);
        end
    end
    function k = primeFourierPos (obj, z)
        % calculates the expected 1st Fourier component according to the
        % resolution and number of periods saved and also taking account of
        % the magnification due to beam divergence
        M = (obj.r + z)./obj.r;
        k = (1./(obj.pxperiod.*M)).*length(obj.x1:obj.x2);
    end
    function f = weightedLSE (obj,simu,expo)
        % conducts the normalized weighted LSQ-error and returns the value
        f = ( (simu./mean(simu) - expo./mean(expo)).^2 ).*expo./mean(expo);
        f = sum(f);
    end
    function img     = loadSmallImg (obj,ind)
        % loads "small images" created with Save2img and throws failure if
        % folder is empty. img-s are assumed to be 16bit.
        % ind is the index of the img in the folder (not the name!)
        % TODO: probably obsolete
        files = dir([obj.nameDest '/*.tif']);
        if size(files,1) == 0
            error('No files loaded. Create small imgs first with Save2img')
        end
        file  = sprintf('%s/%s',obj.nameDest, files(ind).name)
        img   = double(imread(file,'tif'));
        img   = img./(2^16-1);
    end
	function Save2img (obj,img,n)
        % saves img to <obj.nameDest> under <000n.tif> at 16bit tif images
        zspace = '0000';
        zer=zspace(1:end-length(num2str(n)));
        file=sprintf('%s/%s%d.tif',obj.nameDest,zer,n)
        img = uint16((2^16-1).*(img));       % save as 16bit
        imwrite(img,file,'tif');
    end
    function [Fh Fv] = visCalc (obj,img,n)
        % TODO: - remove "set" dependency
        %       - change order of V and H (to be consistent with others)
        if ~exist('set','var')
            set = 1;        % will be obsolete...
        end
        % analyses visibility of <img> in vertical and horizontal direction
        per = obj.per_approx(set,n);  % approximated period in [px]
        mittel=mean(mean(img));

        % --- vertical
        meanv = mean(img,2);          % visibility in vertical   direction
        Fv = obj.vis1D(meanv,per)./mittel;

        % --- horizontal
        meanh = mean(img,1);          % visibility in horizontal direction
        Fh = obj.vis1D(meanh,per)./mittel;
    end
    function f = vis1D (obj,vec,per)
        % visibility calculation by identifying 1st Fourier component and
        % giving it's value.
        % TODO: include magnification (like in paper)
        dim = fn_leak1(vec,per);            % find correct window
        vec = vec(1:dim);                   % set
        if obj.usewin == 1
            vec=vec(:).*tukeywin(dim,0.75);
        end
        four = abs(fft(vec));
        first = round(length(vec)/per);     % expected position of 1st coef
        four_crp = four(2:round(dim/2));
        first = first-1;                    % because cropping
        [peak pos] = max(four_crp((first-6):(first+6)));
        f = peak/four(1);
    end
	function plotWaveAtGrating(obj,f,gra)
        % plot grating <gra> and grating-wave-front <f> with absorption and
        % phase part
        fig1=figure(1);
        set(fig1,'Position',[100 600 1024 400],'Color','white')
        area(obj.y0,gra);colormap summer;
        hold on;
        plot(obj.y0,angle(f),'ro');
        plot(obj.y0,abs(f).^2,'b*');
        line([0 obj.a],[1.02 1.02],'Color','g','LineWidth',3);
        legend('grating','phase shift','absorbtion')
        hold off
        title('grating')
    end
    function plotTalbotProfile(obj,calc)
        % plots talbot profile (should be called from the loop)
        fig2 = figure(2);
        set(fig2,'Position',[100 100 1024 400],'Color','white')
        sb1  = subplot (1,2,1);
        plot(calc);
        title('vertical mean')
    end
    function plotTalbotCarpet(obj,z,pwav)
        sb2 = subplot (1,2,2);
        figure(2)
        imagesc(z*1000,[],pwav);
        colormap(gray);
        line(obj.D_R.*1e3.*[1 1],[0 obj.padding],'Color','r','LineWidth',2);
        line(obj.D_T.*1e3.*[1 1],[0 obj.padding],'Color','r','LineWidth',2);
        line(obj.D_defr.*1e3.*[1 1],[0 obj.padding],'Color','b','LineWidth',2);
        line(obj.D_def.*1e3.*[1 1],[0 obj.padding],'Color','b','LineWidth',2);
        xlabel('propagation distance [mm]');
        legend('D_R','D_T','D_R (defocussed)','D_T (defocussed)');
    end
end
end
%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------
function n = fn_leak1(x,p)
% leakage function by R. Mokso
% x is the array, p is the approximate period

for k=1:1:round(p)
    t=x(1:end-k+1);
    ft=fft(t);
    aft=abs(ft);
    naft=aft/(length(t));
    pos=round(length(t)/p);
    [peak,peakpos]=max(naft(pos-2:pos+2));
    peakpos=pos-2+peakpos-1;
    f1(k)=peak;
end

[peak,peakpos]=max(f1);
n=length(x)-peakpos+1;
end