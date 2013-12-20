classdef simulation < handle
% contains methods for simulating a spherical wave propagation after a
% grating in 1D
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
        duty       % [1] duty cycle
        E          % [keV]
        gAngle     % [°] angle of grating's bump slopes
        gHeight    % [m] height of gratings bump
        N          % [1] number of particles --> 2^n !!!
        padding    % [1] total number (with zero-padding)
        periods    % [1] grating-size (in terms of periods)
        phShift    % [1] grating phase shift
        plotper    % [1] number of periods for plotting
        r          % [m] radius of incidient wave curvature
        srcsz      % [m] source size of beam for simulation
        rr         % [m] radius of incident wave curvature used for ML fit
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
    function set.gHeight(obj,value)
        % when setting "gHeight" set correct phase-shift introduced by
        % grating and absorption level. here we make use of external
        % "xray-interaction-constants" library provided by Zhang Jiang
        obj.gHeight = value;
        result=refrac('Au',obj.E,17);
        delta = result.dispersion;
        bet   = result.absorption;
        obj.phShift = obj.k*delta*obj.gHeight;
        obj.absorb  = (-bet.*obj.k.*obj.gHeight);
    end
    function G = talbotGrid1D (obj)
        % construction of 1D grid from all parameters. if "gAngle" is
        % non-zero, then the grating bumps are assumed to be trapezoidal
        % and the grating is constructed.
        sca    = obj.N./obj.periods;
        xx     = round(sca.*obj.duty);
        ha_xx  = round(xx/2);
        if obj.gAngle ~= 0 | isempty(obj.gAngle) ~= 0
            angle = obj.gAngle .* (pi./180);
            obj.dy0 = obj.periods.*(obj.a/obj.N); % TODO: remove redundancy with (waveFieldGrat)
            dmax   = ((sca-xx)*obj.dy0/2)/obj.gHeight; % max length of bump's slope
            sinval = sin(angle);
            if sinval >= dmax
                error('angle of grating bumps is too big');
            end
            % length of slopes
            n1 = round(2*obj.gHeight*sinval/obj.dy0);
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
        % passing the grating. (1) we calculate the phase shift along the
        % grating due to the beam divergence. (2) we consider the
        % absorption where the grating is. (3) we consider the phase shift
        % along the grating (induced by the grating)
        
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
        
        % Grating Operators
        ds = sqrt(obj.rr.^2+XX.^2+YY.^2)-obj.rr;   % beam diverg. phase shift
        U  = exp(i.*obj.k.*(ds + obj.rr));         % (1) impinging spherical wave
        U  = U.*((exp(obj.absorb.*G)));            % (2) absorption
        U  = U.*exp(-i.*obj.phShift.*G);           % (3) phase-shift
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
        U = abs(U).^2;
    end
    function U = waveFieldPropMutual (obj, z, U0)
        % calculates the itensity of a propagated wave field U0 along the
        % z-axis and taking account of the mutual coherence (given by
        % finite source size) --> principle from Weitkamp-paper with
        % Gaussian src
        w    = (z + (z==0)*1e-9) * (obj.srcsz/obj.r); % if z=0, set z=1e-9
        sigm_sq = w^2/(8*log(2));
        srcgauss = exp( -(1/2).* (obj.y0.^2)./ sigm_sq);
        srcgauss = srcgauss./sum(srcgauss);           % normalized Gauss
        
        % if no wave-field at z=0 is given, it is calculated
        if isempty(U0)
            U0 = obj.waveFieldGrat(obj.talbotGrid1D);
        end
        
        ff  = fftshift(fft(ifftshift(U0)));          % FFT of wave field
        H   = exp( -i.*pi.*obj.lambda.*z.*obj.u.^2); % Fresnel kernel
        C   = exp(i.*obj.k.*z)./(2*obj.rr);          % const
        C = C/abs(C);                                % normalized amplitude
        U   = C .* fftshift(ifft(ifftshift(ff.*H))); % convolution
        U   = abs(U).^2;                             % intensity
        gam = fftshift(fft(ifftshift(srcgauss)));    % damping factor
        Uf  = fftshift(fft(ifftshift(U)));           % FFT of intensity
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
        C = C/abs(C);                                % normalized amplitude
        U   = C .* fftshift(ifft2(ifftshift(ff.*H))); % convolution
        U   = abs(U).^2;                             % intensity
        gam = fftshift(fft2(ifftshift(srcgauss)));    % damping factor
        Uf  = fftshift(fft2(ifftshift(U)));           % FFT of intensity
        U   = abs(fftshift(ifft2(ifftshift(Uf.*gam))));
    end
    function U = scale2Det (obj, U0, psize)
        % scales Talbot-carpet to the pixel size (psize) of the detector.
        % if the images are not cropped, then obj.plotper is set with
        % obj.periods
        if isempty(obj.plotper)
            obj.plotper = obj.periods;
        end
        l     = obj.plotper.*obj.a;               % [m] size of FOV
        res_o = length(obj.x1:obj.x2);            % [1] original resoultion
        res_n = round(obj.plotper.*obj.a./psize); % [1] new resolution

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
    function [del bet]   = IndexRef (obj,file,meth)
        % fits NIST data for delta and beta to the given energy
        file    = sprintf('%s.dat',file);
        Evec    = obj.E;
        newdata = importdata(file);
        Enist   = newdata(:,1)./1000;  % [keV]
        delta   = newdata(:,2);
        beta    = newdata(:,3);
        meth   = lower(meth);
        [M,m,n]=unique(Enist);
        dupind = setdiff([1:size(Enist)],m); % position of edges
        if isempty(dupind) == 1
            vec1 = interp1(log(Enist),log(delta),log(Evec), ...
                    meth,'extrap');
            vec2 = interp1(log(Enist),log(beta),log(Evec), ...
                    meth,'extrap');
            l=length(Evec)+1;
            del = exp(vec1);
            bet = exp(vec2);
        else
            dE = 1e-12;
            Enist(dupind) = Enist(dupind)-dE;
            l=1;
            ll=1;
            del=[];
            bet=[];
            dupind= [dupind length(Enist)];
            for i=1:length(dupind)
                tmpE = Evec(Evec <= Enist(dupind(i)));
                tmpE = tmpE(l:end);
                vec1  = interp1(log(Enist(ll:(dupind(i)))), ...
                       log(delta(ll:(dupind(i)))),log(tmpE), meth,'extrap');
                vec2  = interp1(log(Enist(ll:(dupind(i)))), ...
                       log(beta(ll:(dupind(i)))),log(tmpE), meth,'extrap');
                del = [del exp(vec)];
                bet = [bet exp(vec)];
                l  = length(attnew)+1;
                ll = dupind(i)+1;
            end
        end
%         figure;
%         plot(log(Enist),log(delta));
%         hold on;
%         plot(log(Evec),log(del),'ro')
%         hold off;
        %flux   = F .* exp(-attnew.*thick.*rho); % multiplication with exp
    end
    function f = weightedLSQ (obj,simu,expo)
        % conducts the normalized weighted LSQ-error and returns the value
        f = ( (simu./mean(simu) - expo./mean(expo)).^2 ).*expo./mean(expo) ;
        f = sum(f);
    end
    function f = weightedLSQ2 (obj,simu,expo)
        % conducts the weighted LSQ-error and returns the value
        f = ( (simu - expo).^2 ).*expo;
        f = sum(f);
    end
    function y = LSQ(obj,a,b)
        % calculates LSQ normalized with mean value (non-used!!)
        y = sum( (a./mean(a) - b./mean(b)).^2 );
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