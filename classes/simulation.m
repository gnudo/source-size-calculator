classdef simulation < handle
% contains methods for simulating a spherical wave propagation after a
% grating in 1D and 2D and analysis tools for loading experimental data
% and visibility processing for both experimental and simulation data
    properties (Constant)
        % natural constants
        c        = 299792458;              % [m/s]
        eV       = 1.6021764874e-19;       % [J] transfer to Joule
        h_planck = 6.6260689633e-34;       % [Js]
    end
    properties
        % variables with same nomenclature as in the paper
        a          % [m] grating period
        alpha      % [°] angle of grating's pillar slopes
        dc         % [1] duty cycle
        E          % [keV]
        h          % [m] height of grating's pillars
        k          % [1/m ]wave vector
        lambda     % [m] wave length of impinging wave-front
        M          % magnification due to beam divergence
        psize      % [m] px size of detector
        R          % [m] radius of incidient wave curvature
        RR         % [m] radius of incident wave curvature used for ML fit
        sigma      % [m] source size of beam for simulation
        z          % [m] array of all propagation distances
    end
    properties
        % variables from numerical implementation
        absorb     % [1] grating absorption (values between: 0-1)
        N          % [1] number of particles --> 2^n !!!
        nameDest   % destination directory (full path)
        periods    % [1] grating-size (in terms of periods)
        phShift    % [1] grating phase shift
        usewin     % switch whether to use window function
    end
    properties
        % calculated values from parameters
        D_def      % [m] defocussed D_T
        D_defr     % [m] defocussed D_R
        D_T        % [m] Talbot distance
        D_R        % [m] Replication distance
        du         % [1/m] distance in k-space
        dy0        % [m] distance between sampling points (from parameters)
        pxperiod   % [px] period of grating in px-values
        u          % [1/m] grating coordinates in k-space (zero-padded)
        x1         % [px] px-value for plotting-ROI
        x2         % [px] px-value for plotting-ROI
        y0         % [m] grating coordinates (zero-padded)
        per        % [px] number of px per period
        per_approx % [px] approx. period due to beam divergence
    end
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
methods
    function set.E(obj,value)
        % method that is called when energy E [keV] is "set". it
        % calculates wave length [m] of impinging wavefront, Talbot
        % distances [m] (defocussed, replication etc.) and the wavevector
        obj.E      = value;                           % set energy
        obj.lambda = obj.c.*obj.h_planck ./ (obj.E.*obj.eV.*1000);
        obj.D_T    = 2.*(obj.a.^2)/obj.lambda;        % Talbot D
        obj.D_R    = obj.D_T/2;                       % replication D
        obj.D_def  = (obj.R*obj.D_T)/(obj.R-obj.D_T); % defocussed D_T
        obj.D_defr = (obj.R*obj.D_R)/(obj.R-obj.D_R); % defocussed D_R
        obj.k      = 2.*pi./obj.lambda;               % wave vector
    end
    function set.nameDest(obj,value)
        % set destination directory and mkdir if folder doesn't exist
        if exist(value,'dir') ~= 7
            mkdir(value);
        end
        obj.nameDest = value;
    end
    function set.z(obj,value)
        % method that is called when propagation distances are set.
        % magnification M as well as the period in [px] is calculated
        obj.z   = value;
        
        % Check whether curvature of impinging wavefront was differently
        % set
        if isempty(obj.RR)
            obj.RR = obj.R;
        end
        obj.M   = ( obj.RR + obj.z ) ./ obj.RR;
        obj.per = obj.a/obj.psize;
        obj.per_approx = obj.per.*obj.M;
    end
    function calcRefracAbsorb(obj,material,density)
        % when setting "h" set correct phase-shift introduced by
        % grating and absorption level. here we make use of external
        % "xray-interaction-constants" library provided by Zhang Jiang
        result = refrac(material,obj.E,density);
        delta  = result.dispersion;
        bet    = result.absorption;
        obj.phShift = -obj.k*delta*obj.h;
        obj.absorb  = -bet.*obj.k.*obj.h;
    end
    function G = talbotGrid1D (obj)
        % construction of 1D grid from all parameters. if "alpha" is
        % non-zero, then the grating bumps are assumed to be trapezoidal
        % and the grating is constructed.
        sca    = round(obj.N./obj.periods);
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
    function f = waveFieldGrat (obj, G)
        % calculates the wavefield of the impinging spherical wave after
        % passing the grating.
        
        % Coordinate-specific values
        obj.dy0 = obj.periods.*(obj.a/obj.N);
        obj.y0  = [-(obj.N/2):(obj.N/2-1)].*obj.dy0;
        obj.du  = 1./ (obj.N.*obj.dy0);
        obj.u   = [-(obj.N/2):(obj.N/2-1)].*obj.du;

        % Check whether grating is 1D or 2D
        if size(G,1) > 1 & size(G,2) > 1
            [XX YY] = meshgrid(obj.y0,obj.y0);
        else
            XX = obj.y0;
            YY = 0;
        end
        
        % Check whether "phase shift" and "absorption" level are set. Run
        % e.g. obj.calcRefracAbsorb(17) before obj.waveFieldGrat.
        if isempty(obj.phShift) | isempty(obj.absorb)
            error('obj.phShift and obj.absorb must be set before running obj.waveFieldGrat(G)!');
        end 
        
        % Eq. (11) from paper
        f  = exp( i.*obj.phShift.*G ) .*exp( obj.absorb.*G ) .* ...
                              exp( i.*obj.k.*sqrt(obj.RR.^2+XX.^2+YY.^2) );
    end
    function psi = waveFieldProp (obj, z, f)
        % calculates the itensity of a propagated wave field f along the
        % z-axis. "psi" is the absolute square root of Eq. (8).
        if isempty(f)
            f = obj.waveFieldGrat(obj.talbotGrid1D);
        end
        ff  = fftshift(fft(ifftshift(f)));           % FFT of wave field
        H   = exp( -i.*pi.*obj.lambda.*z.*obj.u.^2); % Fresnel kernel
        C   = exp(i.*obj.k.*z)./(2*obj.R);           % const
        C   = C/abs(C);                              % normalized amplitude
        psi = C .* fftshift(ifft(ifftshift(ff.*H))); % convolution
        psi = abs(psi).^2;                           % intensity
    end
    function I = waveFieldPropMutual (obj, z, f)
        % calculates the itensity of a propagated wave field f along the
        % z-axis and taking account of the mutual coherence (principle from
        % Weitkamp-paper with Gaussian src). implements Eqs. (9) and (10).
        w    = (z + (z==0)*1e-9) * (obj.sigma/obj.R); % if z=0, set z=1e-9
        sigm_sq = w^2/(8*log(2));
        srcgauss = exp( -(1/2).* (obj.y0.^2)./ sigm_sq);
        srcgauss = srcgauss./sum(srcgauss);           % normalized Gauss
        
        I   = obj.waveFieldProp(z,f);
        gam = fftshift(fft(ifftshift(srcgauss)));     % damping factor
        Uf  = fftshift(fft(ifftshift(I)));            % FFT of intensity
        I   = abs(fftshift(ifft(ifftshift(Uf.*gam))));
    end
    function I = waveFieldPropMutual2D (obj, z, f)
        % calculates the itensity of a propagated wave field f along the
        % z-axis and taking account of the mutual coherence (principle from
        % Weitkamp-paper with Gaussian src) --> all in 2D
        
        % Gaussian of Source sizes in both directions
        w    = (z + (z==0).*1e-9) * (obj.sigma/obj.R); % if z=0, set z=1e-9
        sigm_sq = w.^2/(8.*log(2));
        
        [XX YY] = meshgrid(obj.y0,obj.y0);
        
        gaussX = exp(-(1/2).*(XX.^2)./sigm_sq(1));
        gaussY = exp(-(1/2).*(YY.^2)./sigm_sq(2));
        srcgauss = gaussX .* gaussY;
        srcgauss = srcgauss./sum(sum(srcgauss));       % normalized Gauss
        
        % if no wave-field at z=0 is given, it is calculated
        if isempty(f)
            f = obj.waveFieldGrat(obj.talbotGrid2D);
        end

        [UU VV] = meshgrid(obj.u,obj.u);
        
        ff  = fftshift(fft2(ifftshift(f)));         % FFT of wave field
        H   = exp(-i.*pi.*obj.lambda.*z.*(UU.^2+VV.^2)); % Fresnel kernel
        C   = exp(i.*obj.k.*z)./(2*obj.R);           % const
        C   = C/abs(C);                              % normalized amplitude
        I   = C .* fftshift(ifft2(ifftshift(ff.*H)));% convolution
        I   = abs(I).^2;                             % intensity
        gam = fftshift(fft2(ifftshift(srcgauss)));   % damping factor
        If  = fftshift(fft2(ifftshift(I)));          % FFT of intensity
        I   = abs(fftshift(ifft2(ifftshift(If.*gam))));
    end
    function F = calculateFsim(obj)
        % calculates 1st Fourier coefficients in 1D from all class
        % properties

        % GRID creation & Wave field @ Grid
        gra = obj.talbotGrid1D;
        f = obj.waveFieldGrat(gra);
        
        % Correction for border areas. We reduce the number of simulated
        % periods by 10%
        obj.pxperiod=obj.N/obj.periods;           % px per period
        periods_crop = round(0.9 .* obj.periods);
        middle = obj.N/2;
        x1 = middle - round( (periods_crop/2).*obj.pxperiod ) + 1;
        x2 = middle + round( (periods_crop/2).*obj.pxperiod );

        % Propagation along z-axis
        resolu    = 512;     % desired resolution
        for ii=1:length(obj.z)
            lineTalbot = obj.waveFieldPropMutual(obj.z(ii),f); % Fresnel-pro + mutual coh.fnct
            crop = lineTalbot(x1:x2);                % crop to <plotper> periods
            crop = interp1(crop,linspace(1,length(crop),resolu));
            pwav(:,ii) = crop;
        end
        
        % Visibility calculation
        per = obj.M .* size(pwav,1)/periods_crop; % period in [px]

        for jj=1:length(obj.z)
            vec = pwav(:,jj);
            F(jj,1) = obj.vis1D (vec,per(jj));
        end 
    end
    function [F_H F_V] = calculateFsim2D(obj,E,sigma,dc,alpha)
        % calculates 1st Fourier coefficients for all given parameters for
        % both the horizontal and vertical direction
        for ii=1:2
            obj.sigma  = sigma(ii);
            obj.dc     = dc(ii);
            obj.alpha  = alpha(ii);
            obj.E      = E;
            F(:,ii)    = obj.calculateFsim;
        end
        F_H = F(:,1);
        F_V = F(:,2);
    end
    function f = scale2Det (obj, f)
        % scales Talbot-carpet to the pixel size (psize) of the detector.
        l     = obj.periods.*obj.a;       % [m] size of FOV
        res_o = length(f);                % [1] orig. resoultion
        res_n = round(l./obj.psize);      % [1] new resolution

        x_o   = linspace(0,l,res_o);      % original resolution from simu
        x_n   = linspace(0,l,res_n);      % resolution from detector
        
        if size(f,1) > 1 & size(f,2) > 1
            x = linspace(1,res_o,res_n);
            f = interp2(f,x,x.');
        else
            f = interp1(x_o,f,x_n);
        end
    end
    function p = weightedLSE (obj,simu,expo)
        % conducts the normalized weighted LSQ-error and returns the value
        p = ( (simu./mean(simu) - expo./mean(expo)).^2 ).*expo./mean(expo);
        p = sum(p);
    end
    function img = loadSmallImg (obj,ind)
        % loads "small images" created with Save2img and throws failure if
        % folder is empty. img-s are assumed to be 16bit.
        % ind is the index of the img in the folder (not the name!)
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
        img = uint16((2^16-1).*(img)); % save as 16bit
        imwrite(img,file,'tif');
    end
    function [Fh Fv] = visCalc (obj,img,n)
        % analyses visibility of <img> in vertical and horizontal direction
        per = obj.per_approx(n);       % approximated period in [px]
        img_mean =mean(img(:));

        % --- horizontal
        meanh = mean(img,1);           % visibility in horizontal direction
        Fh = obj.vis1D(meanh,per)./img_mean;

        % --- vertical
        meanv = mean(img,2);           % visibility in vertical   direction
        Fv = obj.vis1D(meanv,per)./img_mean;
    end
    function f = vis1D (obj,vec,per)
        % visibility calculation by identifying 1st Fourier component and
        % giving it's value.
        dim = fn_leak1(vec,per);            % find correct window
        vec = vec(1:dim);                   % set
        if obj.usewin == 1
            vec=vec(:).*tukeywin(dim,0.75);
        end
        four = abs(fft(vec));
        kx = round(length(vec)/per);     % expected position of 1st coef
        four_crp = four(2:round(dim/2));
        kx = kx-1;                    % because cropping
        [peak pos] = max(four_crp((kx-4):(kx+4)));
        f = peak/four(1);
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