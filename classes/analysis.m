classdef analysis < handle
% contains all analysis tools for loading experimental data and visibility 
% processing for both experimental and simulation data
    properties
        % constants from calling script (TOMCAT secific things)
        nameDest    % destination directory (full path)
        z           % propagation distances (lines --> different data sets)
        usewin      % switch whether to use window function
        psize       % [m] px size of detector
        r           % [m] distance source <-> grid (exp. setup)
        a           % [m] grating period
    end
    properties
        % calculated values from parameters
        per         % [px] number of px per period
        per_approx  % [px] approx. period due to beam divergence
        M           % magnification due to beam divergence
    end
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
methods
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
    function [x1 x2] = middleROI(obj,reso,size)
        % places an small array (reso) in the middle of a big array (size)
        % and gives the start and ending coordinates
        space = size-reso;
        x1 = round(space/2);
        x2 = x1 + reso-1;
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
    function [Fv Fh] = visCalc (obj,img,n)
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
    function [Fv Fh] = visCalcF (obj,img,set,n)
        % analyses visibility of <img> with method in F-space along x and y
        % coordinates of the image
        per = obj.per_approx(set,n);  % approximated period in [px]
        mittel=mean(mean(img));
        
        % --- vertical
        for jj = 1:length(img)
            vec = img(:,jj);
            Fv(jj) = obj.vis1D(vec(:),per);
        end
        Fv = mean(Fv)./mittel;
                
        % --- horizontal
        for jj = 1:length(img)
            vec = img(jj,:);
             Fh(jj) = obj.vis1D(vec(:),per);
        end

        Fh = mean(Fh)./mittel;
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
    function Save2img (obj,img,n)
        % saves img to <obj.nameDest> under <000n.tif> at 16bit tif images
        zspace = '0000';
        zer=zspace(1:end-length(num2str(n)));
        file=sprintf('%s/%s%d.tif',obj.nameDest,zer,n)
        img = uint16((2^16-1).*(img));       % save as 16bit
        imwrite(img,file,'tif');
    end
end
%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------
end
