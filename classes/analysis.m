classdef analysis < handle
% contains all analysis tools for loading experimental data and visibility 
% processing for both experimental and simulation data
    properties (Constant)
        % constants
        psize   = 0.38e-6;  % [m] px size of detector
        r       = 26.3;     % [m] distance source <-> grid (exp. setup)
        a       = 6.84e-6;  % [m] grating period
        range   = 1.1;      % [m] range of detector movement
    end
    properties
        % constants from calling script (TOMCAT secific things)
        nameData    % name of data
        nameDest    % destination directory (full path)
        ndarks      % number of dark img-s in beginning
        noffset     % offset to 1st img
        nexp        % factor of multiple imgages for longer exposure time
        z           % propagation distances (lines --> different data sets)
        usewin      % switch whether to use window function
    end
    properties
        % calculated values from parameters
        src         % source path
        dark        % dark image
        per         % [px] number of px per period
        per_approx  % [px] approx. period due to beam divergence
        M           % magnification due to beam divergence
    end
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
methods
    function set.nameData(obj,value)
        % when setting "nameData" the file "dataInfo.csv" is loaded where
        % the first line represents the data names with corresponding
        % columns containing the explicit path. if an incorrect name is
        % given, the while loop doesn't throw a return and the error is
        % given. (TOMCAT specific)
        obj.nameData = value;
        data  = fn_read_mixed_csv('dataInfo.csv',',');
        k = 1;
        while k <= size(data,2)
            if strcmp(value, data{1,k}) == 1
                obj.src{1} = data{2,k};
                if length(data{3,k}) > 0
                    obj.src{2} = data{3,k};
                end
                return
            end
            k = k+1;
        end
        error('NO such datafile! check "dataInfo.csv" or data-src')
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
    function img     = loadDarks (obj)
        % loads dark img-s from given parameters. if file set is divided
        % into 2 subfolders, then img is 3D-matrix. (TOMCAT specific)
        for ii = 1:length(obj.src)      % loop through subsets of data
            files = dir([obj.src{ii} '/*.tif']);
            for kk = 1:obj.ndarks       % loop through number of darks
                file = sprintf('%s/%s',obj.src{ii}, files(kk).name)
                darks(:,:,kk)=double(imread(file));
            end
        img(:,:,ii) = double(median(darks,3));
        end
        obj.dark = img;
    end
    function [x1 x2] = middleROI(obj,reso,size)
        % places an small array (reso) in the middle of a big array (size)
        % and gives the start and ending coordinates
        space = size-reso;
        x1 = round(space/2);
        x2 = x1 + reso-1;
    end
	function img     = loadAndCorrect (obj,set,n)
        % loads imgs and (dark+flat field corrected). "even": img, "odd":
        % flat (TOMCAT specific)
        files = dir([obj.src{set} '/*.tif']);
        ind   = obj.noffset + 2*n;
        file1 = sprintf('%s/%s',obj.src{set}, files(ind-1).name)
        file2 = sprintf('%s/%s',obj.src{set}, files(ind).name)
        
        img_tmp1 = double(imread(file1,'tif'));      % flat
        img_tmp2 = double(imread(file2,'tif'));      % image
        img = (img_tmp2-obj.dark(:,:,set)) ./ (img_tmp1-obj.dark(:,:,set));
    end
    function img     = loadMultiple (obj,set,n)
        % load multiple img-s w/o flat-/darkfield-correction
        % TODO: probably obsolete
        files = dir([obj.src{set} '/*.tif']);
        strt  = (obj.noffset+1+(n-1)*(obj.nexp+1));
        ende  = strt+obj.nexp;
        file  = sprintf('%s/%s',obj.src{set}, files(strt).name)
        img = 0;
        for ind = (strt+1):ende
            file     = sprintf('%s/%s',obj.src{set}, files(ind).name)
            img = img + double(imread(file,'tif'))-obj.dark(:,:,set);      % image
        end
        %img = img./mean((mean(img)));
    end
    function img     = loadMultipleCorr (obj,set,n)
        % load multiple img-s with flat-/darkfield-correction
        % TODO: probably obsolete
        files = dir([obj.src{set} '/*.tif']);
        strt  = (obj.noffset+1+(n-1)*(obj.nexp+1));
        ende  = strt+obj.nexp;
        %file  = sprintf('%s/%s',obj.src{set}, files(strt).name)
        img = 0;
        for ind = (strt+1):ende
            file     = sprintf('%s/%s',obj.src{set}, files(ind).name)
            flat     = sprintf('%s/%s',obj.src{set}, files(ende+1).name);
            img = img + (double(imread(file,'tif'))-obj.dark(:,:,set))./(double(imread(flat,'tif'))-obj.dark(:,:,set));      % image
        end
        %img = img./mean((mean(img)));
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
    function [Fv Fh] = visCalc (obj,img,set,n)
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
        % saves img to <dest> under <000n.tif>
        % TODO: probably obsolete
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
