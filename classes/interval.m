classdef interval < handle
	% interval class: used for creating "interval" objects and contains
	% methods for "nesting" new intervals und checking that fit parameters
	% are within the given margins.
    properties
        DN       % lower margin for fitting parameter (is changed upon k++)
        UP       % upper margin for fitting parameter (is changed upon k++)
        DN_lim   % absolute lower limit for fitting parameter
        UP_lim   % absolute upper limit for fitting parameter
        del      % 2*del is the size of the current parted interval
    end
%--------------------------------------------------------------------------
% Methods
%-------------------------------------------------------------------------- 
    methods
    function obj = interval(dn,up)
        % sets the lower and upper limits when the object is created
        obj.DN_lim = dn;
        obj.UP_lim = up;
        obj.DN = dn;
        obj.UP = up;
    end
	function n_i = nestIntervals(obj,n_max,s)
        % nests n_max intervals with streching factor s (see first
        % subroutine "Nest intervals" from Fig. 3. for k>1, checklimits()
        % is run to make sure that the results stay within the given
        % margins
        if ~isempty(obj.del)
            obj.checklimits()
        end
        obj.del  = ( (obj.UP-obj.DN)/2 ) / n_max * s;
        n_i = linspace(obj.DN,obj.UP,2*n_max+1);
        n_i = n_i(2:2:end);
    end
	function setNewMinMax(obj,val)
        % sets new "min" and "max" margins according to the point "val" and
        % it's margin "obj.del"
        obj.DN = val-obj.del;
        obj.UP = val+obj.del;
    end
    function checklimits(obj,s)
        % if the lower or upper margins are exceeding the absolute margin,
        % then the parted interval is shifted towards the middle
        if obj.DN < obj.DN_lim
            corr_fac = obj.DN_lim-obj.DN;
            obj.DN = obj.DN_lim;
            obj.UP = obj.UP+corr_fac;
        end
        if obj.UP > obj.UP_lim
            corr_fac = obj.UP-obj.UP_lim;
            obj.UP = obj.UP_lim;
            obj.DN = obj.DN-corr_fac;
        end
    end
    end
end

