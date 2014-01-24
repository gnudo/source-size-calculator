classdef interval < handle
% interval class: used for creating "interval" objects and contains
% methods for "nesting" new intervals und checking that fit parameters
% are within the given margins.
    properties
        min       % lower margin for fitting parameter (is changed upon k++)
        max       % upper margin for fitting parameter (is changed upon k++)
        min_lim   % absolute lower limit for fitting parameter
        max_lim   % absolute upper limit for fitting parameter
        del       % 2*del is the size of the current parted interval
        val       % current parameter value with min <= val <= max
    end
%--------------------------------------------------------------------------
% Methods
%-------------------------------------------------------------------------- 
    methods
    function obj = interval(dn,up)
        % sets the lower and upper limits when the object is created
        obj.min_lim = dn;
        obj.max_lim = up;
        obj.min = dn;
        obj.max = up;
    end
    function n_i = nestIntervals(obj,n_max,s)
        % nests n_max intervals with streching factor s (see first
        % subroutine "Nest intervals" from Fig. 3.
        % For k>1, checklimits() is run to make sure that the results stay
        % within the given (absolute) margins
        if ~isempty(obj.del)  % equivalent condition as: if k > 1
            obj.checklimits()
        end
        obj.del  = ( (obj.max-obj.min)/2 ) / n_max * s;
        n_i = linspace(obj.min,obj.max,2*n_max+1);
        n_i = n_i(2:2:end);
    end
    function setNewMinMax(obj,val)
        % stores the current parameter "val" and sets new "min" and "max"
        % margins according to the point "val" and it's margin "obj.del"
        obj.val = val;
        obj.min = val-obj.del;
        obj.max = val+obj.del;
    end
    function checklimits(obj,s)
        % if the lower or upper margins are exceeding the absolute margin,
        % then the parted interval is shifted towards the middle
        if obj.min < obj.min_lim
            corr_fac = obj.min_lim-obj.min;
            obj.min = obj.min_lim;
            obj.max = obj.max+corr_fac;
        end
        if obj.max > obj.max_lim
            corr_fac = obj.max-obj.max_lim;
            obj.max = obj.max_lim;
            obj.min = obj.min-corr_fac;
        end
    end
    end
end

