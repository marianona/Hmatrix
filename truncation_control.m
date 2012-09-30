classdef truncation_control
    properties
        rel_eps
        abs_eps
        %adaptive
    end
    
    methods
        function TC = truncation_control()
            TC.rel_eps = eps;
            TC.abs_eps = eps;
            %TC.adaptive = 1;
        end
    end
    
end
