classdef fullmatrix
    %inadmissible nodes
    properties
        cols
        rows
        e
    end
    
    methods
        function FM = fullmatrix(table)
            if nargin == 0
                FM.rows = 0;
                FM.cols = 0;
            else
                [FM.rows,FM.cols] = size(table);
                FM.e = table;
            end
        end
        
        
        %******************************************************************
        %****************************utilities*****************************
        %******************************************************************
        function table = getTable(FM)
            table = FM.e;
        end
        
        function [ret1, ret2] = getsize(FM)
            if nargout == 2
                ret1 = FM.rows;
                ret2 = FM.cols;
            else
                ret1 = [FM.rows FM.cols];
            end
        end
        
        function SM = getsupermatrix(FM)
            SM =  supermatrix();
            [SM.rows SM.cols] = FM.getsize();
            SM.f = FM;
        end
        
        %******************************************************************
        %****************************operators*****************************
        %******************************************************************
        function FM = plus(a,b)
            if isa(a, 'fullmatrix') && isa(b, 'fullmatrix')
                FM = fullmatrix(a.e + b.e);
            else
                FM = fullmatrix(a.getTable + b.getTable);
            end
        end
        
        function FM = minus(a,b)
            if isa(a, 'fullmatrix') && isa(b, 'fullmatrix')
                FM = fullmatrix(a.e - b.e);
            else
                FM = fullmatrix(a.getTable - b.getTable);
            end
        end
        
        function FM = uminus(FM)
            FM.e = -FM.e;
        end
        
        function FM = uplus(FM)
        end
        
        function FM = mtimes(a,b)
            if isa(a, 'fullmatrix') && isa(b, 'fullmatrix')
                FM = fullmatrix(a.e * b.e);
            end
        end
    end
    
end
