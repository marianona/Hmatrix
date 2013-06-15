classdef fullmatrix < matrixAbs
    %Representation of inadmissible nodes
    
    properties
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
        %***************************base methods***************************
        %******************************************************************
        function table = getTable(FM)
            table = FM.e;
        end
        
        function SM = getsupermatrix(FM)
            SM = supermatrix();
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
