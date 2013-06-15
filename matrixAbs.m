classdef matrixAbs
    %MATRIXABS An abstract base class for rk- full- and supermatrices
    
    properties
        rows
        cols
    end
    
    methods (Abstract)
        table = getTable(MA)
        SM = getsupermatrix(MA)
    end
    
    methods
        function [ret1, ret2] = getsize(MA)
            if nargout == 2
                ret1 = MA.rows;
                ret2 = MA.cols;
            else
                ret1 = [MA.rows MA.cols];
            end
        end
        
        function rank = getRank(MA)
            table = MA.getTable();
            if issparse(table)
                s = svd(full(table));
            else
                s = svd(table);
            end
            rank = nnz(s);
        end
    end
    
end
