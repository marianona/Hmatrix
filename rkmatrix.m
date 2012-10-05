classdef rkmatrix
    %admissible nodes
    properties
        k
        kt
        tc
        rows
        cols
        a
        b
    end
    
    methods
        function RK = rkmatrix(table, err)
            if nargin == 0
                RK.k = 0;
                RK.kt = 0;
                RK.rows = 0;
                RK.cols = 0;
                RK.a = [];
                RK.b = [];
            else
                [RK.rows, RK.cols] = size(table);
                
                if issparse(table)
                    [RK.a s v] = svd(full(table));
                else
                    [RK.a s v] = svd(table);
                end
                RK.k = nnz(diag(s));
                % svd's s size is rows x cols
                % if needed reduce the dimensions of a and b
                if  RK.rows > RK.k || RK.cols > RK.k
                    RK.a = RK.a(:, 1:RK.k);
                    RK.b = v(:, 1:RK.k) * s(1:RK.k, 1:RK.k)';
                else
                    RK.b = v * s';
                end
                
                RK.kt = RK.k;
                % low rank approximation
                if nargin == 2 && err > 0
                    for i = (RK.k-1):-1:0
                        tmpa = RK.a(:, 1:i);
                        tmpb = RK.b(:, 1:i);
                        r = tmpa * tmpb';
                        if issparse(table)
                            n = norm(full(r - table)) / norm(full(table));
                        else
                            n = norm(r - table) / norm(table);
                        end
                        if n > err
                            if i + 1 < RK.kt
                                RK.kt = i + 1;
                                RK.a = RK.a(:, 1:RK.kt);
                                RK.b = RK.b(:, 1:RK.kt);
                            end
                            break;
                        end
                    end
                end
                
            end
        end
        
        function [u s v] = rsvd(RK)
            [u_a, v_a] = qr(RK.a);
            [u_b, v_b] = qr(RK.b);
            usv = v_a(:, 1:RK.kt) * v_b(:, 1:RK.kt)';
            
            if issparse(usv)
                [u_s, s, v_s] =  svd(full(usv));
            else
                [u_s, s, v_s] =  svd(usv);
            end
            snnz = nnz(diag(s));
            % svd's s size is rows x cols
            % if needed reduce the dimensions of u_s and v_s
            if  RK.rows > snnz || RK.cols > snnz
                s = s(1:snnz, 1:snnz);
                u = u_a * u_s(:, 1:snnz);
                v = u_b * v_s(:, 1:snnz);
            else
                u = u_a * u_s;
                v = u_b * v_s;
            end
        end
        
        %******************************************************************
        %****************************utilities*****************************
        %******************************************************************
        function table = getTable(RK)
            table = RK.a * RK.b';
        end
        
        function [ret1, ret2] = getsize(RK)
            if nargout == 2
                ret1 = RK.rows;
                ret2 = RK.cols;
            else
                ret1 = [RK.rows RK.cols];
            end
        end
        
        function SM = getsupermatrix(RK)
            SM =  supermatrix();
            [SM.rows SM.cols] = RK.getsize();
            SM.r = RK;
        end
        
        %******************************************************************
        %****************************operators*****************************
        %******************************************************************
        function RK = plus(a,b)
            if isa(a, 'rkmatrix') && isa(b, 'rkmatrix')
                RK = rkmatrix();
                [RK.rows RK.cols] = a.getsize();
                RK.tc = truncation_control();
                RK = add_rkmatrix(RK, a, b);
            end
        end

        function RK = minus(a,b)
            if isa(a, 'rkmatrix') && isa(b, 'rkmatrix')
                RK = a + (-b);
            end
        end
        
        function RK = uminus(RK)
            RK.a = -RK.a;
        end
        
        function RK = uplus(RK)
        end
        
        function RK = mtimes(a,b)
            if isa(a, 'rkmatrix') && isa(b, 'rkmatrix')
                %if ((a.rows - b.cols) * a.kt * b.kt <= (a.kt - b.kt) * a.rows * b.cols)
                % assuming a.rows almost equal to b.cols
                if (b.kt < a.kt)
                    RK = rkmatrix((a.a * (a.b' * b.a)) * b.b');
                else
                    RK = rkmatrix(a.a * ((a.b' * b.a) * b.b'));
                end
            end
        end
    end
    
end
