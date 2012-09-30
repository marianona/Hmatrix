classdef supermatrix
    properties
        rows
        cols
        block_rows
        block_cols
        r%rkmatrix - admissible
        f%fullmatrix - inadmissible
        s%tree nodes - supermatrixes
    end
    
    methods
        function SM = supermatrix(table, admissibility, domainx, domainy)
            SM.block_rows = 0;
            SM.block_cols = 0;
            if nargin == 0
                SM.rows = 0;
                SM.cols = 0;
            else
                [SM.rows, SM.cols] = size(table);
                if nargin == 1
                    SM.f = fullmatrix(table);
                elseif nargin == 4
                    if admissibility(domainx, domainy)
                        SM.r = rkmatrix(table);
                    else
                        SM.f = fullmatrix(table);
                    end
                end
            end
        end
        
        function SM = fulliterate(SM, admissibility, maxiters, minblocksize, err)
            SM = SM.fulliterate_(admissibility, maxiters, minblocksize, [0,1], [0,1], err);
        end
        
        function SM = fulliterate_(SM, admissibility, maxiters, minblocksize, domainx, domainy, err)
            if maxiters == -1 && minblocksize == -1
                error('Max iterations and Max size must  not be -1 at the same time\n');
            elseif SM.hasF ~= 1
                error('Supermatrix must be in the initial Fullmatrix form\n');
            elseif maxiters == 0 || (minblocksize ~= -1 && SM.rows * SM.cols <= minblocksize)
                if admissibility(domainx, domainy)
                    SM.r = rkmatrix(SM.f.e, err);
                    SM.f = [];
                end
            else
                if admissibility(domainx, domainy)
                    SM.r = rkmatrix(SM.f.e, err);
                    SM.f = [];
                else
                    newmaxiters = maxiters - 1;
                    if newmaxiters < 0
                        newmaxiters = -1;
                    end
                    r1 = ceil(SM.rows/2);
                    c1 = ceil(SM.cols/2);
                    tmpx = domainx(1)+(domainx(2)-domainx(1))*r1/SM.rows;
                    tmpy = domainy(1)+(domainy(2)-domainy(1))*c1/SM.cols;
                    SM.s = [
                        supermatrix(SM.f.e(1:r1, 1:c1)).fulliterate_(admissibility, newmaxiters, minblocksize, [domainx(1),tmpx], [domainy(1),tmpy], err) ...
                        supermatrix(SM.f.e(1:r1, c1+1:SM.cols)).fulliterate_(admissibility, newmaxiters, minblocksize, [domainx(1),tmpx], [tmpy,domainy(2)], err) ;
                        supermatrix(SM.f.e(r1+1:SM.rows, 1:c1)).fulliterate_(admissibility, newmaxiters, minblocksize, [tmpx,domainx(2)], [domainy(1),tmpy], err) ...
                        supermatrix(SM.f.e(r1+1:SM.rows, c1+1:SM.cols)).fulliterate_(admissibility, newmaxiters, minblocksize, [tmpx,domainx(2)], [tmpy,domainy(2)], err)
                        ];
                    SM.block_rows = 2;
                    SM.block_cols = 2;
                    SM.f = [];
                end
            end
        end
        
        function SM = iterate(SM, admissibility)
            SM = SM.iterate_(admissibility, [0,1], [0,1]);
        end
        
        function SM = iterate_(SM, admissibility, domainx, domainy)%eswteriki sinartisi Iterate me parametro to domain
            if isa(SM.r, 'rkmatrix')
                return;
            elseif isa(SM.f, 'fullmatrix')
                if admissibility(domainx, domainy)
                    SM.r = rkmatrix(SM.f.e);
                    SM.f = [];
                else
                    r1 = ceil(SM.rows/2);
                    c1 = ceil(SM.cols/2);
                    tmpx = domainx(1)+(domainx(2)-domainx(1))*r1/SM.rows;
                    tmpy = domainy(1)+(domainy(2)-domainy(1))*c1/SM.cols;
                    SM.s = [
                        supermatrix(SM.f.e(1:r1, 1:c1), admissibility, [domainx(1),tmpx], [domainy(1),tmpy]) ...
                        supermatrix(SM.f.e(1:r1, c1+1:SM.cols), admissibility, [domainx(1),tmpx], [tmpy,domainy(2)]) ;
                        supermatrix(SM.f.e(r1+1:SM.rows, 1:c1), admissibility, [tmpx,domainx(2)], [domainy(1),tmpy]) ...
                        supermatrix(SM.f.e(r1+1:SM.rows, c1+1:SM.cols), admissibility, [tmpx,domainx(2)], [tmpy,domainy(2)])
                        ];
                    SM.block_rows = 2;
                    SM.block_cols = 2;
                    SM.f = [];
                end
            elseif isa(SM.s, 'supermatrix')
                r1 = SM.s(1,1).rows;
                c1 = SM.s(1,1).cols;
                tmpx = domainx(1)+(domainx(2)-domainx(1))*r1/SM.rows;
                tmpy = domainy(1)+(domainy(2)-domainy(1))*c1/SM.cols;
                SM.s = [
                    SM.s(1,1).iterate_(admissibility, [domainx(1),tmpx], [domainy(1),tmpy]) ...
                    SM.s(1,2).iterate_(admissibility, [domainx(1),tmpx], [tmpy,domainy(2)]) ;
                    SM.s(2,1).iterate_(admissibility, [tmpx,domainx(2)], [domainy(1),tmpy]) ...
                    SM.s(2,2).iterate_(admissibility, [tmpx,domainx(2)], [tmpy,domainy(2)])
                    ];
            end
        end
        
        function SM = invert(SM)
            if SM.hasS()
                SM.s(1,1) = SM.s(1,1).invert();
                SM.s(1,2) = -(SM.s(1,1) * SM.s(1,2));
                SM.s(2,2) = muladd_supermatrix(SM.s(2,2), SM.s(2,1), SM.s(1,2));
                SM.s(2,1) = -(SM.s(2,1) * SM.s(1,1));
                
                SM.s(2,2) = SM.s(2,2).invert();
                SM.s(1,2) = SM.s(1,2) * SM.s(2,2);
                SM.s(1,1) = muladd_supermatrix(SM.s(1,1), SM.s(1,2), SM.s(2,1));
                SM.s(2,1) = SM.s(2,2) * SM.s(2,1);
            else
                SM = supermatrix(inv(SM.getTable()));
            end
        end
        
        %******************************************************************
        %****************************utilities*****************************
        %******************************************************************
        function table = getTable(SM)
            if isa(SM.r, 'rkmatrix')
                table = SM.r.getTable();
            elseif isa(SM.f, 'fullmatrix')
                table = SM.f.getTable();
            elseif isa(SM.s, 'supermatrix')
                table = [
                    SM.s(1,1).getTable()    SM.s(1,2).getTable()
                    SM.s(2,1).getTable()    SM.s(2,2).getTable()
                ];
            end
        end
        
        function [ret1, ret2] = getsize(SM)
            if nargout == 2
                ret1 = SM.rows;
                ret2 = SM.cols;
            else
                ret1 = [SM.rows SM.cols];
            end
        end
        
        function SM = getsupermatrix(SM)
        end
        
        function ret = hasR(SM)
            ret = isa(SM.r, 'rkmatrix');
        end
        
        function ret = hasF(SM)
            ret = isa(SM.f, 'fullmatrix');
        end
        
        function ret = hasS(SM)
            ret = isa(SM.s, 'supermatrix');
        end
        
        %******************************************************************
        %****************************operators*****************************
        %******************************************************************
        function SM = plus(a,b)
            if isa(a, 'supermatrix') && isa(b, 'supermatrix')
                if isa(a.r, 'rkmatrix') && isa(b.r, 'rkmatrix')
                    sum = a.r + b.r;
                    SM = sum.getsupermatrix;
                elseif isa(a.f, 'fullmatrix') && isa(b.f, 'fullmatrix')
                    sum = a.f + b.f;
                    SM = sum.getsupermatrix;
                elseif isa(a.s, 'supermatrix') && isa(b.s, 'supermatrix') && a.block_rows == b.block_rows && a.block_cols == b.block_cols
                    SM = supermatrix();
                    [SM.rows SM.cols] = a.getsize();
                    SM.block_rows = a.block_rows;
                    SM.block_cols = a.block_cols;
                    SM.s = supermatrix;
                    SM.s(SM.block_rows, SM.block_cols) = supermatrix;
                    for i = 1:numel(a.s)
                        SM.s(i) = a.s(i) + b.s(i);
                    end
                else
                    SM = supermatrix(a.getTable() + b.getTable());
                end
            end
        end

        function SM = minus(a,b)
            SM = a + (-b);
        end
        
        function SM = uminus(SM)
            if isa(SM.r, 'rkmatrix')
                SM.r = -SM.r;
            elseif isa(SM.f, 'fullmatrix')
                SM.f = -SM.f;
            elseif isa(SM.s, 'supermatrix')
                for i = 1:numel(SM.s)
                    SM.s(i) = -SM.s(i);
                end
            end
        end
        
        function SM = uplus(SM)
        end
        
        function SM = mtimes(a,b)
            if isa(a, 'supermatrix') && isa(b, 'supermatrix')
                SM = mul_supermatrix(a, b);
            end
        end
    end
    
end
