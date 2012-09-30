function [numOfnnzRK, numOfRK, total, xx, yy] = calc_rkeigenvalues_rec( SM )
    numOfnnzRK = 0;    
    numOfRK = 0;
    total = 0;
    xx = [];
    yy = [];
    if SM.hasR
        ss = svd(full(SM.getTable))';
        if ss(1) ~= 0
            numOfnnzRK = 1;
            ss = ss / ss(1);
            
            lenss = length(ss);
            xx = (1:lenss)/lenss;
            yy = ss;
        end
        
        numOfRK = 1;
        total = 1;
    elseif SM.hasS
        for i = 1 : numel(SM.s)
            [nnzrk, rk, tot, xxt, yyt] = calc_rkeigenvalues_rec(SM.s(i));
            numOfnnzRK = numOfnnzRK + nnzrk;
            numOfRK = numOfRK + rk;
            total = total + tot;
            xx = [xx xxt]; %#ok<AGROW>
            yy = [yy yyt]; %#ok<AGROW>
        end
    else
        total = 1;
    end
end

